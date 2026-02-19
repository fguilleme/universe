#include "presets.h"

#include <ctype.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>

#include "textures.h"

static uint32_t xorshift32(uint32_t *state) {
  uint32_t x = *state;
  x ^= x << 13;
  x ^= x >> 17;
  x ^= x << 5;
  *state = x;
  return x;
}

static double frand01(uint32_t *state) {
  return (double)(xorshift32(state) & 0xFFFFFF) / (double)0x1000000;
}

static double au_from_km(double km) { return km / 149597870.7; }

static double au_per_year_from_km_per_s(double km_s) {
  return km_s * 31557600.0 / 149597870.7;
}

static double visual_radius_from_physical(double r_au, double extra_scale) {
  // Exaggerate sizes while keeping relative ordering.
  // Use a low exponent to keep smaller bodies visible without making the Sun
  // unrealistically large compared to orbital distances.
  return pow(fmax(r_au, 1e-12), 0.30) * 0.60 * extra_scale;
}

static bool file_exists(const char *path) {
  struct stat st;
  return path && path[0] && stat(path, &st) == 0 && S_ISREG(st.st_mode);
}

static const char *skip_ws(const char *s) {
  while (s && *s && isspace((unsigned char)*s))
    s++;
  return s;
}

static void rstrip(char *s) {
  if (!s)
    return;
  size_t n = strlen(s);
  while (n > 0 && isspace((unsigned char)s[n - 1]))
    s[--n] = 0;
}

static void strip_quotes_inplace(char *s) {
  if (!s)
    return;
  const size_t n = strlen(s);
  if (n >= 2 && ((s[0] == '"' && s[n - 1] == '"') ||
                 (s[0] == '\'' && s[n - 1] == '\''))) {
    memmove(s, s + 1, n - 2);
    s[n - 2] = 0;
  }
}

typedef struct {
  char name[32];
  double mass_solar;
  double radius_km;
  uint32_t layer;
  double tilt_deg;
  double rotation_period_days;
  double x_km, y_km, z_km;
  double vx_km_s, vy_km_s, vz_km_s;
} BodyInit;

typedef enum {
  PRESET_YAML_KIND_UNKNOWN = 0,
  PRESET_YAML_KIND_BODIES,
  PRESET_YAML_KIND_TWO_BODY,
  PRESET_YAML_KIND_DISK_GALAXY,
} PresetYamlKind;

typedef struct {
  bool found;
  PresetYamlKind kind;

  bool has_sim_G;
  double sim_G;
  bool has_sim_softening;
  double sim_softening;
  bool has_sim_density;
  double sim_density;
  bool has_sim_merge_on_collision;
  bool sim_merge_on_collision;

  // Generator params.
  bool has_primary_mass;
  double primary_mass;
  bool has_secondary_mass;
  double secondary_mass;
  bool has_orbit_radius;
  double orbit_radius;

  bool has_n;
  int n;
  bool has_disk_radius;
  double disk_radius;
  bool has_central_mass;
  double central_mass;
  bool has_seed;
  uint32_t seed;
  bool has_noise;
  double noise;

  bool has_primary_radius;
  double primary_radius;
  bool has_secondary_radius;
  double secondary_radius;
  bool has_center_radius;
  double center_radius;

  bool has_time_scale;
  double time_scale;

  bool has_cam_cx;
  double cam_cx;
  bool has_cam_cy;
  double cam_cy;
  bool has_cam_cz;
  double cam_cz;
  bool has_cam_yaw;
  double cam_yaw;
  bool has_cam_pitch;
  double cam_pitch;
  bool has_cam_zoom_world_h;
  double cam_zoom_world_h;

  bool has_selected_name;
  char selected_name[32];

  bool has_follow_selected;
  bool follow_selected;

  // Bodies list.
  BodyInit *bodies;
  size_t body_count;
  size_t body_cap;
} PresetYaml;

static void add_body_init_list(Sim *sim, const BodyInit *bodies, size_t count);

static void preset_yaml_destroy(PresetYaml *p) {
  if (!p)
    return;
  free(p->bodies);
  *p = (PresetYaml){0};
}

static uint32_t parse_tex_layer(const char *s) {
  if (!s)
    return TEX_GENERIC;

  while (*s && isspace((unsigned char)*s))
    s++;
  if (!*s)
    return TEX_GENERIC;

  // Allow numeric layers.
  if ((*s >= '0' && *s <= '9') || *s == '-') {
    char *end = NULL;
    const long v = strtol(s, &end, 10);
    if (end && end != s && v >= 0)
      return (uint32_t)v;
  }

  // Allow symbolic names used in code.
  if (strcmp(s, "TEX_SUN") == 0)
    return TEX_SUN;
  if (strcmp(s, "TEX_MERCURY") == 0)
    return TEX_MERCURY;
  if (strcmp(s, "TEX_VENUS") == 0)
    return TEX_VENUS;
  if (strcmp(s, "TEX_EARTH") == 0)
    return TEX_EARTH;
  if (strcmp(s, "TEX_MARS") == 0)
    return TEX_MARS;
  if (strcmp(s, "TEX_JUPITER") == 0)
    return TEX_JUPITER;
  if (strcmp(s, "TEX_SATURN") == 0)
    return TEX_SATURN;
  if (strcmp(s, "TEX_URANUS") == 0)
    return TEX_URANUS;
  if (strcmp(s, "TEX_NEPTUNE") == 0)
    return TEX_NEPTUNE;
  return TEX_GENERIC;
}

static bool parse_double(const char *s, double *out) {
  if (!s || !out)
    return false;
  char *end = NULL;
  const double v = strtod(s, &end);
  if (!end || end == s)
    return false;
  *out = v;
  return true;
}

static bool parse_bool(const char *s, bool *out) {
  if (!s || !out)
    return false;
  if (strcmp(s, "true") == 0 || strcmp(s, "True") == 0 || strcmp(s, "1") == 0 ||
      strcmp(s, "yes") == 0 || strcmp(s, "on") == 0) {
    *out = true;
    return true;
  }
  if (strcmp(s, "false") == 0 || strcmp(s, "False") == 0 ||
      strcmp(s, "0") == 0 || strcmp(s, "no") == 0 || strcmp(s, "off") == 0) {
    *out = false;
    return true;
  }
  return false;
}

static bool parse_int(const char *s, int *out) {
  if (!s || !out)
    return false;
  char *end = NULL;
  const long v = strtol(s, &end, 10);
  if (!end || end == s)
    return false;
  *out = (int)v;
  return true;
}

static bool parse_u32(const char *s, uint32_t *out) {
  if (!s || !out)
    return false;
  char *end = NULL;
  const unsigned long v = strtoul(s, &end, 0);
  if (!end || end == s)
    return false;
  *out = (uint32_t)v;
  return true;
}

static int count_indent(const char *s) {
  int n = 0;
  while (s && s[n] == ' ')
    n++;
  return n;
}

static PresetYamlKind parse_preset_kind(const char *s) {
  if (!s || !s[0])
    return PRESET_YAML_KIND_UNKNOWN;
  if (strcmp(s, "bodies") == 0)
    return PRESET_YAML_KIND_BODIES;
  if (strcmp(s, "two_body") == 0)
    return PRESET_YAML_KIND_TWO_BODY;
  if (strcmp(s, "disk_galaxy") == 0)
    return PRESET_YAML_KIND_DISK_GALAXY;
  return PRESET_YAML_KIND_UNKNOWN;
}

static bool preset_yaml_push_body(PresetYaml *py, const BodyInit *b) {
  if (!py || !b)
    return false;
  if (py->body_count == py->body_cap) {
    const size_t next = py->body_cap ? py->body_cap * 2 : 16;
    BodyInit *nb = (BodyInit *)realloc(py->bodies, next * sizeof(BodyInit));
    if (!nb)
      return false;
    py->bodies = nb;
    py->body_cap = next;
  }
  py->bodies[py->body_count++] = *b;
  return true;
}

static const char *resolve_presets_yaml_path(void) {
  const char *env = getenv("UNIVERSE_PRESETS_YAML");
  if (env && env[0] && file_exists(env))
    return env;
  if (file_exists("data/presets.yaml"))
    return "data/presets.yaml";
  if (file_exists("../data/presets.yaml"))
    return "../data/presets.yaml";
  return NULL;
}

static bool load_preset_yaml(const char *preset_name, PresetYaml *out,
                            const char **out_used_path) {
  if (!preset_name || !out)
    return false;
  *out = (PresetYaml){0};
  if (out_used_path)
    *out_used_path = NULL;

  const char *path = resolve_presets_yaml_path();
  if (!path)
    return false;

  FILE *f = fopen(path, "rb");
  if (!f)
    return false;

  bool in_presets = false;
  bool in_target = false;
  bool in_bodies = false;

  bool have_cur = false;
  BodyInit cur = {0};

  char line[1024];
  while (fgets(line, (int)sizeof(line), f)) {
    // Strip comments.
    char *hash = strchr(line, '#');
    if (hash)
      *hash = 0;
    rstrip(line);

    const int indent = count_indent(line);
    const char *t = skip_ws(line);
    if (!t || !*t)
      continue;

    if (!in_presets) {
      if (indent == 0 && strcmp(t, "presets:") == 0)
        in_presets = true;
      continue;
    }

    // New preset header.
    if (indent == 2) {
      const size_t tl = strlen(t);
      if (tl >= 2 && t[tl - 1] == ':') {
        char name[64];
        const size_t nl = tl - 1;
        if (nl < sizeof(name)) {
          memcpy(name, t, nl);
          name[nl] = 0;
          const bool match = (strcmp(name, preset_name) == 0);

          if (in_target && in_bodies && have_cur && cur.name[0]) {
            (void)preset_yaml_push_body(out, &cur);
          }
          in_target = match;
          in_bodies = false;
          have_cur = false;
          cur = (BodyInit){0};
          if (match) {
            out->found = true;
            if (out_used_path)
              *out_used_path = path;
          }
        }
        continue;
      }
    }

    if (!in_target)
      continue;

    // Leaving bodies block.
    if (in_bodies && indent <= 4) {
      if (have_cur && cur.name[0]) {
        (void)preset_yaml_push_body(out, &cur);
      }
      in_bodies = false;
      have_cur = false;
      cur = (BodyInit){0};
      // Keep parsing this line in case it's another preset key.
    }

    if (!in_bodies && indent == 4 && strcmp(t, "bodies:") == 0) {
      in_bodies = true;
      continue;
    }

    // Start of a new body item.
    if (in_bodies && indent == 6 && t[0] == '-' &&
        (t[1] == 0 || isspace((unsigned char)t[1]))) {
      if (have_cur && cur.name[0]) {
        (void)preset_yaml_push_body(out, &cur);
      }
      cur = (BodyInit){0};
      have_cur = true;

      const char *p = skip_ws(t + 1);
      if (!p || !*p)
        continue;

      // Allow "- name: ..." on the same line.
      t = p;
      // fallthrough to key/value parsing
    }

    const bool can_parse_kv =
        (!in_bodies && indent == 4) || (in_bodies && indent >= 8);
    if (!can_parse_kv)
      continue;

    const char *colon = strchr(t, ':');
    if (!colon)
      continue;

    char key[64];
    const size_t klen = (size_t)(colon - t);
    if (klen == 0 || klen >= sizeof(key))
      continue;
    memcpy(key, t, klen);
    key[klen] = 0;
    rstrip(key);

    char val[256];
    const char *v = skip_ws(colon + 1);
    strncpy(val, v ? v : "", sizeof(val) - 1);
    val[sizeof(val) - 1] = 0;
    rstrip(val);
    strip_quotes_inplace(val);

    if (!in_bodies) {
      if (strcmp(key, "kind") == 0) {
        out->kind = parse_preset_kind(val);
      } else if (strcmp(key, "sim_G") == 0) {
        out->has_sim_G = parse_double(val, &out->sim_G);
      } else if (strcmp(key, "sim_softening") == 0) {
        out->has_sim_softening = parse_double(val, &out->sim_softening);
      } else if (strcmp(key, "sim_density") == 0) {
        out->has_sim_density = parse_double(val, &out->sim_density);
      } else if (strcmp(key, "sim_merge_on_collision") == 0) {
        out->has_sim_merge_on_collision =
            parse_bool(val, &out->sim_merge_on_collision);
      } else if (strcmp(key, "primary_mass") == 0) {
        out->has_primary_mass = parse_double(val, &out->primary_mass);
      } else if (strcmp(key, "secondary_mass") == 0) {
        out->has_secondary_mass = parse_double(val, &out->secondary_mass);
      } else if (strcmp(key, "orbit_radius") == 0) {
        out->has_orbit_radius = parse_double(val, &out->orbit_radius);
      } else if (strcmp(key, "primary_radius") == 0) {
        out->has_primary_radius = parse_double(val, &out->primary_radius);
      } else if (strcmp(key, "secondary_radius") == 0) {
        out->has_secondary_radius = parse_double(val, &out->secondary_radius);
      } else if (strcmp(key, "center_radius") == 0) {
        out->has_center_radius = parse_double(val, &out->center_radius);
      } else if (strcmp(key, "n") == 0) {
        out->has_n = parse_int(val, &out->n);
      } else if (strcmp(key, "disk_radius") == 0) {
        out->has_disk_radius = parse_double(val, &out->disk_radius);
      } else if (strcmp(key, "central_mass") == 0) {
        out->has_central_mass = parse_double(val, &out->central_mass);
      } else if (strcmp(key, "seed") == 0) {
        out->has_seed = parse_u32(val, &out->seed);
      } else if (strcmp(key, "noise") == 0) {
        out->has_noise = parse_double(val, &out->noise);
      } else if (strcmp(key, "time_scale") == 0) {
        out->has_time_scale = parse_double(val, &out->time_scale);
      } else if (strcmp(key, "cam_cx") == 0) {
        out->has_cam_cx = parse_double(val, &out->cam_cx);
      } else if (strcmp(key, "cam_cy") == 0) {
        out->has_cam_cy = parse_double(val, &out->cam_cy);
      } else if (strcmp(key, "cam_cz") == 0) {
        out->has_cam_cz = parse_double(val, &out->cam_cz);
      } else if (strcmp(key, "cam_yaw") == 0) {
        out->has_cam_yaw = parse_double(val, &out->cam_yaw);
      } else if (strcmp(key, "cam_pitch") == 0) {
        out->has_cam_pitch = parse_double(val, &out->cam_pitch);
      } else if (strcmp(key, "cam_zoom_world_h") == 0) {
        out->has_cam_zoom_world_h = parse_double(val, &out->cam_zoom_world_h);
      } else if (strcmp(key, "selected_name") == 0) {
        out->has_selected_name = true;
        strncpy(out->selected_name, val, sizeof(out->selected_name) - 1);
        out->selected_name[sizeof(out->selected_name) - 1] = 0;
      } else if (strcmp(key, "follow_selected") == 0) {
        out->has_follow_selected = parse_bool(val, &out->follow_selected);
      }
    } else {
      if (strcmp(key, "name") == 0) {
        strncpy(cur.name, val, sizeof(cur.name) - 1);
        cur.name[sizeof(cur.name) - 1] = 0;
      } else if (strcmp(key, "mass_solar") == 0) {
        (void)parse_double(val, &cur.mass_solar);
      } else if (strcmp(key, "radius_km") == 0) {
        (void)parse_double(val, &cur.radius_km);
      } else if (strcmp(key, "tex_layer") == 0) {
        cur.layer = parse_tex_layer(val);
      } else if (strcmp(key, "tilt_deg") == 0) {
        (void)parse_double(val, &cur.tilt_deg);
      } else if (strcmp(key, "rotation_period_days") == 0) {
        (void)parse_double(val, &cur.rotation_period_days);
      } else if (strcmp(key, "x_km") == 0) {
        (void)parse_double(val, &cur.x_km);
      } else if (strcmp(key, "y_km") == 0) {
        (void)parse_double(val, &cur.y_km);
      } else if (strcmp(key, "z_km") == 0) {
        (void)parse_double(val, &cur.z_km);
      } else if (strcmp(key, "vx_km_s") == 0) {
        (void)parse_double(val, &cur.vx_km_s);
      } else if (strcmp(key, "vy_km_s") == 0) {
        (void)parse_double(val, &cur.vy_km_s);
      } else if (strcmp(key, "vz_km_s") == 0) {
        (void)parse_double(val, &cur.vz_km_s);
      }
    }
  }

  if (in_target && in_bodies && have_cur && cur.name[0]) {
    (void)preset_yaml_push_body(out, &cur);
  }

  fclose(f);
  return out->found;
}

static void preset_apply_sim_overrides(Sim *sim, const PresetYaml *py) {
  if (!sim || !py)
    return;
  if (py->has_sim_G)
    sim->G = py->sim_G;
  if (py->has_sim_softening)
    sim->softening = py->sim_softening;
  if (py->has_sim_density)
    sim->density = py->sim_density;
  if (py->has_sim_merge_on_collision)
    sim->merge_on_collision = py->sim_merge_on_collision;
}

static uint32_t preset_select_by_name(const Sim *sim, const char *name) {
  if (!sim || !name || !name[0])
    return 0;
  for (size_t i = 0; i < sim->count; i++) {
    if (strcmp(sim->bodies[i].name, name) == 0)
      return sim->bodies[i].id;
  }
  return 0;
}

bool preset_apply_from_yaml(const char *preset_name, Sim *sim, Camera *cam,
                            double *time_scale, uint32_t *selected_id,
                            bool *follow_selected) {
  if (!preset_name || !sim)
    return false;

  PresetYaml py;
  const char *used_path = NULL;
  if (!load_preset_yaml(preset_name, &py, &used_path))
    return false;

  sim_reset(sim);
  // Reset sim parameters to known defaults (sim_reset() does not).
  sim->G = 1.0;
  sim->softening = 0.01;
  sim->density = 1.0;
  sim->merge_on_collision = true;
  preset_apply_sim_overrides(sim, &py);

  if (py.kind == PRESET_YAML_KIND_BODIES) {
    if (py.body_count == 0) {
      preset_yaml_destroy(&py);
      return false;
    }
    add_body_init_list(sim, py.bodies, py.body_count);
  } else if (py.kind == PRESET_YAML_KIND_TWO_BODY) {
    if (!py.has_primary_mass || !py.has_secondary_mass || !py.has_orbit_radius ||
        !py.has_primary_radius || !py.has_secondary_radius) {
      preset_yaml_destroy(&py);
      return false;
    }
    const double r = py.orbit_radius;
    const double r_safe = fmax(r, 1e-9);
    const double v = sqrt(sim->G * py.primary_mass / r_safe);

    (void)sim_add_body(sim, (Body){.name = "Primary",
                                   .x = 0.0,
                                   .y = 0.0,
                                   .vx = 0.0,
                                   .vy = 0.0,
                                   .mass = py.primary_mass,
                                   .radius = py.primary_radius,
                                   .render_radius = py.primary_radius,
                                   .r = 1.0f,
                                   .g = 0.9f,
                                   .b = 0.6f});
    (void)sim_add_body(sim, (Body){.name = "Secondary",
                                   .x = r,
                                   .y = 0.0,
                                   .vx = 0.0,
                                   .vy = v,
                                   .mass = py.secondary_mass,
                                   .radius = py.secondary_radius,
                                   .render_radius = py.secondary_radius,
                                   .r = 0.5f,
                                   .g = 0.8f,
                                   .b = 1.0f});
  } else if (py.kind == PRESET_YAML_KIND_DISK_GALAXY) {
    if (!py.has_n || !py.has_disk_radius || !py.has_central_mass ||
        !py.has_center_radius || !py.has_seed || !py.has_noise) {
      preset_yaml_destroy(&py);
      return false;
    }

    (void)sim_add_body(sim, (Body){.name = "Center",
                                   .x = 0.0,
                                   .y = 0.0,
                                   .vx = 0.0,
                                   .vy = 0.0,
                                   .mass = py.central_mass,
                                   .radius = py.center_radius,
                                   .render_radius = py.center_radius,
                                   .r = 1.0f,
                                   .g = 0.85f,
                                   .b = 0.5f});

    uint32_t rng = py.seed;
    for (int i = 0; i < py.n; i++) {
      const double u = frand01(&rng);
      const double vv = frand01(&rng);
      const double rr = py.disk_radius * sqrt(u);
      const double a = vv * 2.0 * M_PI;
      const double x = rr * cos(a);
      const double y = rr * sin(a);

      const double speed = sqrt(sim->G * py.central_mass / fmax(rr, 1e-3));
      const double tx = -sin(a);
      const double ty = cos(a);
      const double noise = (frand01(&rng) - 0.5) * py.noise;

      Body b = {0};
      b.x = x;
      b.y = y;
      b.vx = tx * speed * (1.0 + noise);
      b.vy = ty * speed * (1.0 + noise);
      b.mass = 1.0 + 15.0 * frand01(&rng);
      b.r = (float)(0.4 + 0.6 * frand01(&rng));
      b.g = (float)(0.4 + 0.6 * frand01(&rng));
      b.b = 1.0f;
      (void)sim_add_body(sim, b);
    }
  } else {
    preset_yaml_destroy(&py);
    return false;
  }

  if (cam) {
    if (py.has_cam_cx)
      cam->cx = py.cam_cx;
    if (py.has_cam_cy)
      cam->cy = py.cam_cy;
    if (py.has_cam_cz)
      cam->cz = py.cam_cz;
    if (py.has_cam_yaw)
      cam->yaw = py.cam_yaw;
    if (py.has_cam_pitch)
      cam->pitch = py.cam_pitch;
    if (py.has_cam_zoom_world_h)
      cam->zoom_world_h = py.cam_zoom_world_h;
  }

  if (time_scale && py.has_time_scale)
    *time_scale = py.time_scale;

  if (follow_selected && py.has_follow_selected)
    *follow_selected = py.follow_selected;

  if (selected_id) {
    if (py.has_selected_name) {
      *selected_id = preset_select_by_name(sim, py.selected_name);
      if (!*selected_id && follow_selected)
        *follow_selected = false;
    } else {
      *selected_id = 0;
    }
  }

  fprintf(stderr, "Loaded preset %s from %s\n", preset_name,
          used_path ? used_path : "(unknown)");
  preset_yaml_destroy(&py);
  return true;
}

static void add_body_init_list(Sim *sim, const BodyInit *bodies, size_t count) {
  for (size_t i = 0; i < count; i++) {
    const BodyInit *in = &bodies[i];
    Body b = {0};
    strncpy(b.name, in->name, sizeof(b.name) - 1);
    b.name[sizeof(b.name) - 1] = 0;

    b.x = au_from_km(in->x_km);
    b.y = au_from_km(in->y_km);
    b.z = au_from_km(in->z_km);
    b.vx = au_per_year_from_km_per_s(in->vx_km_s);
    b.vy = au_per_year_from_km_per_s(in->vy_km_s);
    b.vz = au_per_year_from_km_per_s(in->vz_km_s);
    b.mass = in->mass_solar;

    const double r_au = au_from_km(in->radius_km);
    b.radius = r_au;
    b.render_radius = visual_radius_from_physical(r_au, 1.0);
    b.tex_layer = in->layer;
    b.tilt_rad = (float)(in->tilt_deg * (M_PI / 180.0));
    {
      const double pd = in->rotation_period_days;
      if (pd != 0.0) {
        const double sign = (pd >= 0.0) ? 1.0 : -1.0;
        const double pd_abs = fabs(pd);
        b.spin_rate_rad_per_time =
            (float)(sign * (2.0 * M_PI * 365.25) / pd_abs);
      }
    }
    b.spin_phase_rad = 0.0f;
    b.r = 1.0f;
    b.g = 1.0f;
    b.b = 1.0f;
    (void)sim_add_body(sim, b);
  }
}

void preset_seed_two_body_orbit(Sim *sim) {
  if (preset_apply_from_yaml("two_body", sim, NULL, NULL, NULL, NULL))
    return;

  // Fallback if YAML is missing.
  sim_reset(sim);
  sim->G = 1.0;
  sim->softening = 0.02;
  sim->density = 2.0;
  sim->merge_on_collision = true;

  double sun_mass = 1.0e6;
  double planet_mass = 10.0;
  double r = 30.0;

  PresetYaml py;
  const char *used_path = NULL;
  if (load_preset_yaml("two_body", &py, &used_path)) {
    if (py.has_sim_G)
      sim->G = py.sim_G;
    if (py.has_sim_softening)
      sim->softening = py.sim_softening;
    if (py.has_sim_density)
      sim->density = py.sim_density;
    if (py.has_sim_merge_on_collision)
      sim->merge_on_collision = py.sim_merge_on_collision;

    if (py.kind == PRESET_YAML_KIND_BODIES && py.body_count) {
      fprintf(stderr, "Loaded preset two_body (bodies) from %s\n", used_path);
      add_body_init_list(sim, py.bodies, py.body_count);
      preset_yaml_destroy(&py);
      return;
    }

    if (py.has_primary_mass)
      sun_mass = py.primary_mass;
    if (py.has_secondary_mass)
      planet_mass = py.secondary_mass;
    if (py.has_orbit_radius)
      r = py.orbit_radius;
    fprintf(stderr, "Loaded preset two_body from %s\n", used_path);
    preset_yaml_destroy(&py);
  }

  const double r_safe = fmax(r, 1e-9);
  const double v = sqrt(sim->G * sun_mass / r_safe);

  (void)sim_add_body(sim, (Body){.name = "Primary",
                                 .x = 0.0,
                                 .y = 0.0,
                                 .vx = 0.0,
                                 .vy = 0.0,
                                 .mass = sun_mass,
                                 .radius = 1.0,
                                 .render_radius = 1.0,
                                 .r = 1.0f,
                                 .g = 0.9f,
                                 .b = 0.6f});
  (void)sim_add_body(sim, (Body){.name = "Secondary",
                                 .x = r,
                                 .y = 0.0,
                                 .vx = 0.0,
                                 .vy = v,
                                 .mass = planet_mass,
                                 .radius = 0.5,
                                 .render_radius = 0.5,
                                 .r = 0.5f,
                                 .g = 0.8f,
                                 .b = 1.0f});
}

void preset_seed_disk_galaxy(Sim *sim, int n, double disk_radius,
                             double central_mass) {
  if (preset_apply_from_yaml("disk_galaxy", sim, NULL, NULL, NULL, NULL))
    return;

  // Fallback if YAML is missing.
  sim_reset(sim);
  sim->G = 1.0;
  sim->softening = 0.05;
  sim->density = 1.5;
  sim->merge_on_collision = true;

  int use_n = n;
  double use_disk_radius = disk_radius;
  double use_central_mass = central_mass;
  uint32_t seed = 0x12345678u;
  double noise_scale = 0.08;

  PresetYaml py;
  const char *used_path = NULL;
  if (load_preset_yaml("disk_galaxy", &py, &used_path)) {
    if (py.has_sim_G)
      sim->G = py.sim_G;
    if (py.has_sim_softening)
      sim->softening = py.sim_softening;
    if (py.has_sim_density)
      sim->density = py.sim_density;
    if (py.has_sim_merge_on_collision)
      sim->merge_on_collision = py.sim_merge_on_collision;

    if (py.kind == PRESET_YAML_KIND_BODIES && py.body_count) {
      fprintf(stderr, "Loaded preset disk_galaxy (bodies) from %s\n",
              used_path);
      add_body_init_list(sim, py.bodies, py.body_count);
      preset_yaml_destroy(&py);
      return;
    }

    if (py.has_n)
      use_n = py.n;
    if (py.has_disk_radius)
      use_disk_radius = py.disk_radius;
    if (py.has_central_mass)
      use_central_mass = py.central_mass;
    if (py.has_seed)
      seed = py.seed;
    if (py.has_noise)
      noise_scale = py.noise;
    fprintf(stderr, "Loaded preset disk_galaxy from %s\n", used_path);
    preset_yaml_destroy(&py);
  }

  (void)sim_add_body(sim, (Body){.name = "Center",
                                 .x = 0.0,
                                 .y = 0.0,
                                 .vx = 0.0,
                                 .vy = 0.0,
                                 .mass = use_central_mass,
                                 .radius = 0.5,
                                 .render_radius = 1.0,
                                 .r = 1.0f,
                                 .g = 0.85f,
                                 .b = 0.5f});

  uint32_t rng = seed;
  for (int i = 0; i < use_n; i++) {
    const double u = frand01(&rng);
    const double v = frand01(&rng);
    const double r = use_disk_radius * sqrt(u);
    const double a = v * 2.0 * M_PI;
    const double x = r * cos(a);
    const double y = r * sin(a);

    // Circular-ish orbit around center with slight noise.
    const double speed = sqrt(sim->G * use_central_mass / fmax(r, 1e-3));
    const double tx = -sin(a);
    const double ty = cos(a);
    const double noise = (frand01(&rng) - 0.5) * noise_scale;

    Body b = {0};
    b.x = x;
    b.y = y;
    b.vx = tx * speed * (1.0 + noise);
    b.vy = ty * speed * (1.0 + noise);
    b.mass = 1.0 + 15.0 * frand01(&rng);
    b.r = (float)(0.4 + 0.6 * frand01(&rng));
    b.g = (float)(0.4 + 0.6 * frand01(&rng));
    b.b = 1.0f;
    (void)sim_add_body(sim, b);
  }
}

void preset_seed_solar_system(Sim *sim) {
  if (preset_apply_from_yaml("solar_system", sim, NULL, NULL, NULL, NULL))
    return;

  // Fallback if YAML is missing.
  // Accurate-ish 3D initial conditions from JPL Horizons (DE441),
  // barycentric ecliptic-of-J2000 frame at 2000-01-01 00:00 TDB.
  // Units in this simulation:
  // - distance: AU
  // - mass: solar mass
  // - time: years
  // In these units, G = 4*pi^2.
  sim_reset(sim);
  sim->G = 4.0 * M_PI * M_PI;
  sim->softening = 1e-6;
  sim->density = 1.0;
  sim->merge_on_collision = false;

  PresetYaml py;
  const char *used_path = NULL;
  if (load_preset_yaml("solar_system", &py, &used_path)) {
    if (py.has_sim_G)
      sim->G = py.sim_G;
    if (py.has_sim_softening)
      sim->softening = py.sim_softening;
    if (py.has_sim_density)
      sim->density = py.sim_density;
    if (py.has_sim_merge_on_collision)
      sim->merge_on_collision = py.sim_merge_on_collision;

    if (py.kind == PRESET_YAML_KIND_BODIES && py.body_count) {
      fprintf(stderr, "Loaded preset solar_system from %s (%zu bodies)\n",
              used_path, py.body_count);
      add_body_init_list(sim, py.bodies, py.body_count);
      preset_yaml_destroy(&py);
      return;
    }
    preset_yaml_destroy(&py);
  }

  fprintf(stderr,
          "Presets YAML not found/invalid; using built-in solar_system preset. "
          "(Set UNIVERSE_PRESETS_YAML=...)\n");

  const BodyInit built_in[] = {
      {"Sun", 1.0, 695700.0, TEX_SUN, 7.25, 25.38, -1.068108951496322E+06,
       -4.177210908491462E+05, 3.086887010002915E+04, 9.305302656256911E-03,
       -1.283177282717393E-02, -1.631700118015769E-04},
      {"Mercury", 1.660e-7, 2439.7, TEX_MERCURY, 0.034, 58.646,
       -2.212073002393702E+07, -6.682435921338345E+07, -3.461577076477692E+06,
       3.666229234452722E+01, -1.230266984222893E+01, -4.368336206255391E+00},
      {"Venus", 2.447e-6, 6051.8, TEX_VENUS, 177.36, -243.025,
       -1.075068040813442E+08, -1.053381041265540E+07, 6.024630234834461E+06,
       2.949586892816910E+00, -3.520915385418408E+01, -5.945882137071832E-01},
      {"Earth", 3.003e-6, 6371.0, TEX_EARTH, 23.439, 0.99726968,
       -2.521092392536846E+07, 1.449177150228508E+08, -2.394650324223201E+04,
       -2.983862964912880E+01, -5.218619426624726E+00, 1.193210383927141E-03},
      {"Mars", 3.227e-7, 3389.5, TEX_MARS, 25.19, 1.02595675,
       -2.155838897834890E+08, -3.539753056080124E+07, 4.532312700768938E+06,
       3.886195690590875E+00, -2.315847231020980E+01, -5.777378581341476E-01},
      {"Jupiter", 9.545e-4, 69911.0, TEX_JUPITER, 3.13, 0.41354,
       5.990175449652451E+08, 4.394132557364160E+08, -1.520962039552424E+07,
       -7.882322868602207E+00, 1.134831390779115E+01, 1.246425392362529E-01},
      {"Saturn", 2.858e-4, 58232.0, TEX_SATURN, 26.73, 0.44401,
       9.587267807903839E+08, -9.825068848392826E+08, -1.695064654623866E+07,
       5.908189281781105E+00, 6.296823371892003E+00, -3.704231008630100E-01},
      {"Uranus", 4.366e-5, 25362.0, TEX_URANUS, 97.77, -0.71833,
       2.158708519561221E+09, 2.058822953381595E+09, -1.785026690391928E+07,
       -4.813504855508626E+00, 4.665693498302723E+00, 8.548558001407977E-02},
      {"Neptune", 5.151e-5, 24622.0, TEX_NEPTUNE, 28.32, 0.67125,
       2.510985726269384E+09, -3.732999840511448E+09, 1.873087787685826E+07,
       4.445846366126391E+00, 3.067772238252208E+00, -1.638828425643102E-01},
  };
  add_body_init_list(sim, built_in, sizeof(built_in) / sizeof(built_in[0]));
}
