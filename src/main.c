#include <SDL2/SDL.h>

#include <GL/glew.h>

#include <math.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <sys/stat.h>

#include "glutil.h"
#include "render_bg.h"
#include "render_lines.h"
#include "render_particles.h"
#include "render_rings.h"
#include "scene_fbo.h"
#include "presets.h"
#include "input_mouse.h"
#include "shaders.h"
#include "sim.h"
#include "textures.h"
#include "ui_overlay.h"
#include "ui_labels.h"
#include "ui_text.h"
#include "util.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include "camera.h"
#include "mat4.h"
#include "dump_yaml.h"

static bool camera_changed(const Camera *a, const Camera *b) {
  const double eps_pos = 1e-6;
  const double eps_ang = 1e-6;
  const double eps_zoom = 1e-9;
  if (fabs(a->cx - b->cx) > eps_pos)
    return true;
  if (fabs(a->cy - b->cy) > eps_pos)
    return true;
  if (fabs(a->cz - b->cz) > eps_pos)
    return true;
  if (fabs(a->yaw - b->yaw) > eps_ang)
    return true;
  if (fabs(a->pitch - b->pitch) > eps_ang)
    return true;
  if (fabs(a->zoom_world_h - b->zoom_world_h) > eps_zoom)
    return true;
  return false;
}

static bool dir_exists(const char *path) {
  struct stat st;
  if (stat(path, &st) != 0)
    return false;
  return S_ISDIR(st.st_mode);
}

static char *path_join2(const char *a, const char *b) {
  const size_t al = strlen(a);
  const size_t bl = strlen(b);
  const bool need_slash = (al > 0 && a[al - 1] != '/');
  char *out = (char *)malloc(al + bl + (need_slash ? 2 : 1));
  if (!out)
    return NULL;
  memcpy(out, a, al);
  size_t p = al;
  if (need_slash)
    out[p++] = '/';
  memcpy(out + p, b, bl);
  out[p + bl] = '\0';
  return out;
}

static char *resolve_assets_dir(void) {
  const char *env = getenv("UNIVERSE_ASSETS_DIR");
  if (env && env[0] && dir_exists(env))
    return strdup(env);

  // If launched from repo root.
  if (dir_exists("assets"))
    return strdup("assets");

  // If launched from build directory.
  if (dir_exists("../assets"))
    return strdup("../assets");

  // Try relative to the executable.
  char *base = SDL_GetBasePath();
  if (base) {
    char *p1 = path_join2(base, "../assets");
    if (p1 && dir_exists(p1)) {
      SDL_free(base);
      return p1;
    }
    free(p1);
    char *p2 = path_join2(base, "assets");
    if (p2 && dir_exists(p2)) {
      SDL_free(base);
      return p2;
    }
    free(p2);
    SDL_free(base);
  }

  return strdup("assets");
}

static char *resolve_bodies_yaml_path(void) {
  const char *env = getenv("UNIVERSE_BODIES_YAML");
  if (env && env[0])
    return strdup(env);

  // Prefer dumping into build/ since it's already gitignored.
  if (dir_exists("build"))
    return strdup("build/bodies.yaml");

  // If launched from build directory.
  if (dir_exists("../build"))
    return strdup("../build/bodies.yaml");

  return strdup("bodies.yaml");
}

enum {
  TRAIL_LEN = 4096,
};

static size_t trail_len_scale_up(size_t n) {
  if (n < 2)
    n = 2;
  double next = ceil((double)n * 1.5);
  if (next < 2.0)
    next = 2.0;
  if (next > (double)TRAIL_LEN)
    next = (double)TRAIL_LEN;
  return (size_t)next;
}

static size_t trail_len_scale_down(size_t n) {
  if (n < 2)
    n = 2;
  double next = floor((double)n / 1.5);
  if (next < 2.0)
    next = 2.0;
  return (size_t)next;
}

typedef struct {
  uint32_t id;
  float x[TRAIL_LEN];
  float y[TRAIL_LEN];
  float z[TRAIL_LEN];
  size_t head;
  size_t count;
  float last_x;
  float last_y;
  float last_z;
  bool initialized;
} Trail;

static const Body *find_body_by_id(const Sim *sim, uint32_t id) {
  if (!id)
    return NULL;
  for (size_t i = 0; i < sim->count; i++) {
    if (sim->bodies[i].id == id)
      return &sim->bodies[i];
  }
  return NULL;
}

static const Body *find_heaviest_body(const Sim *sim) {
  const Body *best = NULL;
  for (size_t i = 0; i < sim->count; i++) {
    const Body *b = &sim->bodies[i];
    if (!best || b->mass > best->mass)
      best = b;
  }
  return best;
}

static Trail *trails_find(Trail *trails, size_t trail_count, uint32_t id) {
  for (size_t i = 0; i < trail_count; i++) {
    if (trails[i].id == id)
      return &trails[i];
  }
  return NULL;
}

static Trail *trails_get_or_add(Trail **trails, size_t *trail_count,
                                size_t *trail_cap, uint32_t id) {
  Trail *t = trails_find(*trails, *trail_count, id);
  if (t)
    return t;

  if (*trail_count == *trail_cap) {
    size_t next = *trail_cap ? (*trail_cap * 2) : 256;
    Trail *nt = (Trail *)realloc(*trails, next * sizeof(Trail));
    if (!nt)
      return NULL;
    *trails = nt;
    *trail_cap = next;
  }

  Trail *out = &(*trails)[(*trail_count)++];
  memset(out, 0, sizeof(*out));
  out->id = id;
  return out;
}

static void trails_compact_live(Trail *trails, size_t *trail_count,
                                const Sim *sim) {
  size_t out = 0;
  for (size_t i = 0; i < *trail_count; i++) {
    if (find_body_by_id(sim, trails[i].id)) {
      if (out != i)
        trails[out] = trails[i];
      out++;
    }
  }
  *trail_count = out;
}

static void trail_push(Trail *t, float x, float y, float z) {
  if (!t->initialized) {
    t->last_x = x;
    t->last_y = y;
    t->last_z = z;
    t->initialized = true;
  }

  const float dx = x - t->last_x;
  const float dy = y - t->last_y;
  const float dz = z - t->last_z;
  const float d2 = dx * dx + dy * dy + dz * dz;
  // Avoid over-sampling when zoomed out.
  if (d2 < 1e-6f)
    return;

  t->last_x = x;
  t->last_y = y;
  t->last_z = z;

  t->x[t->head] = x;
  t->y[t->head] = y;
  t->z[t->head] = z;
  t->head = (t->head + 1) % TRAIL_LEN;
  if (t->count < TRAIL_LEN)
    t->count++;
}

static void print_controls(void) {
  fprintf(stderr, "Controls:\n"
                  "  F1: toggle help\n"
                  "  F2: toggle body labels\n"
                  "  F3: toggle realistic scaling\n"
                  "  Y: dump bodies to YAML (UNIVERSE_BODIES_YAML=...)\n"
                  "  RMB drag: spawn body with velocity\n"
                  "    (hold Ctrl for circular orbit; Alt reverses)\n"
                  "  LMB click: select nearest body\n"
                  "  LMB drag: rotate camera\n"
                  "  F: follow selected body\n"
                  "  MMB drag / WASD: pan camera\n"
                  "  Mouse wheel: zoom\n"
                  "  Space: pause/resume\n"
                  "  R: reset to solar system\n"
                  "  1: two-body orbit demo\n"
                  "  2: disk galaxy demo\n"
                  "  3: solar system demo\n"
                  "  C: clear all bodies\n"
                  "  M: toggle merge-on-collision\n"
                  "  T: toggle trails\n"
                  "  PgUp/PgDn: trail length\n"
                  "  V: toggle velocity vectors\n"
                  "  [ / ]: decrease/increase spawn mass\n"
                  "  - / =: slow down / speed up time\n"
                  "  G / H: decrease/increase gravity\n"
                  "  Esc: quit\n");
}

int main(int argc, char **argv) {
  (void)argc;
  (void)argv;

  if (SDL_Init(SDL_INIT_VIDEO | SDL_INIT_TIMER) != 0) {
    fprintf(stderr, "SDL_Init failed: %s\n", SDL_GetError());
    return 1;
  }

  // macOS officially supports 4.1 core (not 3.3 core). If we request 3.3 on
  // macOS, SDL may silently fall back to an older context, breaking shaders.
#if defined(__APPLE__)
  SDL_GL_SetAttribute(SDL_GL_CONTEXT_MAJOR_VERSION, 4);
  SDL_GL_SetAttribute(SDL_GL_CONTEXT_MINOR_VERSION, 1);
#else
  SDL_GL_SetAttribute(SDL_GL_CONTEXT_MAJOR_VERSION, 3);
  SDL_GL_SetAttribute(SDL_GL_CONTEXT_MINOR_VERSION, 3);
#endif
  SDL_GL_SetAttribute(SDL_GL_DEPTH_SIZE, 24);
  SDL_GL_SetAttribute(SDL_GL_CONTEXT_PROFILE_MASK, SDL_GL_CONTEXT_PROFILE_CORE);
  SDL_GL_SetAttribute(SDL_GL_DOUBLEBUFFER, 1);
  SDL_GL_SetAttribute(SDL_GL_DEPTH_SIZE, 0);
  SDL_GL_SetAttribute(SDL_GL_STENCIL_SIZE, 0);
  SDL_GL_SetAttribute(SDL_GL_MULTISAMPLEBUFFERS, 0);

  SDL_Window *win = SDL_CreateWindow(
      "Universe Sandbox (C / SDL2 / OpenGL)", SDL_WINDOWPOS_CENTERED,
      SDL_WINDOWPOS_CENTERED, 1280, 720,
      SDL_WINDOW_OPENGL | SDL_WINDOW_RESIZABLE | SDL_WINDOW_ALLOW_HIGHDPI);
  if (!win) {
    fprintf(stderr, "SDL_CreateWindow failed: %s\n", SDL_GetError());
    SDL_Quit();
    return 1;
  }

  SDL_GLContext ctx = SDL_GL_CreateContext(win);
  if (!ctx) {
#if defined(__APPLE__)
    // Try 3.3 core as a fallback (non-mac systems).
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_MAJOR_VERSION, 3);
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_MINOR_VERSION, 3);
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_PROFILE_MASK,
                        SDL_GL_CONTEXT_PROFILE_CORE);
    ctx = SDL_GL_CreateContext(win);
#endif
  }
  if (!ctx) {
    fprintf(stderr, "SDL_GL_CreateContext failed: %s\n", SDL_GetError());
    SDL_DestroyWindow(win);
    SDL_Quit();
    return 1;
  }

  SDL_GL_SetSwapInterval(1);

  glewExperimental = GL_TRUE;
  GLenum glew_err = glewInit();
  // GLEW may emit a harmless GL error during init on core profiles.
  glGetError();
  if (glew_err != GLEW_OK) {
    fprintf(stderr, "glewInit failed: %s\n", glewGetErrorString(glew_err));
    SDL_GL_DeleteContext(ctx);
    SDL_DestroyWindow(win);
    SDL_Quit();
    return 1;
  }

  GLint default_depth_bits = 0;

  // Verify context is new enough for GLSL 330 + texture arrays.
  {
    const char *ver = (const char *)glGetString(GL_VERSION);
    const char *glsl = (const char *)glGetString(GL_SHADING_LANGUAGE_VERSION);
    glGetIntegerv(GL_DEPTH_BITS, &default_depth_bits);
    int maj = 0, min = 0;
    if (ver)
      sscanf(ver, "%d.%d", &maj, &min);
    fprintf(stderr, "GL_VERSION: %s\n", ver ? ver : "(null)");
    fprintf(stderr, "GLSL_VERSION: %s\n", glsl ? glsl : "(null)");
    fprintf(stderr, "DEPTH_BITS: %d\n", (int)default_depth_bits);
    if (maj < 3 || (maj == 3 && min < 3)) {
      fprintf(stderr, "OpenGL 3.3+ required (got %d.%d).\n", maj, min);
      SDL_GL_DeleteContext(ctx);
      SDL_DestroyWindow(win);
      SDL_Quit();
      return 1;
    }
  }

  // Some backends (notably the Metal shim) may expose a default framebuffer
  // without a depth buffer. In that case we render the 3D scene into an
  // offscreen FBO with an explicit depth attachment, then blit to the window.
  bool use_offscreen_scene_fbo = (default_depth_bits == 0);
  SceneFbo scene = {0};

  ParticleRenderer particles;
  if (!particles_init(&particles))
    return 1;

  BgRenderer bg;
  if (!bg_init(&bg))
    return 1;

  UiTextRenderer ui;
  if (!ui_text_renderer_init(&ui))
    return 1;

  LineRenderer lines;
  if (!lines_init(&lines))
    return 1;

  RingRenderer rings;
  if (!rings_init(&rings))
    return 1;

  (void)ui;

  glEnable(GL_DEPTH_TEST);
  glDepthFunc(GL_LESS);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

  char *assets_dir = resolve_assets_dir();
  fprintf(stderr, "Assets dir: %s\n", assets_dir);

  GLuint tex_array = textures_create_planet_array(assets_dir);
  if (!tex_array) {
    fprintf(stderr, "Failed to create texture array\n");
    return 1;
  }

  // Ensure unit 0 has the correct texture target bound for sampler2DArray.
  glActiveTexture(GL_TEXTURE0);
  glBindTexture(GL_TEXTURE_2D_ARRAY, tex_array);

  char *hdr_path = path_join2(assets_dir, "HDR_blue_nebulae-1.hdr");
  GLuint env_tex = hdr_path ? textures_load_hdr_equirect(hdr_path) : 0;
  free(hdr_path);

  // ui_text_renderer_init already created/bound the font texture.

  Sim sim;
  sim_init(&sim);
  {
    int threads = SDL_GetCPUCount();
    const char *env = getenv("UNIVERSE_THREADS");
    if (env && env[0]) {
      threads = atoi(env);
    }

    const char *metal_env = getenv("UNIVERSE_METAL");
    const bool want_metal =
        metal_env && metal_env[0] && strcmp(metal_env, "0") != 0;
    sim_enable_metal(&sim, want_metal);
    if (want_metal) {
      threads = 1;
    }

    if (threads < 1)
      threads = 1;
    if (threads > 16)
      threads = 16;
    sim_set_num_threads(&sim, threads);
    fprintf(stderr, "Gravity threads: %d\n", threads);
    fprintf(stderr, "Metal gravity: %s\n", sim.use_metal ? "on" : "off");
  }

  enum {
    PRESET_TWO_BODY = 1,
    PRESET_GALAXY = 2,
    PRESET_SOLAR = 3,
  };
  int preset = PRESET_SOLAR;
  // Solar system initial conditions are for 2000-01-01 00:00 (JPL/Horizons).
  const double solar_epoch_jd = 2451544.5;

  Camera cam = {.cx = 0.0,
                .cy = 0.0,
                .cz = 0.0,
                .yaw = 0.0,
                .pitch = 0.0,
                .zoom_world_h = 10.0};

  bool running = true;
  bool paused = false;
  bool show_help = false;
  bool show_labels = false;
  bool realistic_scale = false;

  MouseState mouse;
  mouse_state_init(&mouse);
  double spawn_mass = 100.0;

  uint32_t selected_id = 0;
  uint32_t last_selected_printed = 0;
  bool follow_selected = false;
  bool show_trails = true;
  bool show_vectors = true;
  double vector_scale = 0.02;
  double render_radius_scale = 1.6;
  size_t trail_draw_len = TRAIL_LEN;

  Trail *trails = NULL;
  size_t trail_count = 0;
  size_t trail_cap = 0;

  ParticleGPU *particle_cpu = NULL;
  size_t particle_cpu_cap = 0;

  LineVertex *line_cpu = NULL;
  size_t line_cpu_cap = 0;

  TextVertex *text_cpu = NULL;
  size_t text_cpu_cap = 0;

  // Interpreted as "simulation time units per real second".
  // For the solar-system preset, units are years.
  double time_scale = 1.0;
  const double fixed_dt = 1.0e-4;
  double accumulator = 0.0;
  double sim_time = 0.0;

  uint64_t prev_ticks = (uint64_t)SDL_GetPerformanceCounter();
  const uint64_t freq = (uint64_t)SDL_GetPerformanceFrequency();

  uint64_t title_ticks = prev_ticks;
  int frames_since_title = 0;

  print_controls();

  // Apply initial preset from YAML.
  if (!preset_apply_from_yaml("solar_system", &sim, &cam, &time_scale,
                              &selected_id, &follow_selected)) {
    preset_seed_solar_system(&sim);
    preset = PRESET_SOLAR;
    cam = (Camera){.cx = 0.0,
                   .cy = 0.0,
                   .cz = 0.0,
                   .yaw = 0.0,
                   .pitch = -1.24,
                   .zoom_world_h = 4.0};
    time_scale = 0.03;
    selected_id = sim.count ? sim.bodies[0].id : 0;
    follow_selected = false;
  }

  Camera last_cam_printed = cam;
  bool last_cam_printed_init = false;
  uint32_t last_cam_print_ms = 0;

  while (running) {
    SDL_Event e;
    while (SDL_PollEvent(&e)) {
      if (e.type == SDL_QUIT)
        running = false;
      if (e.type == SDL_WINDOWEVENT && e.window.event == SDL_WINDOWEVENT_CLOSE)
        running = false;
      if (e.type == SDL_KEYDOWN && !e.key.repeat) {
        switch (e.key.keysym.sym) {
        case SDLK_F1:
          show_help = !show_help;
          break;
        case SDLK_F2:
          show_labels = !show_labels;
          break;
        case SDLK_F3:
          realistic_scale = !realistic_scale;
          fprintf(stderr, "realistic_scale: %s\n",
                  realistic_scale ? "on" : "off");
          break;
        case SDLK_y: {
          char *out_path = resolve_bodies_yaml_path();
          const bool ok =
              out_path ? dump_bodies_yaml(out_path, &sim, sim_time) : false;
          if (ok) {
            fprintf(stderr, "Wrote YAML: %s\n", out_path);
          } else {
            fprintf(stderr, "Failed to write YAML: %s\n",
                    out_path ? out_path : "(null)");
          }
          free(out_path);
        } break;
        case SDLK_ESCAPE:
          running = false;
          break;
        case SDLK_SPACE:
          paused = !paused;
          break;
        case SDLK_r:
          if (!preset_apply_from_yaml("solar_system", &sim, &cam, &time_scale,
                                      &selected_id, &follow_selected)) {
            preset_seed_solar_system(&sim);
            cam.cx = 0.0;
            cam.cy = 0.0;
            cam.cz = 0.0;
            cam.yaw = 0.0;
            cam.pitch = -1.24, cam.zoom_world_h = 4.0;
            time_scale = 0.03;
            selected_id = sim.count ? sim.bodies[0].id : 0;
            follow_selected = false;
          }
          preset = PRESET_SOLAR;
          sim_time = 0.0;
          break;
        case SDLK_1:
          if (!preset_apply_from_yaml("two_body", &sim, &cam, &time_scale,
                                      &selected_id, &follow_selected)) {
            preset_seed_two_body_orbit(&sim);
            cam.cx = 0.0;
            cam.cy = 0.0;
            cam.cz = 0.0;
            cam.yaw = 0.0;
            cam.pitch = 0.0;
            cam.zoom_world_h = 40.0;
            time_scale = 1.0;
            selected_id = 0;
            follow_selected = false;
          }
          preset = PRESET_TWO_BODY;
          sim_time = 0.0;
          break;
        case SDLK_2:
          if (!preset_apply_from_yaml("disk_galaxy", &sim, &cam, &time_scale,
                                      &selected_id, &follow_selected)) {
            preset_seed_disk_galaxy(&sim, 1500, 140.0, 2.0e6);
            cam.cx = 0.0;
            cam.cy = 0.0;
            cam.cz = 0.0;
            cam.yaw = 0.0;
            cam.pitch = 0.0;
            cam.zoom_world_h = 320.0;
            time_scale = 1.0;
            selected_id = 0;
            follow_selected = false;
          }
          preset = PRESET_GALAXY;
          sim_time = 0.0;
          break;
        case SDLK_3:
          if (!preset_apply_from_yaml("solar_system", &sim, &cam, &time_scale,
                                      &selected_id, &follow_selected)) {
            preset_seed_solar_system(&sim);
            cam.cx = 0.0;
            cam.cy = 0.0;
            cam.cz = 0.0;
            cam.yaw = 0.0;
            cam.pitch = -1.24, cam.zoom_world_h = 4.0;
            time_scale = 0.03;
            selected_id = sim.count ? sim.bodies[0].id : 0;
            follow_selected = false;
          }
          preset = PRESET_SOLAR;
          sim_time = 0.0;
          break;
        case SDLK_c:
          sim_reset(&sim);
          sim_time = 0.0;
          selected_id = 0;
          follow_selected = false;
          break;
        case SDLK_PERIOD:
          if (paused) {
            sim_step(&sim, fixed_dt);
            sim_time += fixed_dt;
          }
          break;
        case SDLK_m:
          sim.merge_on_collision = !sim.merge_on_collision;
          fprintf(stderr, "merge_on_collision: %s\n",
                  sim.merge_on_collision ? "on" : "off");
          break;
        case SDLK_f:
          follow_selected = !follow_selected;
          fprintf(stderr, "follow_selected: %s\n",
                  follow_selected ? "on" : "off");
          break;
        case SDLK_t:
          show_trails = !show_trails;
          fprintf(stderr, "trails: %s\n", show_trails ? "on" : "off");
          break;
        case SDLK_PAGEUP:
          trail_draw_len = trail_len_scale_up(trail_draw_len);
          fprintf(stderr, "trail_draw_len: %zu\n", trail_draw_len);
          break;
        case SDLK_PAGEDOWN:
          trail_draw_len = trail_len_scale_down(trail_draw_len);
          fprintf(stderr, "trail_draw_len: %zu\n", trail_draw_len);
          break;
        case SDLK_v:
          show_vectors = !show_vectors;
          fprintf(stderr, "velocity vectors: %s\n",
                  show_vectors ? "on" : "off");
          break;
        case SDLK_LEFTBRACKET:
          spawn_mass = fmax(1.0, spawn_mass / 1.7);
          fprintf(stderr, "spawn_mass: %.3g\n", spawn_mass);
          break;
        case SDLK_RIGHTBRACKET:
          spawn_mass = fmin(1e12, spawn_mass * 1.7);
          fprintf(stderr, "spawn_mass: %.3g\n", spawn_mass);
          break;
        case SDLK_MINUS:
          time_scale = fmax(1e-7, time_scale / 1.3);
          fprintf(stderr, "time_scale: %.3g\n", time_scale);
          break;
        case SDLK_EQUALS:
          time_scale = fmin(1e6, time_scale * 1.3);
          fprintf(stderr, "time_scale: %.3g\n", time_scale);
          break;
        case SDLK_g:
          sim.G = fmax(1e-6, sim.G / 1.2);
          fprintf(stderr, "G: %.6g\n", sim.G);
          break;
        case SDLK_h:
          sim.G = fmin(1e6, sim.G * 1.2);
          fprintf(stderr, "G: %.6g\n", sim.G);
          break;
        default:
          break;
        }
      }

      if (e.type == SDL_MOUSEWHEEL) {
        mouse_handle_wheel(&e.wheel, &cam);
      }

      if (e.type == SDL_MOUSEBUTTONDOWN) {
        mouse_handle_button_down(&mouse, win, &e.button, &cam, &sim);
      }

      if (e.type == SDL_MOUSEBUTTONUP) {
        mouse_handle_button_up(&mouse, win, &e.button, &cam, &sim,
                               realistic_scale, render_radius_scale, spawn_mass,
                               &selected_id, &follow_selected);
      }

      if (e.type == SDL_MOUSEMOTION) {
        mouse_handle_motion(&mouse, win, &e.motion, &cam);
      }
    }

    const uint8_t *keys = SDL_GetKeyboardState(NULL);
    const float pan_step = (float)(cam.zoom_world_h * 0.01);
    float pvx = 0.0f, pvy = 0.0f;
    if (keys[SDL_SCANCODE_W])
      pvy += pan_step;
    if (keys[SDL_SCANCODE_S])
      pvy -= pan_step;
    if (keys[SDL_SCANCODE_A])
      pvx -= pan_step;
    if (keys[SDL_SCANCODE_D])
      pvx += pan_step;
    if (pvx != 0.0f || pvy != 0.0f) {
      float wx, wy, wz;
      view_to_world_dir(&cam, pvx, pvy, 0.0f, &wx, &wy, &wz);
      cam.cx += (double)wx;
      cam.cy += (double)wy;
      cam.cz += (double)wz;
    }

    // Log camera when it changes (rate-limited so dragging doesn't spam too
    // much).
    {
      const uint32_t now_ms = (uint32_t)SDL_GetTicks();
      const bool changed =
          !last_cam_printed_init || camera_changed(&cam, &last_cam_printed);
      if (changed && (uint32_t)(now_ms - last_cam_print_ms) >= 100) {
        fprintf(stderr,
                "cam: pos=(%.6f %.6f %.6f) yaw=%.6f pitch=%.6f zoom_h=%.6f\n",
                cam.cx, cam.cy, cam.cz, cam.yaw, cam.pitch, cam.zoom_world_h);
        last_cam_printed = cam;
        last_cam_printed_init = true;
        last_cam_print_ms = now_ms;
      }
    }

    const uint64_t now_ticks = (uint64_t)SDL_GetPerformanceCounter();
    const double frame_dt = (double)(now_ticks - prev_ticks) / (double)freq;
    prev_ticks = now_ticks;

    // Avoid giant time steps after debugging pauses, etc.
    accumulator += clampd(frame_dt, 0.0, 0.05) * time_scale;

    if (!paused) {
      int steps = 0;
      while (accumulator >= fixed_dt && steps < 200) {
        sim_step(&sim, fixed_dt);
        sim_time += fixed_dt;
        accumulator -= fixed_dt;
        steps++;
      }
    } else {
      accumulator = 0.0;
    }

    if (selected_id && !find_body_by_id(&sim, selected_id)) {
      selected_id = 0;
      follow_selected = false;
    }

    // Print selection changes.
    if (selected_id != last_selected_printed) {
      if (selected_id) {
        const Body *b = find_body_by_id(&sim, selected_id);
        fprintf(stderr, "selected: %s (id=%u mass=%.6g)\n",
                b ? b->name : "(unknown)", (unsigned)selected_id,
                b ? b->mass : 0.0);
      } else {
        fprintf(stderr, "selection cleared\n");
      }
      last_selected_printed = selected_id;
    }

    if (follow_selected && selected_id) {
      const Body *b = find_body_by_id(&sim, selected_id);
      if (b) {
        cam.cx = b->x;
        cam.cy = b->y;
        cam.cz = b->z;
      }
    }

    trails_compact_live(trails, &trail_count, &sim);
    for (size_t i = 0; i < sim.count; i++) {
      const Body *b = &sim.bodies[i];
      Trail *t = trails_get_or_add(&trails, &trail_count, &trail_cap, b->id);
      if (t)
        trail_push(t, (float)b->x, (float)b->y, (float)b->z);
    }

    int dw = 1, dh = 1;
    SDL_GL_GetDrawableSize(win, &dw, &dh);
    if (use_offscreen_scene_fbo) {
      if (!scene_fbo_ensure(&scene, dw, dh)) {
        use_offscreen_scene_fbo = false;
      }
    }

    glBindFramebuffer(GL_FRAMEBUFFER, use_offscreen_scene_fbo ? scene.fbo : 0);
    glViewport(0, 0, dw, dh);

    glClearColor(0.02f, 0.02f, 0.03f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    Mat4 world_to_clip = camera_world_to_clip(&cam, dw, dh);
    const float pixels_per_world = (float)((double)dh / cam.zoom_world_h);
    float cam_right_x = 1.0f, cam_right_y = 0.0f, cam_right_z = 0.0f;
    float cam_up_x = 0.0f, cam_up_y = 1.0f, cam_up_z = 0.0f;
    view_to_world_dir(&cam, 1.0f, 0.0f, 0.0f, &cam_right_x, &cam_right_y,
                      &cam_right_z);
    view_to_world_dir(&cam, 0.0f, 1.0f, 0.0f, &cam_up_x, &cam_up_y, &cam_up_z);

    // Background (translation-independent)
    bg_draw(&bg, &cam, dw, dh, env_tex);

    // Build line vertex buffer (trails + vectors + spawn preview).
    size_t line_count = 0;
    if (show_trails) {
      const Trail *follow_t = NULL;
      const Body *follow_now = NULL;
      float follow_base_x = 0.0f;
      float follow_base_y = 0.0f;
      float follow_base_z = 0.0f;
      if (follow_selected && selected_id) {
        follow_t = trails_find(trails, trail_count, selected_id);
        follow_now = find_body_by_id(&sim, selected_id);
        if (follow_now) {
          // Trails are computed in the followed body's inertial frame at each
          // sample time, but rendering still uses a camera centered on the
          // followed body *now*. Add the followed body's current position back
          // so the relative trail stays spatially attached in view.
          follow_base_x = (float)follow_now->x;
          follow_base_y = (float)follow_now->y;
          follow_base_z = (float)follow_now->z;
        }
      }
      for (size_t ti = 0; ti < trail_count; ti++) {
        const Trail *t = &trails[ti];
        const Body *b = find_body_by_id(&sim, t->id);
        if (!b || t->count < 2)
          continue;
        if (follow_t && t->id == selected_id)
          continue;

        const size_t t_count =
            (size_t)fmin((double)t->count, (double)trail_draw_len);
        const size_t count =
            follow_t ? (size_t)fmin((double)t_count, (double)follow_t->count)
                     : t_count;
        if (count < 2)
          continue;

        float tr = 1.0f, tg = 1.0f, tb = 1.0f;
        trail_color_from_id(b->id, &tr, &tg, &tb);

        const size_t start = (t->head + TRAIL_LEN - count) % TRAIL_LEN;
        const size_t fstart =
            follow_t ? ((follow_t->head + TRAIL_LEN - count) % TRAIL_LEN) : 0;
        for (size_t k = 0; k + 1 < count; k++) {
          const size_t i0 = (start + k) % TRAIL_LEN;
          const size_t i1 = (start + k + 1) % TRAIL_LEN;
          const size_t f0 = follow_t ? ((fstart + k) % TRAIL_LEN) : 0;
          const size_t f1 = follow_t ? ((fstart + k + 1) % TRAIL_LEN) : 0;
          const float age0 = (float)k / (float)(count - 1);
          const float age1 = (float)(k + 1) / (float)(count - 1);
          const float a0 = 0.55f * age0 * age0;
          const float a1 = 0.55f * age1 * age1;

          const float x0 = follow_t
                               ? (t->x[i0] - follow_t->x[f0] + follow_base_x)
                               : t->x[i0];
          const float y0 = follow_t
                               ? (t->y[i0] - follow_t->y[f0] + follow_base_y)
                               : t->y[i0];
          const float z0 = follow_t
                               ? (t->z[i0] - follow_t->z[f0] + follow_base_z)
                               : t->z[i0];
          const float x1 = follow_t
                               ? (t->x[i1] - follow_t->x[f1] + follow_base_x)
                               : t->x[i1];
          const float y1 = follow_t
                               ? (t->y[i1] - follow_t->y[f1] + follow_base_y)
                               : t->y[i1];
          const float z1 = follow_t
                               ? (t->z[i1] - follow_t->z[f1] + follow_base_z)
                               : t->z[i1];

          if (line_count + 2 > line_cpu_cap) {
            size_t next = line_cpu_cap ? (line_cpu_cap * 2) : 16384;
            while (next < line_count + 2)
              next *= 2;
            LineVertex *nl =
                (LineVertex *)realloc(line_cpu, next * sizeof(LineVertex));
            if (!nl)
              break;
            line_cpu = nl;
            line_cpu_cap = next;
          }

          line_cpu[line_count++] = (LineVertex){
              .x = x0, .y = y0, .z = z0, .r = tr, .g = tg, .b = tb, .a = a0};
          line_cpu[line_count++] = (LineVertex){
              .x = x1, .y = y1, .z = z1, .r = tr, .g = tg, .b = tb, .a = a1};
        }
      }
    }

    if (show_vectors) {
      for (size_t i = 0; i < sim.count; i++) {
        const Body *b = &sim.bodies[i];
        const double x0 = b->x;
        const double y0 = b->y;
        const double z0 = b->z;
        const double x1 = b->x + b->vx * vector_scale;
        const double y1 = b->y + b->vy * vector_scale;
        const double z1 = b->z + b->vz * vector_scale;
        const float cr = b->r;
        const float cg = b->g;
        const float cb = b->b;
        const float ca = 0.4f;

        if (line_count + 2 > line_cpu_cap) {
          size_t next = line_cpu_cap ? (line_cpu_cap * 2) : 16384;
          while (next < line_count + 2)
            next *= 2;
          LineVertex *nl =
              (LineVertex *)realloc(line_cpu, next * sizeof(LineVertex));
          if (!nl)
            break;
          line_cpu = nl;
          line_cpu_cap = next;
        }
        line_cpu[line_count++] = (LineVertex){.x = (float)x0,
                                              .y = (float)y0,
                                              .z = (float)z0,
                                              .r = cr,
                                              .g = cg,
                                              .b = cb,
                                              .a = ca};
        line_cpu[line_count++] = (LineVertex){.x = (float)x1,
                                              .y = (float)y1,
                                              .z = (float)z1,
                                              .r = cr,
                                              .g = cg,
                                              .b = cb,
                                              .a = ca};
      }
    }

    if (selected_id) {
      const Body *sel = find_body_by_id(&sim, selected_id);
      if (sel) {
        // Camera-facing selection ring, slightly bigger than the sprite.
        float right_x, right_y, right_z;
        float up_x, up_y, up_z;
        view_to_world_dir(&cam, 1.0f, 0.0f, 0.0f, &right_x, &right_y, &right_z);
        view_to_world_dir(&cam, 0.0f, 1.0f, 0.0f, &up_x, &up_y, &up_z);

        const int segs = 128;
        const double body_r_world =
            body_visual_radius_world(sel, realistic_scale, render_radius_scale);
        const double body_r_px_raw = body_r_world * (double)pixels_per_world;
        const double body_r_px = body_r_px_raw;
        // Keep the ring just a few pixels larger than the sprite.
        // Avoid a large minimum radius which makes small bodies look
        // "over-selected".
        const double ring0_px = clampd(body_r_px + 3.0, 3.0, body_r_px + 10.0);
        const double ring1_px = ring0_px + 2.0;
        const double rr0 = ring0_px / (double)pixels_per_world;
        const double rr1 = ring1_px / (double)pixels_per_world;

        const float cr = 1.0f;
        const float cg = 1.0f;
        const float cb = 1.0f;

        for (int pass = 0; pass < 2; pass++) {
          const double rr = (pass == 0) ? rr0 : rr1;
          const float ca = (pass == 0) ? 0.85f : 0.25f;
          for (int s = 0; s < segs; s++) {
            const double a0 = (double)s * (2.0 * M_PI) / (double)segs;
            const double a1 = (double)(s + 1) * (2.0 * M_PI) / (double)segs;

            const float ox0 = (float)(cos(a0) * rr);
            const float oy0 = (float)(sin(a0) * rr);
            const float ox1 = (float)(cos(a1) * rr);
            const float oy1 = (float)(sin(a1) * rr);

            const float x0 = (float)sel->x + right_x * ox0 + up_x * oy0;
            const float y0 = (float)sel->y + right_y * ox0 + up_y * oy0;
            const float z0 = (float)sel->z + right_z * ox0 + up_z * oy0;
            const float x1 = (float)sel->x + right_x * ox1 + up_x * oy1;
            const float y1 = (float)sel->y + right_y * ox1 + up_y * oy1;
            const float z1 = (float)sel->z + right_z * ox1 + up_z * oy1;

            if (line_count + 2 > line_cpu_cap) {
              size_t next = line_cpu_cap ? (line_cpu_cap * 2) : 16384;
              while (next < line_count + 2)
                next *= 2;
              LineVertex *nl =
                  (LineVertex *)realloc(line_cpu, next * sizeof(LineVertex));
              if (!nl)
                break;
              line_cpu = nl;
              line_cpu_cap = next;
            }
            line_cpu[line_count++] = (LineVertex){
                .x = x0, .y = y0, .z = z0, .r = cr, .g = cg, .b = cb, .a = ca};
            line_cpu[line_count++] = (LineVertex){
                .x = x1, .y = y1, .z = z1, .r = cr, .g = cg, .b = cb, .a = ca};
          }
        }
      }
    }

    if (mouse.spawning) {
      int mx = 0, my = 0;
      SDL_GetMouseState(&mx, &my);
      int mdx = 0, mdy = 0;
      window_to_drawable(win, mx, my, &mdx, &mdy);
      double mwx = 0.0, mwy = 0.0, mwz = 0.0;
      screen_to_world_on_z0(&cam, dw, dh, mdx, mdy, &mwx, &mwy, &mwz);
      if (line_count + 2 > line_cpu_cap) {
        size_t next = line_cpu_cap ? (line_cpu_cap * 2) : 16384;
        while (next < line_count + 2)
          next *= 2;
        LineVertex *nl =
            (LineVertex *)realloc(line_cpu, next * sizeof(LineVertex));
        if (nl) {
          line_cpu = nl;
          line_cpu_cap = next;
        }
      }
      if (line_count + 2 <= line_cpu_cap) {
        line_cpu[line_count++] = (LineVertex){
            .x = (float)mouse.spawn_start_wx,
            .y = (float)mouse.spawn_start_wy,
            .z = (float)mouse.spawn_start_wz,
            .r = 1.0f,
            .g = 1.0f,
            .b = 1.0f,
            .a = 0.9f};
        line_cpu[line_count++] = (LineVertex){.x = (float)mwx,
                                              .y = (float)mwy,
                                              .z = (float)mwz,
                                              .r = 1.0f,
                                              .g = 1.0f,
                                              .b = 1.0f,
                                              .a = 0.9f};
      }
    }

    if (line_count) {
      lines_draw(&lines, &world_to_clip, line_cpu, line_count);
    }

    // Build particle buffer.
    if (sim.count > particle_cpu_cap) {
      size_t next = particle_cpu_cap ? particle_cpu_cap : 1024;
      while (next < sim.count)
        next *= 2;
      ParticleGPU *np =
          (ParticleGPU *)realloc(particle_cpu, next * sizeof(ParticleGPU));
      if (np) {
        particle_cpu = np;
        particle_cpu_cap = next;
      }
    }
    for (size_t i = 0; i < sim.count; i++) {
      const Body *b = &sim.bodies[i];
      particle_cpu[i].x = (float)b->x;
      particle_cpu[i].y = (float)b->y;
      particle_cpu[i].z = (float)b->z;
      particle_cpu[i].radius = (float)body_visual_radius_world(
          b, realistic_scale, render_radius_scale);
      particle_cpu[i].r = b->r;
      particle_cpu[i].g = b->g;
      particle_cpu[i].b = b->b;
      particle_cpu[i].layer = b->tex_layer;
      particle_cpu[i].tilt = b->tilt_rad;
      particle_cpu[i].spin_rate = b->spin_rate_rad_per_time;
      particle_cpu[i].spin_phase = b->spin_phase_rad;
    }

    particles_draw(&particles, &world_to_clip, (float)sim_time, cam_right_x,
                   cam_right_y, cam_right_z, cam_up_x, cam_up_y, cam_up_z,
                   tex_array, particle_cpu, sim.count);

    rings_build_and_draw(&rings, &world_to_clip, &sim, assets_dir,
                         realistic_scale, render_radius_scale);

    if (use_offscreen_scene_fbo) {
      glBindFramebuffer(GL_READ_FRAMEBUFFER, scene.fbo);
      glBindFramebuffer(GL_DRAW_FRAMEBUFFER, 0);
      glBlitFramebuffer(0, 0, dw, dh, 0, 0, dw, dh, GL_COLOR_BUFFER_BIT,
                        GL_NEAREST);
      glBindFramebuffer(GL_FRAMEBUFFER, 0);
      glViewport(0, 0, dw, dh);
    }

    // Inertial orbit inset: show selected body's trail in a Sun-centered frame
    // while following (main view becomes body-centered).
    if (follow_selected && selected_id) {
      const Body *sel = find_body_by_id(&sim, selected_id);
      const Body *center = find_heaviest_body(&sim);
      Trail *t = trails_find(trails, trail_count, selected_id);
      if (sel && center && t && t->count >= 2) {
        const int inset_w = (int)clampd((double)dw * 0.28, 220.0, 420.0);
        const int inset_h = inset_w;
        const int inset_x = dw - inset_w - 14;
        const int inset_y = dh - inset_h - 14;

        // Determine a suitable scale from the trail extents.
        double max_r = 1.0;
        for (size_t k = 0; k < t->count; k++) {
          const size_t idx = (t->head + TRAIL_LEN - t->count + k) % TRAIL_LEN;
          const double dx = (double)t->x[idx] - center->x;
          const double dy = (double)t->y[idx] - center->y;
          const double r = sqrt(dx * dx + dy * dy);
          if (r > max_r)
            max_r = r;
        }
        max_r = fmax(max_r, 0.01);
        const double half_world = max_r * 1.15;

        glDisable(GL_DEPTH_TEST);
        glDepthMask(GL_FALSE);
        glEnable(GL_SCISSOR_TEST);
        glScissor(inset_x, inset_y, inset_w, inset_h);
        glClearColor(0.02f, 0.02f, 0.03f, 0.85f);
        glClear(GL_COLOR_BUFFER_BIT);
        glDisable(GL_SCISSOR_TEST);

        glViewport(inset_x, inset_y, inset_w, inset_h);
        Mat4 inset_w2c = mat4_inset_world_to_clip(inset_w, inset_h, center->x,
                                                  center->y, half_world);

        // Build line list: orbit trail + small cross for Sun and current pos.
        size_t inset_line_count = 0;
        const size_t needed = (t->count - 1) * 2 + 8;
        if (needed > line_cpu_cap) {
          size_t next = line_cpu_cap ? line_cpu_cap : 16384;
          while (next < needed)
            next *= 2;
          LineVertex *nl =
              (LineVertex *)realloc(line_cpu, next * sizeof(LineVertex));
          if (nl) {
            line_cpu = nl;
            line_cpu_cap = next;
          }
        }

        if (needed <= line_cpu_cap) {
          for (size_t k = 0; k + 1 < t->count; k++) {
            const size_t i0 = (t->head + TRAIL_LEN - t->count + k) % TRAIL_LEN;
            const size_t i1 =
                (t->head + TRAIL_LEN - t->count + k + 1) % TRAIL_LEN;
            line_cpu[inset_line_count++] = (LineVertex){.x = t->x[i0],
                                                        .y = t->y[i0],
                                                        .z = 0.0f,
                                                        .r = 0.6f,
                                                        .g = 0.8f,
                                                        .b = 1.0f,
                                                        .a = 0.75f};
            line_cpu[inset_line_count++] = (LineVertex){.x = t->x[i1],
                                                        .y = t->y[i1],
                                                        .z = 0.0f,
                                                        .r = 0.6f,
                                                        .g = 0.8f,
                                                        .b = 1.0f,
                                                        .a = 0.75f};
          }

          const float cross = (float)(half_world * 0.03);
          const float cx = (float)center->x;
          const float cy = (float)center->y;
          const float px = (float)sel->x;
          const float py = (float)sel->y;

          // Sun cross
          line_cpu[inset_line_count++] = (LineVertex){.x = cx - cross,
                                                      .y = cy,
                                                      .z = 0.0f,
                                                      .r = 1.0f,
                                                      .g = 0.9f,
                                                      .b = 0.4f,
                                                      .a = 0.9f};
          line_cpu[inset_line_count++] = (LineVertex){.x = cx + cross,
                                                      .y = cy,
                                                      .z = 0.0f,
                                                      .r = 1.0f,
                                                      .g = 0.9f,
                                                      .b = 0.4f,
                                                      .a = 0.9f};
          line_cpu[inset_line_count++] = (LineVertex){.x = cx,
                                                      .y = cy - cross,
                                                      .z = 0.0f,
                                                      .r = 1.0f,
                                                      .g = 0.9f,
                                                      .b = 0.4f,
                                                      .a = 0.9f};
          line_cpu[inset_line_count++] = (LineVertex){.x = cx,
                                                      .y = cy + cross,
                                                      .z = 0.0f,
                                                      .r = 1.0f,
                                                      .g = 0.9f,
                                                      .b = 0.4f,
                                                      .a = 0.9f};

          // Current body cross
          line_cpu[inset_line_count++] = (LineVertex){.x = px - cross,
                                                      .y = py,
                                                      .z = 0.0f,
                                                      .r = 1.0f,
                                                      .g = 1.0f,
                                                      .b = 1.0f,
                                                      .a = 0.95f};
          line_cpu[inset_line_count++] = (LineVertex){.x = px + cross,
                                                      .y = py,
                                                      .z = 0.0f,
                                                      .r = 1.0f,
                                                      .g = 1.0f,
                                                      .b = 1.0f,
                                                      .a = 0.95f};
          line_cpu[inset_line_count++] = (LineVertex){.x = px,
                                                      .y = py - cross,
                                                      .z = 0.0f,
                                                      .r = 1.0f,
                                                      .g = 1.0f,
                                                      .b = 1.0f,
                                                      .a = 0.95f};
          line_cpu[inset_line_count++] = (LineVertex){.x = px,
                                                      .y = py + cross,
                                                      .z = 0.0f,
                                                      .r = 1.0f,
                                                      .g = 1.0f,
                                                      .b = 1.0f,
                                                      .a = 0.95f};
        }

        if (inset_line_count) {
          lines_draw(&lines, &inset_w2c, line_cpu, inset_line_count);
        }

        // Restore main viewport for screen-space overlays.
        glViewport(0, 0, dw, dh);
        glDepthMask(GL_TRUE);
        glEnable(GL_DEPTH_TEST);
      }
    }

    frames_since_title++;
    const double title_dt = (double)(now_ticks - title_ticks) / (double)freq;
    if (title_dt >= 0.25) {
      const double fps =
          (title_dt > 0.0) ? ((double)frames_since_title / title_dt) : 0.0;
      char title[256];
      const Body *sel = selected_id ? find_body_by_id(&sim, selected_id) : NULL;
      if (sel) {
        snprintf(title, sizeof(title),
                 "Universe Sandbox (bodies=%zu fps=%.1f sel=%s m=%.3g)%s",
                 sim.count, fps, sel->name, sel->mass,
                 follow_selected ? " [FOLLOW]" : "");
      } else {
        snprintf(title, sizeof(title),
                 "Universe Sandbox (bodies=%zu fps=%.1f)%s", sim.count, fps,
                 follow_selected ? " [FOLLOW]" : "");
      }
      SDL_SetWindowTitle(win, title);
      title_ticks = now_ticks;
      frames_since_title = 0;
    }

    UiOverlayConfig ucfg = {.solar_epoch_jd = solar_epoch_jd};
    UiOverlayState ust = {.preset = preset,
                          .sim_time = sim_time,
                          .show_help = show_help,
                          .show_labels = show_labels,
                          .follow_selected = follow_selected,
                          .realistic_scale = realistic_scale,
                          .cam = cam,
                          .default_depth_bits = (int)default_depth_bits,
                          .use_offscreen_scene_fbo = use_offscreen_scene_fbo,
                          .scene_depth_bits = scene.depth_bits};
    ui_overlay_draw(&ui, &ucfg, &ust, dw, dh);

    if (selected_id) {
      const Body *sel = find_body_by_id(&sim, selected_id);
      if (sel) {
        ui_draw_selected_hud(&ui, &text_cpu, &text_cpu_cap, dw, dh,
                             &world_to_clip, &cam, sel, realistic_scale,
                             render_radius_scale);
      }
    }

    if (show_labels) {
      ui_draw_body_labels(&ui, &text_cpu, &text_cpu_cap, dw, dh,
                          &world_to_clip, &cam, &sim, realistic_scale,
                          render_radius_scale);
    }

    // Other overlays (help, camera HUD, depth HUD, scaling mode) moved to
    // ui_overlay_draw().

    SDL_GL_SwapWindow(win);
  }

  sim_destroy(&sim);
  free(trails);
  free(particle_cpu);
  free(line_cpu);
  free(text_cpu);
  particles_destroy(&particles);
  rings_destroy(&rings);
  lines_destroy(&lines);
  ui_text_renderer_destroy(&ui);
  glDeleteTextures(1, &tex_array);
  glDeleteTextures(1, &env_tex);
  scene_fbo_destroy(&scene);
  bg_destroy(&bg);
  free(assets_dir);
  SDL_GL_DeleteContext(ctx);
  SDL_DestroyWindow(win);
  SDL_Quit();
  return 0;
}
