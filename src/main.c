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
#include "shaders.h"
#include "sim.h"
#include "textures.h"

#include "font8x8_basic.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

typedef struct {
  double cx;
  double cy;
  double cz;
  double yaw;
  double pitch;
  double zoom_world_h;
} Camera;

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

static double body_visual_radius_world(const Body *b, bool realistic_scale,
                                      double exaggerated_scale) {
  const double base = realistic_scale ? b->radius : b->render_radius;
  const double scale = realistic_scale ? 1.0 : exaggerated_scale;
  return base * scale;
}

typedef struct {
  float x;
  float y;
  float z;
  float r;
  float g;
  float b;
  float a;
} LineVertex;

static double clampd(double x, double lo, double hi) {
  if (x < lo)
    return lo;
  if (x > hi)
    return hi;
  return x;
}

typedef struct {
  float x;
  float y;
  float u;
  float v;
} TextVertex;

static GLuint create_font_texture(void) {
  // 16x8 glyph grid of 8x8 pixels.
  const int gw = 8;
  const int gh = 8;
  const int cols = 16;
  const int rows = 8;
  const int w = gw * cols;
  const int h = gh * rows;
  uint8_t *img = (uint8_t *)calloc((size_t)w * (size_t)h, 1);
  if (!img)
    return 0;

  for (int c = 0; c < 128; c++) {
    const int gx = (c % cols) * gw;
    const int gy = (c / cols) * gh;
    for (int y = 0; y < 8; y++) {
      const uint8_t row = font8x8_basic[c][y];
      for (int x = 0; x < 8; x++) {
        // font8x8_basic stores pixels MSB->LSB left-to-right.
        const int bit = (row >> (7 - x)) & 1;
        img[(gy + y) * w + (gx + x)] = bit ? 255 : 0;
      }
    }
  }

  GLuint tex = 0;
  glGenTextures(1, &tex);
  glBindTexture(GL_TEXTURE_2D, tex);
  glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
  glTexImage2D(GL_TEXTURE_2D, 0, GL_R8, w, h, 0, GL_RED, GL_UNSIGNED_BYTE, img);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

  free(img);
  return tex;
}

static void text_append(TextVertex **verts, size_t *count, size_t *cap, float x,
                        float y, float scale, const char *text) {
  const float gw = 8.0f * scale;
  const float gh = 8.0f * scale;
  const float atlas_w = 128.0f;
  const float atlas_h = 64.0f;

  for (const unsigned char *p = (const unsigned char *)text; *p; p++) {
    const unsigned char ch = *p;
    if (ch == '\n') {
      x = 0.0f;
      y += gh;
      continue;
    }
    if (ch < 32) {
      x += gw;
      continue;
    }

    const int cell_x = (int)(ch % 16);
    const int cell_y = (int)(ch / 16);
    const float u0 = (float)(cell_x * 8) / atlas_w;
    const float v0 = (float)(cell_y * 8) / atlas_h;
    const float u1 = (float)((cell_x + 1) * 8) / atlas_w;
    const float v1 = (float)((cell_y + 1) * 8) / atlas_h;

    if (*count + 6 > *cap) {
      size_t next = (*cap) ? (*cap) * 2 : 4096;
      while (next < *count + 6)
        next *= 2;
      TextVertex *nv = (TextVertex *)realloc(*verts, next * sizeof(TextVertex));
      if (!nv)
        return;
      *verts = nv;
      *cap = next;
    }

    const float x0 = x;
    const float y0 = y;
    const float x1 = x + gw;
    const float y1 = y + gh;
    // Two triangles.
    (*verts)[(*count)++] = (TextVertex){x0, y0, u0, v0};
    (*verts)[(*count)++] = (TextVertex){x1, y0, u1, v0};
    (*verts)[(*count)++] = (TextVertex){x1, y1, u1, v1};
    (*verts)[(*count)++] = (TextVertex){x0, y0, u0, v0};
    (*verts)[(*count)++] = (TextVertex){x1, y1, u1, v1};
    (*verts)[(*count)++] = (TextVertex){x0, y1, u0, v1};

    x += gw;
  }
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

typedef struct {
  float m[16];
} Mat4;

static Mat4 mat4_identity(void) {
  Mat4 out = {.m = {0}};
  out.m[0] = 1.0f;
  out.m[5] = 1.0f;
  out.m[10] = 1.0f;
  out.m[15] = 1.0f;
  return out;
}

static Mat4 mat4_mul(Mat4 a, Mat4 b) {
  Mat4 out = {.m = {0}};
  for (int c = 0; c < 4; c++) {
    for (int r = 0; r < 4; r++) {
      out.m[c * 4 + r] =
          a.m[0 * 4 + r] * b.m[c * 4 + 0] + a.m[1 * 4 + r] * b.m[c * 4 + 1] +
          a.m[2 * 4 + r] * b.m[c * 4 + 2] + a.m[3 * 4 + r] * b.m[c * 4 + 3];
    }
  }
  return out;
}

static Mat4 mat4_translate(float x, float y, float z) {
  Mat4 out = mat4_identity();
  out.m[12] = x;
  out.m[13] = y;
  out.m[14] = z;
  return out;
}

static Mat4 mat4_rot_x(float a) {
  Mat4 out = mat4_identity();
  const float c = cosf(a);
  const float s = sinf(a);
  out.m[5] = c;
  out.m[6] = s;
  out.m[9] = -s;
  out.m[10] = c;
  return out;
}

static Mat4 mat4_rot_z(float a) {
  Mat4 out = mat4_identity();
  const float c = cosf(a);
  const float s = sinf(a);
  out.m[0] = c;
  out.m[1] = s;
  out.m[4] = -s;
  out.m[5] = c;
  return out;
}

static Mat4 mat4_ortho(float left, float right, float bottom, float top,
                       float near_z, float far_z) {
  Mat4 out = {.m = {0}};
  out.m[0] = 2.0f / (right - left);
  out.m[5] = 2.0f / (top - bottom);
  out.m[10] = -2.0f / (far_z - near_z);
  out.m[12] = -(right + left) / (right - left);
  out.m[13] = -(top + bottom) / (top - bottom);
  out.m[14] = -(far_z + near_z) / (far_z - near_z);
  out.m[15] = 1.0f;
  return out;
}

static void mat4_mul_vec4(const Mat4 *m, float x, float y, float z, float w,
                          float out4[4]) {
  out4[0] = m->m[0] * x + m->m[4] * y + m->m[8] * z + m->m[12] * w;
  out4[1] = m->m[1] * x + m->m[5] * y + m->m[9] * z + m->m[13] * w;
  out4[2] = m->m[2] * x + m->m[6] * y + m->m[10] * z + m->m[14] * w;
  out4[3] = m->m[3] * x + m->m[7] * y + m->m[11] * z + m->m[15] * w;
}

static bool world_to_screen_px(const Mat4 *world_to_clip, int w, int h, float x,
                               float y, float z, float *sx, float *sy) {
  float clip[4];
  mat4_mul_vec4(world_to_clip, x, y, z, 1.0f, clip);
  if (fabsf(clip[3]) < 1e-8f)
    return false;
  const float ndc_x = clip[0] / clip[3];
  const float ndc_y = clip[1] / clip[3];
  *sx = (ndc_x * 0.5f + 0.5f) * (float)w;
  *sy = (0.5f - ndc_y * 0.5f) * (float)h;
  return true;
}

static void window_to_drawable(SDL_Window *win, int wx, int wy, int *dx,
                               int *dy) {
  int ww = 1, wh = 1;
  int dw = 1, dh = 1;
  SDL_GetWindowSize(win, &ww, &wh);
  SDL_GL_GetDrawableSize(win, &dw, &dh);
  const float sx = (ww > 0) ? ((float)dw / (float)ww) : 1.0f;
  const float sy = (wh > 0) ? ((float)dh / (float)wh) : 1.0f;
  *dx = (int)lrintf((float)wx * sx);
  *dy = (int)lrintf((float)wy * sy);
}

static Mat4 camera_world_to_clip(const Camera *cam, int w, int h) {
  const double aspect = (h > 0) ? ((double)w / (double)h) : 1.0;
  const float half_h = (float)(cam->zoom_world_h * 0.5);
  const float half_w = (float)(cam->zoom_world_h * aspect * 0.5);
  // Keep the depth range reasonably tight for better depth precision,
  // otherwise small Z differences (e.g. planet behind the Sun) quantize away.
  const float depth_range = (float)fmax(1000.0, cam->zoom_world_h * 100.0);
  Mat4 proj =
      mat4_ortho(-half_w, half_w, -half_h, half_h, -depth_range, depth_range);

  Mat4 t = mat4_translate((float)-cam->cx, (float)-cam->cy, (float)-cam->cz);
  Mat4 rz = mat4_rot_z((float)cam->yaw);
  Mat4 rx = mat4_rot_x((float)cam->pitch);
  Mat4 view = mat4_mul(rx, mat4_mul(rz, t));
  return mat4_mul(proj, view);
}

static void view_to_world_dir(const Camera *cam, float vx, float vy, float vz,
                              float *wx, float *wy, float *wz) {
  // Inverse of view rotation: world = Rz(-yaw)*Rx(-pitch)*view
  const float cp = cosf((float)-cam->pitch);
  const float sp = sinf((float)-cam->pitch);
  float x1 = vx;
  float y1 = cp * vy + -sp * vz;
  float z1 = sp * vy + cp * vz;

  const float cy = cosf((float)-cam->yaw);
  const float sy = sinf((float)-cam->yaw);
  *wx = cy * x1 + -sy * y1;
  *wy = sy * x1 + cy * y1;
  *wz = z1;
}

static void screen_to_world_on_z0(const Camera *cam, int w, int h, int sx,
                                  int sy, double *wx, double *wy, double *wz) {
  const double aspect = (h > 0) ? ((double)w / (double)h) : 1.0;
  const double x_ndc = 2.0 * ((double)sx / (double)w) - 1.0;
  const double y_ndc = 1.0 - 2.0 * ((double)sy / (double)h);
  const float half_h = (float)(cam->zoom_world_h * 0.5);
  const float half_w = (float)(cam->zoom_world_h * aspect * 0.5);

  // View-space position on the z=0 slice.
  const float vx = (float)(x_ndc * half_w);
  const float vy = (float)(y_ndc * half_h);

  // World-space ray: origin at that view point, direction is +viewZ.
  float ox, oy, oz;
  view_to_world_dir(cam, vx, vy, 0.0f, &ox, &oy, &oz);
  ox += (float)cam->cx;
  oy += (float)cam->cy;
  oz += (float)cam->cz;

  float dx, dy, dz;
  view_to_world_dir(cam, 0.0f, 0.0f, 1.0f, &dx, &dy, &dz);
  if (fabsf(dz) < 1e-8f) {
    *wx = ox;
    *wy = oy;
    *wz = 0.0;
    return;
  }
  const float t = -oz / dz;
  *wx = ox + dx * t;
  *wy = oy + dy * t;
  *wz = 0.0;
}

static void seed_two_body_orbit(Sim *sim) {
  sim_reset(sim);
  sim->G = 1.0;
  sim->softening = 0.02;
  sim->density = 2.0;
  sim->merge_on_collision = true;

  const double sun_mass = 1.0e6;
  const double planet_mass = 10.0;
  const double r = 30.0;
  const double v = sqrt(sim->G * sun_mass / r);

  (void)sim_add_body(sim, (Body){.name = "Primary",
                                 .x = 0.0,
                                 .y = 0.0,
                                 .vx = 0.0,
                                 .vy = 0.0,
                                 .mass = sun_mass,
                                 .r = 1.0f,
                                 .g = 0.9f,
                                 .b = 0.6f});
  (void)sim_add_body(sim, (Body){.name = "Secondary",
                                 .x = r,
                                 .y = 0.0,
                                 .vx = 0.0,
                                 .vy = v,
                                 .mass = planet_mass,
                                 .r = 0.5f,
                                 .g = 0.8f,
                                 .b = 1.0f});
}

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

static void seed_disk_galaxy(Sim *sim, int n, double disk_radius,
                             double central_mass) {
  sim_reset(sim);
  sim->G = 1.0;
  sim->softening = 0.05;
  sim->density = 1.5;
  sim->merge_on_collision = true;

  (void)sim_add_body(sim, (Body){.name = "Center",
                                 .x = 0.0,
                                 .y = 0.0,
                                 .vx = 0.0,
                                 .vy = 0.0,
                                 .mass = central_mass,
                                 .r = 1.0f,
                                 .g = 0.85f,
                                 .b = 0.5f});

  uint32_t rng = 0x12345678u;
  for (int i = 0; i < n; i++) {
    const double u = frand01(&rng);
    const double v = frand01(&rng);
    const double r = disk_radius * sqrt(u);
    const double a = v * 2.0 * M_PI;
    const double x = r * cos(a);
    const double y = r * sin(a);

    // Circular-ish orbit around center with slight noise.
    const double speed = sqrt(sim->G * central_mass / fmax(r, 1e-3));
    const double tx = -sin(a);
    const double ty = cos(a);
    const double noise = (frand01(&rng) - 0.5) * 0.08;

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

static void jd_to_gregorian_utc(double jd, int *year, int *month, int *day,
                                int *hour, int *min, int *sec) {
  // Astronomical Julian Day -> proleptic Gregorian calendar.
  // This is sufficient for an on-screen simulation clock.
  const double jd0 = jd + 0.5;
  long z = (long)floor(jd0);
  double f = jd0 - (double)z;

  long a = z;
  if (z >= 2299161) {
    const long alpha = (long)floor(((double)z - 1867216.25) / 36524.25);
    a = z + 1 + alpha - alpha / 4;
  }
  const long b = a + 1524;
  const long c = (long)floor(((double)b - 122.1) / 365.25);
  const long d = (long)floor(365.25 * (double)c);
  const long e = (long)floor(((double)b - (double)d) / 30.6001);

  const double day_f = (double)b - (double)d - floor(30.6001 * (double)e) + f;
  const int day_i = (int)floor(day_f);
  const double frac = day_f - (double)day_i;

  int m = (e < 14) ? (int)(e - 1) : (int)(e - 13);
  int y = (m > 2) ? (int)(c - 4716) : (int)(c - 4715);

  const double h_f = frac * 24.0;
  int h = (int)floor(h_f);
  const double min_f = (h_f - (double)h) * 60.0;
  int mi = (int)floor(min_f);
  int s = (int)lrint((min_f - (double)mi) * 60.0);
  if (s >= 60) {
    s = 0;
    mi++;
  }
  if (mi >= 60) {
    mi = 0;
    h++;
  }
  if (h >= 24) {
    // Rare (rounding) overflow; keep it simple.
    h = 0;
  }

  *year = y;
  *month = m;
  *day = day_i;
  *hour = h;
  *min = mi;
  *sec = s;
}

static void scene_fbo_destroy(GLuint *fbo, GLuint *color_tex, GLuint *depth_rb,
                              int *w, int *h) {
  if (*depth_rb)
    glDeleteRenderbuffers(1, depth_rb);
  if (*color_tex)
    glDeleteTextures(1, color_tex);
  if (*fbo)
    glDeleteFramebuffers(1, fbo);
  *fbo = 0;
  *color_tex = 0;
  *depth_rb = 0;
  *w = 0;
  *h = 0;
}

static bool scene_fbo_ensure(GLuint *fbo, GLuint *color_tex, GLuint *depth_rb,
                             int *cur_w, int *cur_h, int *depth_bits, int w,
                             int h) {
  if (w <= 0 || h <= 0)
    return false;
  if (*fbo && *cur_w == w && *cur_h == h)
    return true;

  scene_fbo_destroy(fbo, color_tex, depth_rb, cur_w, cur_h);

  glGenFramebuffers(1, fbo);
  glBindFramebuffer(GL_FRAMEBUFFER, *fbo);

  glGenTextures(1, color_tex);
  glBindTexture(GL_TEXTURE_2D, *color_tex);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
  glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, w, h, 0, GL_RGBA, GL_UNSIGNED_BYTE,
               NULL);
  glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D,
                         *color_tex, 0);

  glGenRenderbuffers(1, depth_rb);
  glBindRenderbuffer(GL_RENDERBUFFER, *depth_rb);
  glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH_COMPONENT24, w, h);
  GLint depth_size = 0;
  glGetRenderbufferParameteriv(GL_RENDERBUFFER, GL_RENDERBUFFER_DEPTH_SIZE,
                               &depth_size);
  if (depth_bits)
    *depth_bits = (int)depth_size;
  glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT,
                            GL_RENDERBUFFER, *depth_rb);

  const GLenum status = glCheckFramebufferStatus(GL_FRAMEBUFFER);
  glBindFramebuffer(GL_FRAMEBUFFER, 0);
  if (status != GL_FRAMEBUFFER_COMPLETE) {
    fprintf(stderr, "Scene FBO incomplete: 0x%x\n", (unsigned)status);
    scene_fbo_destroy(fbo, color_tex, depth_rb, cur_w, cur_h);
    return false;
  }

  *cur_w = w;
  *cur_h = h;
  return true;
}

static void seed_solar_system(Sim *sim) {
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
  sim->merge_on_collision = false;

  typedef struct {
    const char *name;
    double mass_solar;
    double radius_km;
    uint32_t layer;
    double tilt_deg;
    double rotation_period_days;
    double x_km, y_km, z_km;
    double vx_km_s, vy_km_s, vz_km_s;
  } BodyInit;

  const BodyInit bodies[] = {
      {"Sun", 1.0, 695700.0, TEX_SUN, 7.25, 25.38, -1.068108951496322E+06,
       -4.177210908491462E+05, 3.086887010002915E+04, 9.305302656256911E-03,
       -1.283177282717393E-02, -1.631700118015769E-04},
      {"Mercury", 1.660e-7, 2439.7, TEX_MERCURY, 0.034, 58.646,
       -2.212073002393702E+07, -6.682435921338345E+07, -3.461577076477692E+06,
       3.666229234452722E+01, -1.230266984222893E+01, -4.368336206255391E+00},
      {"Venus", 2.447e-6, 6051.8, TEX_VENUS, 177.4, -243.025,
       -1.085736592234813E+08, -3.784241757371509E+06, 6.190088659339075E+06,
       8.984650886248794E-01, -3.517203951420625E+01, -5.320225928762774E-01},
      {"Earth", 3.003e-6, 6371.0, TEX_EARTH, 23.439281, 0.99726968,
       -2.627903751048988E+07, 1.445101984929515E+08, 3.025245352813601E+04,
       -2.983052803412253E+01, -5.220465675237847E+00, -1.014855999592612E-04},
      {"Mars", 3.213e-7, 3389.5, TEX_MARS, 25.19, 1.025957,
       2.069269460321208E+08, -3.560730804791640E+06, -5.147912373388751E+06,
       1.304308855632342E+00, 2.628158889664317E+01, 5.188465759107714E-01},
      {"Jupiter", 9.545e-4, 69911.0, TEX_JUPITER, 3.13, 0.41354,
       5.978410555886381E+08, 4.387048655696349E+08, -1.520164176015472E+07,
       -7.892632213479861E+00, 1.115034525890079E+01, 1.305100448596264E-01},
      {"Saturn", 2.857e-4, 58232.0, TEX_SATURN, 26.73, 0.44401,
       9.576382282218235E+08, 9.821474893679625E+08, -5.518978744215649E+07,
       -7.419580377753652E+00, 6.725982467906618E+00, 1.775011906748625E-01},
      {"Uranus", 4.366e-5, 25362.0, TEX_URANUS, 97.77, -0.71833,
       2.157706372184191E+09, -2.055243161071827E+09, -3.559278015686727E+07,
       4.646952926205451E+00, 4.614360059359629E+00, -4.301869182469398E-02},
      {"Neptune", 5.151e-5, 24622.0, TEX_NEPTUNE, 28.32, 0.6713,
       2.513785451779509E+09, -3.739265135509532E+09, 1.907027540535474E+07,
       4.475107938022004E+00, 3.062850546988970E+00, -1.667293921151841E-01},
  };

  for (size_t i = 0; i < sizeof(bodies) / sizeof(bodies[0]); i++) {
    const BodyInit *in = &bodies[i];
    const double r_au = au_from_km(in->radius_km);

    Body b = {0};
    snprintf(b.name, sizeof(b.name), "%s", in->name);
    b.x = au_from_km(in->x_km);
    b.y = au_from_km(in->y_km);
    b.z = au_from_km(in->z_km);
    b.vx = au_per_year_from_km_per_s(in->vx_km_s);
    b.vy = au_per_year_from_km_per_s(in->vy_km_s);
    b.vz = au_per_year_from_km_per_s(in->vz_km_s);
    b.mass = in->mass_solar;
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

typedef struct {
  float x;
  float y;
  float z;
  float radius;
  float r;
  float g;
  float b;
  uint32_t layer;
  float tilt;
  float spin_rate;
  float spin_phase;
} ParticleGPU;

typedef struct {
  float x;
  float y;
  float z;
  float outer;
  float inner_ratio;
  float nx;
  float ny;
  float nz;
} RingGPU;

enum {
  TRAIL_LEN = 2048,
};

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

typedef struct {
  char name_lc[32];
  GLuint tex;
  bool radial_strip;
  bool attempted;
} RingTexCache;

static void str_to_lower_ascii(char *dst, size_t cap, const char *src) {
  if (!cap)
    return;
  size_t i = 0;
  for (; src[i] && i + 1 < cap; i++) {
    char c = src[i];
    if (c >= 'A' && c <= 'Z')
      c = (char)(c - 'A' + 'a');
    dst[i] = c;
  }
  dst[i] = 0;
}

static GLuint ring_texture_for_body(const char *assets_dir, RingTexCache **cache,
                                    size_t *count, size_t *cap,
                                    const Body *b, bool *out_radial_strip) {
  char name_lc[32];
  str_to_lower_ascii(name_lc, sizeof(name_lc), b->name);
  if (!name_lc[0])
    return 0;

  for (size_t i = 0; i < *count; i++) {
    if (strcmp((*cache)[i].name_lc, name_lc) == 0) {
      if (out_radial_strip)
        *out_radial_strip = (*cache)[i].radial_strip;
      return (*cache)[i].tex;
    }
  }

  if (*count == *cap) {
    size_t next = *cap ? (*cap * 2) : 16;
    RingTexCache *nc = (RingTexCache *)realloc(*cache, next * sizeof(RingTexCache));
    if (!nc)
      return 0;
    *cache = nc;
    *cap = next;
  }

  RingTexCache *e = &(*cache)[(*count)++];
  memset(e, 0, sizeof(*e));
  strncpy(e->name_lc, name_lc, sizeof(e->name_lc) - 1);
  e->name_lc[sizeof(e->name_lc) - 1] = 0;
  e->attempted = true;

  char c0[64], c1[64], c2[64], c3[64], c4[64];
  snprintf(c0, sizeof(c0), "2k_%s_ring_alpha.png", name_lc);
  snprintf(c1, sizeof(c1), "2k_%s_ring_alpha.jpg", name_lc);
  snprintf(c2, sizeof(c2), "%s_ring_alpha.png", name_lc);
  snprintf(c3, sizeof(c3), "%s_ring.png", name_lc);
  snprintf(c4, sizeof(c4), "2k_%s_ring.png", name_lc);
  const char *cands[] = {c0, c1, c2, c3, c4, NULL};

  char label[80];
  snprintf(label, sizeof(label), "%s rings", b->name);
  int w = 0, h = 0;
  e->tex = textures_try_load_rgba2d(assets_dir, label, cands, &w, &h);
  e->radial_strip = (e->tex && h > 0 && w >= h * 2);
  if (out_radial_strip)
    *out_radial_strip = e->radial_strip;
  return e->tex;
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
  GLuint scene_fbo = 0;
  GLuint scene_color = 0;
  GLuint scene_depth = 0;
  int scene_w = 0;
  int scene_h = 0;
  int scene_depth_bits = 0;

  // Particle program
  GLuint pvs = glutil_compile_shader(GL_VERTEX_SHADER, kParticlesVS);
  GLuint pfs = glutil_compile_shader(GL_FRAGMENT_SHADER, kParticlesFS);
  if (!pvs || !pfs)
    return 1;
  GLuint particle_prog = glutil_link_program(pvs, pfs);
  glDeleteShader(pvs);
  glDeleteShader(pfs);
  if (!particle_prog)
    return 1;

  GLint puWorldToClip = glGetUniformLocation(particle_prog, "uWorldToClip");
  GLint puCamRight = glGetUniformLocation(particle_prog, "uCamRight");
  GLint puCamUp = glGetUniformLocation(particle_prog, "uCamUp");
  GLint puTex = glGetUniformLocation(particle_prog, "uTex");
  GLint puTime = glGetUniformLocation(particle_prog, "uTime");

  // Background program
  GLuint bvs = glutil_compile_shader(GL_VERTEX_SHADER, kBgVS);
  GLuint bfs = glutil_compile_shader(GL_FRAGMENT_SHADER, kBgFS);
  if (!bvs || !bfs)
    return 1;
  GLuint bg_prog = glutil_link_program(bvs, bfs);
  glDeleteShader(bvs);
  glDeleteShader(bfs);
  if (!bg_prog)
    return 1;
  GLint buEnv = glGetUniformLocation(bg_prog, "uEnv");
  GLint buHasEnv = glGetUniformLocation(bg_prog, "uHasEnv");
  GLint buInvViewRot = glGetUniformLocation(bg_prog, "uInvViewRot");
  GLint buAspect = glGetUniformLocation(bg_prog, "uAspect");
  GLint buExposure = glGetUniformLocation(bg_prog, "uExposure");

  // Text program
  GLuint tvs = glutil_compile_shader(GL_VERTEX_SHADER, kTextVS);
  GLuint tfs = glutil_compile_shader(GL_FRAGMENT_SHADER, kTextFS);
  if (!tvs || !tfs)
    return 1;
  GLuint text_prog = glutil_link_program(tvs, tfs);
  glDeleteShader(tvs);
  glDeleteShader(tfs);
  if (!text_prog)
    return 1;
  GLint tuViewport = glGetUniformLocation(text_prog, "uViewport");
  GLint tuOffsetPx = glGetUniformLocation(text_prog, "uOffsetPx");
  GLint tuFont = glGetUniformLocation(text_prog, "uFont");
  GLint tuColor = glGetUniformLocation(text_prog, "uColor");

  // Ring program
  GLuint rvs = glutil_compile_shader(GL_VERTEX_SHADER, kRingsVS);
  GLuint rfs = glutil_compile_shader(GL_FRAGMENT_SHADER, kRingsFS);
  if (!rvs || !rfs)
    return 1;
  GLuint ring_prog = glutil_link_program(rvs, rfs);
  glDeleteShader(rvs);
  glDeleteShader(rfs);
  if (!ring_prog)
    return 1;
  GLint ruWorldToClip = glGetUniformLocation(ring_prog, "uWorldToClip");
  GLint ruRingTex = glGetUniformLocation(ring_prog, "uRingTex");
  GLint ruRadialStrip = glGetUniformLocation(ring_prog, "uRadialStrip");

  // Line program (trails, vectors)
  GLuint lvs = glutil_compile_shader(GL_VERTEX_SHADER, kLinesVS);
  GLuint lfs = glutil_compile_shader(GL_FRAGMENT_SHADER, kLinesFS);
  if (!lvs || !lfs)
    return 1;
  GLuint line_prog = glutil_link_program(lvs, lfs);
  glDeleteShader(lvs);
  glDeleteShader(lfs);
  if (!line_prog)
    return 1;
  GLint luWorldToClip = glGetUniformLocation(line_prog, "uWorldToClip");

  GLuint particle_vao = 0, particle_vbo = 0;
  glGenVertexArrays(1, &particle_vao);
  glGenBuffers(1, &particle_vbo);
  glBindVertexArray(particle_vao);
  glBindBuffer(GL_ARRAY_BUFFER, particle_vbo);
  glBufferData(GL_ARRAY_BUFFER, 1, NULL, GL_STREAM_DRAW);
  glEnableVertexAttribArray(0);
  glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(ParticleGPU),
                        (void *)offsetof(ParticleGPU, x));
  glEnableVertexAttribArray(1);
  glVertexAttribPointer(1, 1, GL_FLOAT, GL_FALSE, sizeof(ParticleGPU),
                        (void *)offsetof(ParticleGPU, radius));
  glEnableVertexAttribArray(2);
  glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, sizeof(ParticleGPU),
                        (void *)offsetof(ParticleGPU, r));
  glEnableVertexAttribArray(3);
  glVertexAttribIPointer(3, 1, GL_UNSIGNED_INT, sizeof(ParticleGPU),
                         (void *)offsetof(ParticleGPU, layer));
  glEnableVertexAttribArray(4);
  glVertexAttribPointer(4, 3, GL_FLOAT, GL_FALSE, sizeof(ParticleGPU),
                        (void *)offsetof(ParticleGPU, tilt));
  glVertexAttribDivisor(0, 1);
  glVertexAttribDivisor(1, 1);
  glVertexAttribDivisor(2, 1);
  glVertexAttribDivisor(3, 1);
  glVertexAttribDivisor(4, 1);
  glBindVertexArray(0);

  GLuint line_vao = 0, line_vbo = 0;
  glGenVertexArrays(1, &line_vao);
  glGenBuffers(1, &line_vbo);
  glBindVertexArray(line_vao);
  glBindBuffer(GL_ARRAY_BUFFER, line_vbo);
  glBufferData(GL_ARRAY_BUFFER, 1, NULL, GL_STREAM_DRAW);
  glEnableVertexAttribArray(0);
  glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(LineVertex),
                        (void *)offsetof(LineVertex, x));
  glEnableVertexAttribArray(1);
  glVertexAttribPointer(1, 4, GL_FLOAT, GL_FALSE, sizeof(LineVertex),
                        (void *)offsetof(LineVertex, r));
  glBindVertexArray(0);

  GLuint bg_vao = 0;
  glGenVertexArrays(1, &bg_vao);
  glBindVertexArray(0);

  GLuint text_vao = 0, text_vbo = 0;
  glGenVertexArrays(1, &text_vao);
  glGenBuffers(1, &text_vbo);
  glBindVertexArray(text_vao);
  glBindBuffer(GL_ARRAY_BUFFER, text_vbo);
  glBufferData(GL_ARRAY_BUFFER, 1, NULL, GL_STREAM_DRAW);
  glEnableVertexAttribArray(0);
  glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, sizeof(TextVertex),
                        (void *)offsetof(TextVertex, x));
  glEnableVertexAttribArray(1);
  glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, sizeof(TextVertex),
                        (void *)offsetof(TextVertex, u));
  glBindVertexArray(0);

  GLuint ring_vao = 0, ring_vbo = 0;
  glGenVertexArrays(1, &ring_vao);
  glGenBuffers(1, &ring_vbo);
  glBindVertexArray(ring_vao);
  glBindBuffer(GL_ARRAY_BUFFER, ring_vbo);
  glBufferData(GL_ARRAY_BUFFER, 1, NULL, GL_STREAM_DRAW);
  glEnableVertexAttribArray(0);
  glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(RingGPU), (void *)offsetof(RingGPU, x));
  glEnableVertexAttribArray(1);
  glVertexAttribPointer(1, 1, GL_FLOAT, GL_FALSE, sizeof(RingGPU), (void *)offsetof(RingGPU, outer));
  glEnableVertexAttribArray(2);
  glVertexAttribPointer(2, 1, GL_FLOAT, GL_FALSE, sizeof(RingGPU), (void *)offsetof(RingGPU, inner_ratio));
  glEnableVertexAttribArray(3);
  glVertexAttribPointer(3, 3, GL_FLOAT, GL_FALSE, sizeof(RingGPU), (void *)offsetof(RingGPU, nx));
  glVertexAttribDivisor(0, 1);
  glVertexAttribDivisor(1, 1);
  glVertexAttribDivisor(2, 1);
  glVertexAttribDivisor(3, 1);
  glBindVertexArray(0);

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

  glUseProgram(particle_prog);
  glUniform1i(puTex, 0);
  glUniform1f(puTime, 0.0f);
  glUseProgram(0);

  // Ensure unit 0 has the correct texture target bound for sampler2DArray.
  glActiveTexture(GL_TEXTURE0);
  glBindTexture(GL_TEXTURE_2D_ARRAY, tex_array);

  char *hdr_path = path_join2(assets_dir, "HDR_blue_nebulae-1.hdr");
  GLuint env_tex = hdr_path ? textures_load_hdr_equirect(hdr_path) : 0;
  free(hdr_path);
  glUseProgram(bg_prog);
  glUniform1i(buEnv, 1);
  glUniform1i(buHasEnv, env_tex ? 1 : 0);
  glUniform1f(buExposure, 1.2f);
  glUseProgram(0);

  GLuint font_tex = create_font_texture();
  if (!font_tex) {
    fprintf(stderr, "Failed to create font texture\n");
    return 1;
  }
  glUseProgram(text_prog);
  glUniform1i(tuFont, 2);
  glUseProgram(0);

  glUseProgram(ring_prog);
  glUniform1i(ruRingTex, 3);
  glUniform1i(ruRadialStrip, 0);
  glUseProgram(0);

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
  seed_solar_system(&sim);

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
                .pitch = -1.24,
                .zoom_world_h = 4.0};

  bool running = true;
  bool paused = false;
  bool show_help = false;
  bool show_labels = false;
  bool realistic_scale = false;

  bool panning = false;
  int pan_start_x = 0, pan_start_y = 0;
  double pan_start_cx = 0.0, pan_start_cy = 0.0, pan_start_cz = 0.0;

  bool rotate_down = false;
  bool rotating_cam = false;
  int rotate_start_x = 0, rotate_start_y = 0;
  double rotate_start_yaw = 0.0, rotate_start_pitch = 0.0;

  bool spawning = false;
  int spawn_start_sx = 0, spawn_start_sy = 0;
  double spawn_start_wx = 0.0, spawn_start_wy = 0.0, spawn_start_wz = 0.0;
  double spawn_mass = 100.0;
  double spawn_vel_scale = 1.0;

  uint32_t selected_id = 0;
  bool follow_selected = false;
  bool show_trails = true;
  bool show_vectors = true;
  double vector_scale = 0.02;
  double render_radius_scale = 2.5;

  Trail *trails = NULL;
  size_t trail_count = 0;
  size_t trail_cap = 0;

  ParticleGPU *particle_cpu = NULL;
  size_t particle_cpu_cap = 0;

  RingGPU *ring_cpu = NULL;
  GLuint *ring_tex_cpu = NULL;
  uint8_t *ring_radial_cpu = NULL;
  size_t ring_cpu_cap = 0;
  size_t ring_count = 0;

  RingTexCache *ring_cache = NULL;
  size_t ring_cache_count = 0;
  size_t ring_cache_cap = 0;

  LineVertex *line_cpu = NULL;
  size_t line_cpu_cap = 0;

  TextVertex *text_cpu = NULL;
  size_t text_cpu_cap = 0;

  // Interpreted as "simulation time units per real second".
  // For the solar-system preset, units are years.
  double time_scale = 0.03;
  const double fixed_dt = 1.0e-4;
  double accumulator = 0.0;
  double sim_time = 0.0;

  uint64_t prev_ticks = (uint64_t)SDL_GetPerformanceCounter();
  const uint64_t freq = (uint64_t)SDL_GetPerformanceFrequency();

  uint64_t title_ticks = prev_ticks;
  int frames_since_title = 0;

  print_controls();

  // Select the Sun by default (first body in solar system preset).
  if (sim.count)
    selected_id = sim.bodies[0].id;

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
          fprintf(stderr, "realistic_scale: %s\n", realistic_scale ? "on" : "off");
          break;
        case SDLK_ESCAPE:
          running = false;
          break;
        case SDLK_SPACE:
          paused = !paused;
          break;
        case SDLK_r:
          seed_solar_system(&sim);
          preset = PRESET_SOLAR;
          sim_time = 0.0;
          cam.cx = 0.0;
          cam.cy = 0.0;
          cam.cz = 0.0;
          cam.yaw = 0.0;
          cam.pitch = -1.24, cam.zoom_world_h = 4.0;
          time_scale = 0.03;
          selected_id = sim.count ? sim.bodies[0].id : 0;
          follow_selected = false;
          break;
        case SDLK_1:
          seed_two_body_orbit(&sim);
          preset = PRESET_TWO_BODY;
          sim_time = 0.0;
          cam.cx = 0.0;
          cam.cy = 0.0;
          cam.cz = 0.0;
          cam.yaw = 0.0;
          cam.pitch = 0.0;
          cam.zoom_world_h = 40.0;
          time_scale = 1.0;
          selected_id = 0;
          follow_selected = false;
          break;
        case SDLK_2:
          seed_disk_galaxy(&sim, 1500, 140.0, 2.0e6);
          preset = PRESET_GALAXY;
          sim_time = 0.0;
          cam.cx = 0.0;
          cam.cy = 0.0;
          cam.cz = 0.0;
          cam.yaw = 0.0;
          cam.pitch = 0.0;
          cam.zoom_world_h = 320.0;
          time_scale = 1.0;
          selected_id = 0;
          follow_selected = false;
          break;
        case SDLK_3:
          seed_solar_system(&sim);
          preset = PRESET_SOLAR;
          sim_time = 0.0;
          cam.cx = 0.0;
          cam.cy = 0.0;
          cam.cz = 0.0;
          cam.yaw = 0.0;
          cam.pitch = -1.24, cam.zoom_world_h = 4.0;
          time_scale = 0.03;
          selected_id = sim.count ? sim.bodies[0].id : 0;
          follow_selected = false;
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
        const double z = (e.wheel.y > 0) ? (1.0 / 1.15) : 1.15;
        cam.zoom_world_h = clampd(cam.zoom_world_h * z, 1e-3, 1e9);
      }

      if (e.type == SDL_MOUSEBUTTONDOWN) {
        if (e.button.button == SDL_BUTTON_MIDDLE) {
          int dx = 0, dy = 0;
          window_to_drawable(win, e.button.x, e.button.y, &dx, &dy);
          panning = true;
          pan_start_x = dx;
          pan_start_y = dy;
          pan_start_cx = cam.cx;
          pan_start_cy = cam.cy;
          pan_start_cz = cam.cz;
        } else if (e.button.button == SDL_BUTTON_LEFT) {
          int dx = 0, dy = 0;
          window_to_drawable(win, e.button.x, e.button.y, &dx, &dy);
          rotate_down = true;
          rotating_cam = false;
          rotate_start_x = dx;
          rotate_start_y = dy;
          rotate_start_yaw = cam.yaw;
          rotate_start_pitch = cam.pitch;
        } else if (e.button.button == SDL_BUTTON_RIGHT) {
          int dx = 0, dy = 0;
          window_to_drawable(win, e.button.x, e.button.y, &dx, &dy);
          spawning = true;
          spawn_start_sx = dx;
          spawn_start_sy = dy;
          int dw = 1, dh = 1;
          SDL_GL_GetDrawableSize(win, &dw, &dh);
          screen_to_world_on_z0(&cam, dw, dh, spawn_start_sx, spawn_start_sy,
                                &spawn_start_wx, &spawn_start_wy,
                                &spawn_start_wz);

          const uint16_t mod = SDL_GetModState();
          spawn_vel_scale = (mod & KMOD_SHIFT) ? 5.0 : 1.0;
        }
      }

      if (e.type == SDL_MOUSEBUTTONUP) {
        if (e.button.button == SDL_BUTTON_MIDDLE) {
          panning = false;
        } else if (e.button.button == SDL_BUTTON_RIGHT) {
          if (spawning) {
            int dw = 1, dh = 1;
            SDL_GL_GetDrawableSize(win, &dw, &dh);
            double end_wx = 0.0, end_wy = 0.0, end_wz = 0.0;
            int dx = 0, dy = 0;
            window_to_drawable(win, e.button.x, e.button.y, &dx, &dy);
            screen_to_world_on_z0(&cam, dw, dh, dx, dy, &end_wx, &end_wy,
                                  &end_wz);

            Body b = {0};
            b.x = spawn_start_wx;
            b.y = spawn_start_wy;
            b.z = spawn_start_wz;
            const uint16_t mod = SDL_GetModState();
            if (mod & KMOD_CTRL) {
              const Body *center = selected_id
                                       ? find_body_by_id(&sim, selected_id)
                                       : find_heaviest_body(&sim);
              if (center) {
                const double rx = b.x - center->x;
                const double ry = b.y - center->y;
                const double rz = b.z - center->z;
                const double r = sqrt(rx * rx + ry * ry + rz * rz);
                if (r > 1e-6) {
                  const double invr = 1.0 / r;
                  // Tangent direction in the z=0 plane.
                  double tx = -ry * invr;
                  double ty = rx * invr;
                  if (mod & KMOD_ALT) {
                    tx = -tx;
                    ty = -ty;
                  }
                  const double v_circ = sqrt(sim.G * center->mass / r);
                  b.vx = center->vx + tx * v_circ;
                  b.vy = center->vy + ty * v_circ;
                  b.vz = center->vz;
                }
              }
            } else {
              b.vx = (end_wx - spawn_start_wx) * spawn_vel_scale;
              b.vy = (end_wy - spawn_start_wy) * spawn_vel_scale;
              b.vz = (end_wz - spawn_start_wz) * spawn_vel_scale;
            }
            b.mass = spawn_mass;
            (void)sim_add_body(&sim, b);
          }
          spawning = false;
        } else if (e.button.button == SDL_BUTTON_LEFT) {
          // If we didn't rotate, treat this as a click-selection.
          if (rotate_down && !rotating_cam) {
            int dw = 1, dh = 1;
            SDL_GL_GetDrawableSize(win, &dw, &dh);

            int mdx = 0, mdy = 0;
            window_to_drawable(win, e.button.x, e.button.y, &mdx, &mdy);

            Mat4 world_to_clip = camera_world_to_clip(&cam, dw, dh);
            double best_d2_px = 0.0;
            uint32_t best_id = 0;
            for (size_t i = 0; i < sim.count; i++) {
              const Body *b = &sim.bodies[i];
              float clip[4];
              mat4_mul_vec4(&world_to_clip, (float)b->x, (float)b->y,
                            (float)b->z, 1.0f, clip);
              if (fabsf(clip[3]) < 1e-8f)
                continue;
              const float ndc_x = clip[0] / clip[3];
              const float ndc_y = clip[1] / clip[3];
              const float sx = (ndc_x * 0.5f + 0.5f) * (float)dw;
              const float sy = (0.5f - ndc_y * 0.5f) * (float)dh;

              const float dx = sx - (float)mdx;
              const float dy = sy - (float)mdy;
              const double d2 =
                  (double)dx * (double)dx + (double)dy * (double)dy;
              // Cap pick radius so the Sun doesn't "eat" all clicks.
              const double radius_px_raw =
                  body_visual_radius_world(b, realistic_scale,
                                           render_radius_scale) *
                  ((double)dh / cam.zoom_world_h);
              const double pick_r_px = clampd(radius_px_raw * 1.4, 10.0, 60.0);
              if (d2 > pick_r_px * pick_r_px)
                continue;
              if (!best_id || d2 < best_d2_px) {
                best_d2_px = d2;
                best_id = b->id;
              }
            }

            // Require being reasonably close to something.
            if (best_id && best_d2_px <= 60.0 * 60.0) {
              selected_id = best_id;
            } else {
              selected_id = 0;
            }
            if (selected_id) {
              const Body *b = find_body_by_id(&sim, selected_id);
              fprintf(stderr, "selected: %s (id=%u mass=%.6g)\n",
                      b ? b->name : "(unknown)", (unsigned)selected_id,
                      b ? b->mass : 0.0);
            } else {
              fprintf(stderr, "selection cleared\n");
              follow_selected = false;
            }
          }
          rotate_down = false;
          rotating_cam = false;
        }
      }

      if (e.type == SDL_MOUSEMOTION) {
        if (panning) {
          int dw = 1, dh = 1;
          SDL_GL_GetDrawableSize(win, &dw, &dh);
          int mdx = 0, mdy = 0;
          window_to_drawable(win, e.motion.x, e.motion.y, &mdx, &mdy);
          const int dx = mdx - pan_start_x;
          const int dy = mdy - pan_start_y;
          const double aspect = (dh > 0) ? ((double)dw / (double)dh) : 1.0;
          const double world_per_px_y = cam.zoom_world_h / (double)dh;
          const double world_per_px_x = cam.zoom_world_h * aspect / (double)dw;

          const float vx = (float)(-(double)dx * world_per_px_x);
          const float vy = (float)((double)dy * world_per_px_y);
          float wx, wy, wz;
          view_to_world_dir(&cam, vx, vy, 0.0f, &wx, &wy, &wz);
          cam.cx = pan_start_cx + (double)wx;
          cam.cy = pan_start_cy + (double)wy;
          cam.cz = pan_start_cz + (double)wz;
        }

        if (rotate_down) {
          int mdx = 0, mdy = 0;
          window_to_drawable(win, e.motion.x, e.motion.y, &mdx, &mdy);
          const int dx = mdx - rotate_start_x;
          const int dy = mdy - rotate_start_y;
          if (!rotating_cam && (abs(dx) + abs(dy) > 3))
            rotating_cam = true;
          if (rotating_cam) {
            const double sens = 0.005;
            cam.yaw = rotate_start_yaw + (double)dx * sens;
            cam.pitch =
                clampd(rotate_start_pitch + (double)dy * sens, -1.55, 1.55);
          }
        }
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
      if (!scene_fbo_ensure(&scene_fbo, &scene_color, &scene_depth, &scene_w,
                            &scene_h, &scene_depth_bits, dw, dh)) {
        use_offscreen_scene_fbo = false;
      }
    }

    glBindFramebuffer(GL_FRAMEBUFFER, use_offscreen_scene_fbo ? scene_fbo : 0);
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
    {
      const float cy = cosf((float)-cam.yaw);
      const float sy = sinf((float)-cam.yaw);
      const float cp = cosf((float)-cam.pitch);
      const float sp = sinf((float)-cam.pitch);
      // uInvViewRot = Rz(-yaw) * Rx(-pitch)
      float invR[9];
      invR[0] = cy;
      invR[1] = sy;
      invR[2] = 0.0f;
      invR[3] = -sy * cp;
      invR[4] = cy * cp;
      invR[5] = sp;
      invR[6] = sy * sp;
      invR[7] = -cy * sp;
      invR[8] = cp;

      glBindVertexArray(bg_vao);
      glDisable(GL_DEPTH_TEST);
      glDepthMask(GL_FALSE);
      glDisable(GL_BLEND);
      glUseProgram(bg_prog);
      glUniformMatrix3fv(buInvViewRot, 1, GL_FALSE, invR);
      glUniform1f(buAspect, (dh > 0) ? ((float)dw / (float)dh) : 1.0f);
      glActiveTexture(GL_TEXTURE1);
      glBindTexture(GL_TEXTURE_2D, env_tex);
      glDrawArrays(GL_TRIANGLES, 0, 3);
      glUseProgram(0);
      glBindVertexArray(0);
      glDepthMask(GL_TRUE);
      glEnable(GL_DEPTH_TEST);
      glEnable(GL_BLEND);
    }

    // Build line vertex buffer (trails + vectors + spawn preview).
    size_t line_count = 0;
    if (show_trails) {
      for (size_t ti = 0; ti < trail_count; ti++) {
        const Trail *t = &trails[ti];
        const Body *b = find_body_by_id(&sim, t->id);
        if (!b || t->count < 2)
          continue;

        const size_t start = (t->head + TRAIL_LEN - t->count) % TRAIL_LEN;
        for (size_t k = 0; k + 1 < t->count; k++) {
          const size_t i0 = (start + k) % TRAIL_LEN;
          const size_t i1 = (start + k + 1) % TRAIL_LEN;
          const float age0 = (float)k / (float)(t->count - 1);
          const float age1 = (float)(k + 1) / (float)(t->count - 1);
          const float a0 = 0.55f * age0 * age0;
          const float a1 = 0.55f * age1 * age1;

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

          line_cpu[line_count++] = (LineVertex){.x = t->x[i0],
                                                .y = t->y[i0],
                                                .z = t->z[i0],
                                                .r = b->r,
                                                .g = b->g,
                                                .b = b->b,
                                                .a = a0};
          line_cpu[line_count++] = (LineVertex){.x = t->x[i1],
                                                .y = t->y[i1],
                                                .z = t->z[i1],
                                                .r = b->r,
                                                .g = b->g,
                                                .b = b->b,
                                                .a = a1};
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

    if (spawning) {
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
        line_cpu[line_count++] = (LineVertex){.x = (float)spawn_start_wx,
                                              .y = (float)spawn_start_wy,
                                              .z = (float)spawn_start_wz,
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
      glUseProgram(line_prog);
      glUniformMatrix4fv(luWorldToClip, 1, GL_FALSE, world_to_clip.m);
      glBindVertexArray(line_vao);
      glBindBuffer(GL_ARRAY_BUFFER, line_vbo);
      glBufferData(GL_ARRAY_BUFFER,
                   (GLsizeiptr)(line_count * sizeof(LineVertex)), line_cpu,
                   GL_STREAM_DRAW);
      glDrawArrays(GL_LINES, 0, (GLsizei)line_count);
      glBindVertexArray(0);
      glUseProgram(0);
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

    // Build ring buffer (only bodies with a ring texture).
    ring_count = 0;
    for (size_t i = 0; i < sim.count; i++) {
      const Body *b = &sim.bodies[i];
      bool radial_strip = false;
      GLuint ring_tex = ring_texture_for_body(assets_dir, &ring_cache, &ring_cache_count, &ring_cache_cap, b, &radial_strip);
      if (!ring_tex)
        continue;

      if (ring_count + 1 > ring_cpu_cap) {
        size_t next = ring_cpu_cap ? (ring_cpu_cap * 2) : 16;
        while (next < ring_count + 1)
          next *= 2;
        RingGPU *nr = (RingGPU *)malloc(next * sizeof(RingGPU));
        GLuint *nt = (GLuint *)malloc(next * sizeof(GLuint));
        uint8_t *nm = (uint8_t *)malloc(next * sizeof(uint8_t));
        if (!nr || !nt || !nm) {
          free(nr);
          free(nt);
          free(nm);
          break;
        }
        if (ring_cpu)
          memcpy(nr, ring_cpu, ring_count * sizeof(RingGPU));
        if (ring_tex_cpu)
          memcpy(nt, ring_tex_cpu, ring_count * sizeof(GLuint));
        if (ring_radial_cpu)
          memcpy(nm, ring_radial_cpu, ring_count * sizeof(uint8_t));
        free(ring_cpu);
        free(ring_tex_cpu);
        free(ring_radial_cpu);
        ring_cpu = nr;
        ring_tex_cpu = nt;
        ring_radial_cpu = nm;
        ring_cpu_cap = next;
      }

      const double pr =
          body_visual_radius_world(b, realistic_scale, render_radius_scale);
      const double outer = pr * 2.5;
      const double inner = pr * 1.2;

      // Use the same tilt convention as the planet shader.
      const float tilt = b->tilt_rad;
      float nx = 0.0f;
      float ny = cosf(tilt);
      float nz = sinf(tilt);
      const float nlen = sqrtf(nx * nx + ny * ny + nz * nz);
      if (nlen > 1e-8f) {
        nx /= nlen;
        ny /= nlen;
        nz /= nlen;
      } else {
        nx = 0.0f;
        ny = 1.0f;
        nz = 0.0f;
      }

      ring_cpu[ring_count++] = (RingGPU){.x = (float)b->x,
                                         .y = (float)b->y,
                                         .z = (float)b->z,
                                         .outer = (float)outer,
                                         .inner_ratio = (float)(inner / outer),
                                         .nx = nx,
                                         .ny = ny,
                                         .nz = nz};
      ring_tex_cpu[ring_count - 1] = ring_tex;
      ring_radial_cpu[ring_count - 1] = radial_strip ? 1u : 0u;
    }

    glUseProgram(particle_prog);
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D_ARRAY, tex_array);
    glUniform1f(puTime, (float)sim_time);
    glUniformMatrix4fv(puWorldToClip, 1, GL_FALSE, world_to_clip.m);
    glUniform3f(puCamRight, cam_right_x, cam_right_y, cam_right_z);
    glUniform3f(puCamUp, cam_up_x, cam_up_y, cam_up_z);
    glBindVertexArray(particle_vao);
    glBindBuffer(GL_ARRAY_BUFFER, particle_vbo);
    glBufferData(GL_ARRAY_BUFFER, (GLsizeiptr)(sim.count * sizeof(ParticleGPU)),
                 particle_cpu, GL_STREAM_DRAW);

    // Depth pre-pass so nearer bodies properly occlude farther ones.
    glDisable(GL_BLEND);
    glColorMask(GL_FALSE, GL_FALSE, GL_FALSE, GL_FALSE);
    glDepthMask(GL_TRUE);
    glDepthFunc(GL_LESS);
    glDrawArraysInstanced(GL_TRIANGLE_STRIP, 0, 4, (GLsizei)sim.count);

    // Color pass: draw only the visible fragments.
    glEnable(GL_BLEND);
    glColorMask(GL_TRUE, GL_TRUE, GL_TRUE, GL_TRUE);
    glDepthMask(GL_FALSE);
    glDepthFunc(GL_EQUAL);
    glDrawArraysInstanced(GL_TRIANGLE_STRIP, 0, 4, (GLsizei)sim.count);

    // Restore defaults for other draws.
    glDepthFunc(GL_LESS);
    glDepthMask(GL_TRUE);
    glBindVertexArray(0);
    glUseProgram(0);

    if (ring_count) {
      glUseProgram(ring_prog);
      glUniformMatrix4fv(ruWorldToClip, 1, GL_FALSE, world_to_clip.m);
      glActiveTexture(GL_TEXTURE3);

      glBindVertexArray(ring_vao);
      glBindBuffer(GL_ARRAY_BUFFER, ring_vbo);

      glEnable(GL_BLEND);
      glDepthMask(GL_FALSE);
      glDepthFunc(GL_LEQUAL);
      for (size_t i = 0; i < ring_count; i++) {
        glBindTexture(GL_TEXTURE_2D, ring_tex_cpu[i]);
        glUniform1i(ruRadialStrip, ring_radial_cpu[i] ? 1 : 0);
        glBufferData(GL_ARRAY_BUFFER, (GLsizeiptr)sizeof(RingGPU), &ring_cpu[i],
                     GL_STREAM_DRAW);
        glDrawArraysInstanced(GL_TRIANGLE_STRIP, 0, 4, 1);
      }
      glDepthFunc(GL_LESS);
      glDepthMask(GL_TRUE);

      glBindVertexArray(0);
      glUseProgram(0);
    }

    if (use_offscreen_scene_fbo) {
      glBindFramebuffer(GL_READ_FRAMEBUFFER, scene_fbo);
      glBindFramebuffer(GL_DRAW_FRAMEBUFFER, 0);
      glBlitFramebuffer(0, 0, dw, dh, 0, 0, dw, dh, GL_COLOR_BUFFER_BIT,
                        GL_NEAREST);
      glBindFramebuffer(GL_FRAMEBUFFER, 0);
      glViewport(0, 0, dw, dh);
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

    // Simulation date/time (only meaningful for the solar-system ephemeris).
    if (preset == PRESET_SOLAR) {
      glDisable(GL_DEPTH_TEST);
      glDepthMask(GL_FALSE);
      const double jd = solar_epoch_jd + sim_time * 365.25;
      int yy = 0, mo = 0, dd = 0, hh = 0, mm = 0, ss = 0;
      jd_to_gregorian_utc(jd, &yy, &mo, &dd, &hh, &mm, &ss);
      char dt[64];
      snprintf(dt, sizeof(dt), "%04d-%02d-%02d %02d:%02d:%02d", yy, mo, dd, hh,
               mm, ss);

      const float scale = 3.0f;
      const float text_w = (float)strlen(dt) * 8.0f * scale;
      const float x = (float)clampd(((double)dw - (double)text_w) * 0.5, 4.0,
                                    (double)dw - 4.0 - text_w);
      const float y = 10.0f;

      size_t tv_count = 0;
      text_append(&text_cpu, &tv_count, &text_cpu_cap, x, y, scale, dt);
      if (tv_count) {
        glBindVertexArray(text_vao);
        glBindBuffer(GL_ARRAY_BUFFER, text_vbo);
        glBufferData(GL_ARRAY_BUFFER,
                     (GLsizeiptr)(tv_count * sizeof(TextVertex)), text_cpu,
                     GL_STREAM_DRAW);
        glUseProgram(text_prog);
        glUniform2f(tuViewport, (float)dw, (float)dh);
        glActiveTexture(GL_TEXTURE2);
        glBindTexture(GL_TEXTURE_2D, font_tex);

        glUniform2f(tuOffsetPx, 2.0f, 2.0f);
        glUniform4f(tuColor, 0.0f, 0.0f, 0.0f, 0.65f);
        glDrawArrays(GL_TRIANGLES, 0, (GLsizei)tv_count);

        glUniform2f(tuOffsetPx, 0.0f, 0.0f);
        glUniform4f(tuColor, 1.0f, 1.0f, 1.0f, 0.95f);
        glDrawArrays(GL_TRIANGLES, 0, (GLsizei)tv_count);

        glUseProgram(0);
        glBindVertexArray(0);
      }
      glDepthMask(GL_TRUE);
      glEnable(GL_DEPTH_TEST);
    }

    if (selected_id) {
      const Body *sel = find_body_by_id(&sim, selected_id);
      if (sel) {
        glDisable(GL_DEPTH_TEST);
        glDepthMask(GL_FALSE);
        char hud[256];
        snprintf(hud, sizeof(hud), "Selected: %s", sel->name);

        size_t tv_count = 0;
        text_append(&text_cpu, &tv_count, &text_cpu_cap, 12.0f, 12.0f, 1.0f,
                    hud);
        if (tv_count) {
          glBindVertexArray(text_vao);
          glBindBuffer(GL_ARRAY_BUFFER, text_vbo);
          glBufferData(GL_ARRAY_BUFFER,
                       (GLsizeiptr)(tv_count * sizeof(TextVertex)), text_cpu,
                       GL_STREAM_DRAW);
          glUseProgram(text_prog);
          glUniform2f(tuViewport, (float)dw, (float)dh);
          glActiveTexture(GL_TEXTURE2);
          glBindTexture(GL_TEXTURE_2D, font_tex);
          glUniform2f(tuOffsetPx, 1.0f, 1.0f);
          glUniform4f(tuColor, 0.0f, 0.0f, 0.0f, 0.55f);
          glDrawArrays(GL_TRIANGLES, 0, (GLsizei)tv_count);
          glUniform2f(tuOffsetPx, 0.0f, 0.0f);
          glUniform4f(tuColor, 1.0f, 1.0f, 1.0f, 0.9f);
          glDrawArrays(GL_TRIANGLES, 0, (GLsizei)tv_count);
          glUseProgram(0);
          glBindVertexArray(0);
        }

        // Label above the selected body.
        {
          float sx = 0.0f, sy = 0.0f;
          if (world_to_screen_px(&world_to_clip, dw, dh, (float)sel->x,
                                 (float)sel->y, (float)sel->z, &sx, &sy)) {
            const float scale = 1.5f;
            float radius_px =
                (float)(body_visual_radius_world(sel, realistic_scale,
                                                 render_radius_scale) *
                        (double)dh / cam.zoom_world_h);
            // Prevent the label from flying away when zoomed in.
            radius_px = (float)clampd(radius_px, 6.0, 70.0);
            const float text_w = (float)strlen(sel->name) * 8.0f * scale;
            const float x = clampd((double)(sx - text_w * 0.5f), 4.0,
                                   (double)dw - 4.0 - text_w);
            const float y =
                (float)clampd((double)(sy - radius_px - 18.0f * scale), 4.0,
                              (double)dh - 20.0);

            size_t lcount = 0;
            text_append(&text_cpu, &lcount, &text_cpu_cap, x, y, scale,
                        sel->name);
            if (lcount) {
              glBindVertexArray(text_vao);
              glBindBuffer(GL_ARRAY_BUFFER, text_vbo);
              glBufferData(GL_ARRAY_BUFFER,
                           (GLsizeiptr)(lcount * sizeof(TextVertex)), text_cpu,
                           GL_STREAM_DRAW);
              glUseProgram(text_prog);
              glUniform2f(tuViewport, (float)dw, (float)dh);
              glActiveTexture(GL_TEXTURE2);
              glBindTexture(GL_TEXTURE_2D, font_tex);

              glUniform2f(tuOffsetPx, 1.0f, 1.0f);
              glUniform4f(tuColor, 0.0f, 0.0f, 0.0f, 0.6f);
              glDrawArrays(GL_TRIANGLES, 0, (GLsizei)lcount);

              glUniform2f(tuOffsetPx, 0.0f, 0.0f);
              glUniform4f(tuColor, 1.0f, 1.0f, 1.0f, 0.95f);
              glDrawArrays(GL_TRIANGLES, 0, (GLsizei)lcount);

              glUseProgram(0);
              glBindVertexArray(0);
            }
          }
        }

        glDepthMask(GL_TRUE);
        glEnable(GL_DEPTH_TEST);
      }
    }

    if (show_labels) {
      glDisable(GL_DEPTH_TEST);
      glDepthMask(GL_FALSE);
      size_t tv_count = 0;
      for (size_t i = 0; i < sim.count; i++) {
        const Body *b = &sim.bodies[i];
        float sx = 0.0f, sy = 0.0f;
        if (!world_to_screen_px(&world_to_clip, dw, dh, (float)b->x,
                                (float)b->y, (float)b->z, &sx, &sy))
          continue;

        const float scale = 1.0f;
        float radius_px =
            (float)(body_visual_radius_world(b, realistic_scale,
                                             render_radius_scale) *
                    (double)dh / cam.zoom_world_h);
        radius_px = (float)clampd(radius_px, 6.0, 60.0);
        const float text_w = (float)strlen(b->name) * 8.0f * scale;
        const float x = (float)clampd((double)(sx - text_w * 0.5f), 4.0,
                                      (double)dw - 4.0 - text_w);
        const float y = (float)clampd((double)(sy - radius_px - 16.0f * scale),
                                      4.0, (double)dh - 20.0);
        text_append(&text_cpu, &tv_count, &text_cpu_cap, x, y, scale, b->name);
      }

      if (tv_count) {
        glBindVertexArray(text_vao);
        glBindBuffer(GL_ARRAY_BUFFER, text_vbo);
        glBufferData(GL_ARRAY_BUFFER,
                     (GLsizeiptr)(tv_count * sizeof(TextVertex)), text_cpu,
                     GL_STREAM_DRAW);
        glUseProgram(text_prog);
        glUniform2f(tuViewport, (float)dw, (float)dh);
        glActiveTexture(GL_TEXTURE2);
        glBindTexture(GL_TEXTURE_2D, font_tex);

        glUniform2f(tuOffsetPx, 1.0f, 1.0f);
        glUniform4f(tuColor, 0.0f, 0.0f, 0.0f, 0.55f);
        glDrawArrays(GL_TRIANGLES, 0, (GLsizei)tv_count);

        glUniform2f(tuOffsetPx, 0.0f, 0.0f);
        glUniform4f(tuColor, 1.0f, 1.0f, 1.0f, 0.9f);
        glDrawArrays(GL_TRIANGLES, 0, (GLsizei)tv_count);

        glUseProgram(0);
        glBindVertexArray(0);
      }

      glDepthMask(GL_TRUE);
      glEnable(GL_DEPTH_TEST);
    }

    if (show_help) {
      glDisable(GL_DEPTH_TEST);
      glDepthMask(GL_FALSE);
      const float scale = 2.0f;
      const float x0 = 12.0f;
      float y0 = 34.0f;
      const char *lines[] = {
          "F1: toggle help",
          "F2: toggle body labels",
          "F3: realistic scaling",
          "",
          "RMB drag: spawn body (Shift faster)",
          "  Ctrl: circular orbit around selection",
          "  Alt: reverse orbit direction",
          "LMB click: select body",
          "LMB drag: rotate camera",
          "MMB drag / WASD: pan",
          "Wheel: zoom",
          "F: follow selected",
          "T: trails   V: velocity vectors",
          "Space: pause   .: step",
          "- / =: time scale",
          "G / H: gravity",
          "1: two-body   2: galaxy   3: solar system",
          "R: reset solar system   C: clear",
          "Esc: quit",
          NULL,
      };

      size_t tv_count = 0;
      for (int i = 0; lines[i]; i++) {
        text_append(&text_cpu, &tv_count, &text_cpu_cap, x0, y0, scale,
                    lines[i]);
        y0 += 10.0f * scale;
      }

      if (tv_count) {
        glBindVertexArray(text_vao);
        glBindBuffer(GL_ARRAY_BUFFER, text_vbo);
        glBufferData(GL_ARRAY_BUFFER,
                     (GLsizeiptr)(tv_count * sizeof(TextVertex)), text_cpu,
                     GL_STREAM_DRAW);

        glUseProgram(text_prog);
        glUniform2f(tuViewport, (float)dw, (float)dh);
        glActiveTexture(GL_TEXTURE2);
        glBindTexture(GL_TEXTURE_2D, font_tex);

        // Shadow
        glUniform2f(tuOffsetPx, 1.0f, 1.0f);
        glUniform4f(tuColor, 0.0f, 0.0f, 0.0f, 0.55f);
        glDrawArrays(GL_TRIANGLES, 0, (GLsizei)tv_count);

        // Foreground
        glUniform2f(tuOffsetPx, 0.0f, 0.0f);
        glUniform4f(tuColor, 1.0f, 1.0f, 1.0f, 0.88f);
        glDrawArrays(GL_TRIANGLES, 0, (GLsizei)tv_count);

        glUseProgram(0);
        glBindVertexArray(0);
      }

      glDepthMask(GL_TRUE);
      glEnable(GL_DEPTH_TEST);
    }

    // Camera info overlay.
    {
      glDisable(GL_DEPTH_TEST);
      glDepthMask(GL_FALSE);
      char cam_txt[192];
      snprintf(cam_txt, sizeof(cam_txt),
               "cam: (%.3f %.3f %.3f)  yaw=%.3f  pitch=%.3f  zoom=%.3f", cam.cx,
               cam.cy, cam.cz, cam.yaw, cam.pitch, cam.zoom_world_h);

      const float scale = 2.5f;
      const float text_w = (float)strlen(cam_txt) * 8.0f * scale;
      const float x = (float)clampd(((double)dw - (double)text_w) * 0.5, 4.0,
                                    (double)dw - 4.0 - text_w);
      const float y = (float)dh - (8.0f * scale) - 8.0f - 22.0f;

      size_t tv_count = 0;
      text_append(&text_cpu, &tv_count, &text_cpu_cap, x, y, scale, cam_txt);
      if (tv_count) {
        glBindVertexArray(text_vao);
        glBindBuffer(GL_ARRAY_BUFFER, text_vbo);
        glBufferData(GL_ARRAY_BUFFER,
                     (GLsizeiptr)(tv_count * sizeof(TextVertex)), text_cpu,
                     GL_STREAM_DRAW);

        glUseProgram(text_prog);
        glUniform2f(tuViewport, (float)dw, (float)dh);
        glActiveTexture(GL_TEXTURE2);
        glBindTexture(GL_TEXTURE_2D, font_tex);

        glUniform2f(tuOffsetPx, 2.0f, 2.0f);
        glUniform4f(tuColor, 0.0f, 0.0f, 0.0f, 0.55f);
        glDrawArrays(GL_TRIANGLES, 0, (GLsizei)tv_count);

        glUniform2f(tuOffsetPx, 0.0f, 0.0f);
        glUniform4f(tuColor, 1.0f, 1.0f, 1.0f, 0.85f);
        glDrawArrays(GL_TRIANGLES, 0, (GLsizei)tv_count);

        glUseProgram(0);
        glBindVertexArray(0);
      }

      glDepthMask(GL_TRUE);
      glEnable(GL_DEPTH_TEST);
    }

    // Render backend debug info.
    {
      glDisable(GL_DEPTH_TEST);
      glDepthMask(GL_FALSE);
      char dbg[128];
      const int scene_bits =
          use_offscreen_scene_fbo ? scene_depth_bits : (int)default_depth_bits;
      snprintf(dbg, sizeof(dbg), "Depth: default=%d scene=%d%s",
               (int)default_depth_bits, scene_bits,
               use_offscreen_scene_fbo ? " (FBO)" : "");

      const float scale = 2.0f;
      const float x = 12.0f;
      const float y = (float)dh - (8.0f * scale) - 6.0f;
      size_t tv_count = 0;
      text_append(&text_cpu, &tv_count, &text_cpu_cap, x, y, scale, dbg);
      if (tv_count) {
        glBindVertexArray(text_vao);
        glBindBuffer(GL_ARRAY_BUFFER, text_vbo);
        glBufferData(GL_ARRAY_BUFFER,
                     (GLsizeiptr)(tv_count * sizeof(TextVertex)), text_cpu,
                     GL_STREAM_DRAW);

        glUseProgram(text_prog);
        glUniform2f(tuViewport, (float)dw, (float)dh);
        glActiveTexture(GL_TEXTURE2);
        glBindTexture(GL_TEXTURE_2D, font_tex);

        glUniform2f(tuOffsetPx, 1.0f, 1.0f);
        glUniform4f(tuColor, 0.0f, 0.0f, 0.0f, 0.55f);
        glDrawArrays(GL_TRIANGLES, 0, (GLsizei)tv_count);

        glUniform2f(tuOffsetPx, 0.0f, 0.0f);
        glUniform4f(tuColor, 1.0f, 1.0f, 1.0f, 0.75f);
        glDrawArrays(GL_TRIANGLES, 0, (GLsizei)tv_count);

        glUseProgram(0);
        glBindVertexArray(0);
      }

      glDepthMask(GL_TRUE);
      glEnable(GL_DEPTH_TEST);
    }

    // Scaling mode overlay.
    {
      glDisable(GL_DEPTH_TEST);
      glDepthMask(GL_FALSE);
      const char *mode = realistic_scale ? "Scale: realistic (F3)" : "Scale: exaggerated (F3)";
      const float scale = 2.0f;
      const float text_w = (float)strlen(mode) * 8.0f * scale;
      const float x = (float)clampd(((double)dw - (double)text_w) * 0.5, 4.0, (double)dw - 4.0 - text_w);
      const float y = (float)dh - (8.0f * scale) - 6.0f;

      size_t tv_count = 0;
      text_append(&text_cpu, &tv_count, &text_cpu_cap, x, y, scale, mode);
      if (tv_count) {
        glBindVertexArray(text_vao);
        glBindBuffer(GL_ARRAY_BUFFER, text_vbo);
        glBufferData(GL_ARRAY_BUFFER, (GLsizeiptr)(tv_count * sizeof(TextVertex)), text_cpu, GL_STREAM_DRAW);

        glUseProgram(text_prog);
        glUniform2f(tuViewport, (float)dw, (float)dh);
        glActiveTexture(GL_TEXTURE2);
        glBindTexture(GL_TEXTURE_2D, font_tex);

        glUniform2f(tuOffsetPx, 1.5f, 1.5f);
        glUniform4f(tuColor, 0.0f, 0.0f, 0.0f, 0.55f);
        glDrawArrays(GL_TRIANGLES, 0, (GLsizei)tv_count);

        glUniform2f(tuOffsetPx, 0.0f, 0.0f);
        glUniform4f(tuColor, 1.0f, 1.0f, 1.0f, 0.80f);
        glDrawArrays(GL_TRIANGLES, 0, (GLsizei)tv_count);

        glUseProgram(0);
        glBindVertexArray(0);
      }
      glDepthMask(GL_TRUE);
      glEnable(GL_DEPTH_TEST);
    }

    SDL_GL_SwapWindow(win);
  }

  sim_destroy(&sim);
  free(trails);
  free(particle_cpu);
  free(ring_cpu);
  free(ring_tex_cpu);
  free(ring_radial_cpu);
  if (ring_cache) {
    for (size_t i = 0; i < ring_cache_count; i++) {
      if (ring_cache[i].tex)
        glDeleteTextures(1, &ring_cache[i].tex);
    }
  }
  free(ring_cache);
  free(line_cpu);
  free(text_cpu);
  glDeleteBuffers(1, &particle_vbo);
  glDeleteVertexArrays(1, &particle_vao);
  glDeleteBuffers(1, &ring_vbo);
  glDeleteVertexArrays(1, &ring_vao);
  glDeleteBuffers(1, &line_vbo);
  glDeleteVertexArrays(1, &line_vao);
  glDeleteVertexArrays(1, &bg_vao);
  glDeleteBuffers(1, &text_vbo);
  glDeleteVertexArrays(1, &text_vao);
  glDeleteTextures(1, &tex_array);
  glDeleteTextures(1, &env_tex);
  glDeleteTextures(1, &font_tex);
  scene_fbo_destroy(&scene_fbo, &scene_color, &scene_depth, &scene_w, &scene_h);
  glDeleteProgram(particle_prog);
  glDeleteProgram(line_prog);
  glDeleteProgram(bg_prog);
  glDeleteProgram(text_prog);
  glDeleteProgram(ring_prog);
  free(assets_dir);
  SDL_GL_DeleteContext(ctx);
  SDL_DestroyWindow(win);
  SDL_Quit();
  return 0;
}
