#include "textures.h"

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define STB_IMAGE_IMPLEMENTATION
#include "../third_party/stb_image.h"

static double clampd(double x, double lo, double hi) {
  if (x < lo) return lo;
  if (x > hi) return hi;
  return x;
}

static uint32_t xorshift32(uint32_t *state) {
  uint32_t x = *state;
  x ^= x << 13;
  x ^= x >> 17;
  x ^= x << 5;
  *state = x;
  return x;
}

static void fill_tex_planet(uint8_t *dst, int w, int h, uint32_t seed, float base_r, float base_g, float base_b,
                            float stripe_strength) {
  uint32_t st = seed ? seed : 1u;
  for (int y = 0; y < h; y++) {
    for (int x = 0; x < w; x++) {
      st = xorshift32(&st);
      const float n = (float)(st & 0xFFFF) / 65535.0f;
      const float u = (float)x / (float)(w - 1);
      const float v = (float)y / (float)(h - 1);
      const float stripes = stripe_strength * (0.5f + 0.5f * sinf(v * 18.0f + n * 1.7f));
      const float shade = 0.82f + 0.18f * (n * 0.9f + 0.1f * sinf(u * 30.0f));
      float r = base_r * shade + stripes * 0.06f;
      float g = base_g * shade + stripes * 0.05f;
      float b = base_b * shade + stripes * 0.03f;
      r = (float)clampd(r, 0.0, 1.0);
      g = (float)clampd(g, 0.0, 1.0);
      b = (float)clampd(b, 0.0, 1.0);
      const int i = (y * w + x) * 4;
      dst[i + 0] = (uint8_t)(r * 255.0f);
      dst[i + 1] = (uint8_t)(g * 255.0f);
      dst[i + 2] = (uint8_t)(b * 255.0f);
      dst[i + 3] = 255;
    }
  }
}

static void resample_bilinear_rgba(uint8_t *dst, int dw, int dh, const uint8_t *src, int sw, int sh) {
  for (int y = 0; y < dh; y++) {
    const float v = (dh == 1) ? 0.0f : (float)y / (float)(dh - 1);
    const float sy = v * (float)(sh - 1);
    const int y0 = (int)floorf(sy);
    const int y1 = (y0 + 1 < sh) ? (y0 + 1) : y0;
    const float fy = sy - (float)y0;
    for (int x = 0; x < dw; x++) {
      const float u = (dw == 1) ? 0.0f : (float)x / (float)(dw - 1);
      const float sx = u * (float)(sw - 1);
      const int x0 = (int)floorf(sx);
      const int x1 = (x0 + 1 < sw) ? (x0 + 1) : x0;
      const float fx = sx - (float)x0;

      const uint8_t *p00 = &src[(y0 * sw + x0) * 4];
      const uint8_t *p10 = &src[(y0 * sw + x1) * 4];
      const uint8_t *p01 = &src[(y1 * sw + x0) * 4];
      const uint8_t *p11 = &src[(y1 * sw + x1) * 4];

      for (int c = 0; c < 4; c++) {
        const float a = (float)p00[c] * (1.0f - fx) + (float)p10[c] * fx;
        const float b = (float)p01[c] * (1.0f - fx) + (float)p11[c] * fx;
        const float vout = a * (1.0f - fy) + b * fy;
        dst[(y * dw + x) * 4 + c] = (uint8_t)clampd(vout + 0.5f, 0.0, 255.0);
      }
    }
  }
}

static uint8_t *try_load_rgba(const char *dir, const char *filename, int *w, int *h) {
  char path[1024];
  if (dir && dir[0]) {
    snprintf(path, sizeof(path), "%s/%s", dir, filename);
  } else {
    snprintf(path, sizeof(path), "%s", filename);
  }
  int comp = 0;
  uint8_t *data = stbi_load(path, w, h, &comp, 4);
  return data;
}

static GLuint upload_rgba2d(const uint8_t *rgba, int w, int h) {
  GLuint tex = 0;
  glGenTextures(1, &tex);
  glBindTexture(GL_TEXTURE_2D, tex);
  glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, w, h, 0, GL_RGBA, GL_UNSIGNED_BYTE, rgba);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
  return tex;
}

GLuint textures_try_load_rgba2d(const char *dir, const char *label, const char *const *candidates, int *out_w, int *out_h) {
  if (!candidates) return 0;
  for (size_t i = 0; candidates[i]; i++) {
    int w = 0, h = 0;
    uint8_t *data = try_load_rgba(dir, candidates[i], &w, &h);
    if (!data) continue;
    GLuint tex = upload_rgba2d(data, w, h);
    stbi_image_free(data);
    if (tex) {
      fprintf(stderr, "Loaded %s texture from %s/%s (%dx%d)\n", label ? label : "image", dir ? dir : ".", candidates[i], w,
              h);
    }
    if (out_w) *out_w = w;
    if (out_h) *out_h = h;
    return tex;
  }
  return 0;
}

static void upload_layer(GLuint tex, int layer, const uint8_t *rgba, int w, int h) {
  (void)tex;
  glTexSubImage3D(GL_TEXTURE_2D_ARRAY, 0, 0, 0, layer, w, h, 1, GL_RGBA, GL_UNSIGNED_BYTE, rgba);
}

GLuint textures_create_planet_array(const char *dir) {
  const int tw = 512;
  const int th = 512;
  uint8_t *tmp = (uint8_t *)malloc((size_t)tw * (size_t)th * 4);
  if (!tmp) return 0;

  GLuint tex = 0;
  glGenTextures(1, &tex);
  glBindTexture(GL_TEXTURE_2D_ARRAY, tex);
  glTexImage3D(GL_TEXTURE_2D_ARRAY, 0, GL_RGBA8, tw, th, TEX_LAYER_COUNT, 0, GL_RGBA, GL_UNSIGNED_BYTE, NULL);

  // Layer 0: generic
  fill_tex_planet(tmp, tw, th, 0xA5A5A5A5u, 0.85f, 0.85f, 0.90f, 0.0f);
  upload_layer(tex, TEX_GENERIC, tmp, tw, th);

  typedef struct {
    int layer;
    const char *fallback_name;
    uint32_t seed;
    float r, g, b;
    float stripes;
    const char *candidates[8];
  } LayerSpec;

  const LayerSpec specs[] = {
      {TEX_SUN,
       "sun",
       0xBADC0FFEu,
       1.00f,
       0.92f,
       0.55f,
       0.15f,
       {"2k_sun.jpg", "2k_sun.png", "sun.jpg", "sun.png", NULL}},
      {TEX_MERCURY,
       "mercury",
       0x11111111u,
       0.55f,
       0.53f,
       0.50f,
       0.0f,
       {"2k_mercury.jpg", "2k_mercury.png", "mercury.jpg", "mercury.png", NULL}},
      {TEX_VENUS,
       "venus",
       0x22222222u,
       0.88f,
       0.78f,
       0.55f,
       0.0f,
       {"2k_venus_surface.jpg", "2k_venus_surface.png", "2k_venus_atmosphere.jpg", "venus.jpg", "venus.png", NULL}},
      {TEX_EARTH,
       "earth",
       0x33333333u,
       0.45f,
       0.65f,
       0.95f,
       0.05f,
       {"2k_earth_daymap.jpg", "2k_earth_daymap.png", "earth.jpg", "earth.png", NULL}},
      {TEX_MOON,
       "moon",
       0x3A3A3A3Au,
       0.70f,
       0.70f,
       0.72f,
       0.0f,
       {"2k_moon.jpg", "2k_moon.png", "moon.jpg", "moon.png", NULL}},
      {TEX_MARS,
       "mars",
       0x44444444u,
       0.82f,
       0.46f,
       0.28f,
       0.0f,
       {"2k_mars.jpg", "2k_mars.png", "mars.jpg", "mars.png", NULL}},
      {TEX_JUPITER,
       "jupiter",
       0x55555555u,
       0.90f,
       0.78f,
       0.55f,
       0.35f,
       {"2k_jupiter.jpg", "2k_jupiter.png", "jupiter.jpg", "jupiter.png", NULL}},
      {TEX_SATURN,
       "saturn",
       0x66666666u,
       0.92f,
       0.84f,
       0.62f,
       0.25f,
       {"2k_saturn.jpg", "2k_saturn.png", "saturn.jpg", "saturn.png", NULL}},
      {TEX_URANUS,
       "uranus",
       0x77777777u,
       0.70f,
       0.90f,
       0.92f,
       0.05f,
       {"2k_uranus.jpg", "2k_uranus.png", "uranus.jpg", "uranus.png", NULL}},
      {TEX_NEPTUNE,
       "neptune",
       0x88888888u,
       0.35f,
       0.55f,
       0.95f,
       0.05f,
       {"2k_neptune.jpg", "2k_neptune.png", "neptune.jpg", "neptune.png", NULL}},
  };

  for (size_t i = 0; i < sizeof(specs) / sizeof(specs[0]); i++) {
    const LayerSpec *s = &specs[i];
    int sw = 0, sh = 0;
    uint8_t *src = NULL;
    for (size_t c = 0; s->candidates[c]; c++) {
      src = try_load_rgba(dir, s->candidates[c], &sw, &sh);
      if (src) {
        fprintf(stderr, "Loaded texture layer %d from %s/%s (%dx%d)\n", s->layer, dir ? dir : ".", s->candidates[c], sw,
                sh);
        break;
      }
    }

    if (src) {
      resample_bilinear_rgba(tmp, tw, th, src, sw, sh);
      stbi_image_free(src);
      upload_layer(tex, s->layer, tmp, tw, th);
    } else {
      fprintf(stderr, "Missing texture for %s; using procedural fallback\n", s->fallback_name);
      fill_tex_planet(tmp, tw, th, s->seed, s->r, s->g, s->b, s->stripes);
      upload_layer(tex, s->layer, tmp, tw, th);
    }
  }

  glTexParameteri(GL_TEXTURE_2D_ARRAY, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D_ARRAY, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D_ARRAY, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
  glTexParameteri(GL_TEXTURE_2D_ARRAY, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

  free(tmp);
  return tex;
}

GLuint textures_load_hdr_equirect(const char *path) {
  int w = 0, h = 0, comp = 0;
  float *data = stbi_loadf(path, &w, &h, &comp, 3);
  if (!data) {
    fprintf(stderr, "Failed to load HDR %s: %s\n", path, stbi_failure_reason());
    return 0;
  }

  GLuint tex = 0;
  glGenTextures(1, &tex);
  glBindTexture(GL_TEXTURE_2D, tex);
  glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB16F, w, h, 0, GL_RGB, GL_FLOAT, data);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

  stbi_image_free(data);
  return tex;
}
