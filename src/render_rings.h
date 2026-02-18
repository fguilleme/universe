#pragma once

#include <GL/glew.h>

#include <stdbool.h>

#include "mat4.h"
#include "sim.h"

typedef struct {
  GLuint prog;
  GLint uWorldToClip;
  GLint uRingTex;
  GLint uRadialStrip;
  GLuint vao;
  GLuint vbo;

  // Per-frame instances.
  struct RingGPU *cpu;
  GLuint *tex;
  uint8_t *radial;
  size_t count;
  size_t cap;

  // Texture cache.
  struct RingTexCache *cache;
  size_t cache_count;
  size_t cache_cap;
} RingRenderer;

bool rings_init(RingRenderer *r);
void rings_destroy(RingRenderer *r);

// Builds ring instances for any body with a matching ring texture and draws.
void rings_build_and_draw(RingRenderer *r, const Mat4 *world_to_clip,
                          const Sim *sim, const char *assets_dir,
                          bool realistic_scale, double render_radius_scale);

