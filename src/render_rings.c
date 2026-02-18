#include "render_rings.h"

#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "glutil.h"
#include "shaders.h"
#include "textures.h"
#include "util.h"

typedef struct RingGPU {
  float x;
  float y;
  float z;
  float outer;
  float inner_ratio;
  float nx;
  float ny;
  float nz;
} RingGPU;

typedef struct RingTexCache {
  char name_lc[32];
  GLuint tex;
  bool radial_strip;
  bool attempted;
} RingTexCache;

static GLuint ring_texture_for_body(const char *assets_dir, RingTexCache **cache,
                                    size_t *count, size_t *cap, const Body *b,
                                    bool *out_radial_strip) {
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

bool rings_init(RingRenderer *r) {
  if (!r)
    return false;
  memset(r, 0, sizeof(*r));

  GLuint vs = glutil_compile_shader(GL_VERTEX_SHADER, kRingsVS);
  GLuint fs = glutil_compile_shader(GL_FRAGMENT_SHADER, kRingsFS);
  if (!vs || !fs)
    return false;
  r->prog = glutil_link_program(vs, fs);
  glDeleteShader(vs);
  glDeleteShader(fs);
  if (!r->prog)
    return false;

  r->uWorldToClip = glGetUniformLocation(r->prog, "uWorldToClip");
  r->uRingTex = glGetUniformLocation(r->prog, "uRingTex");
  r->uRadialStrip = glGetUniformLocation(r->prog, "uRadialStrip");

  glGenVertexArrays(1, &r->vao);
  glGenBuffers(1, &r->vbo);
  glBindVertexArray(r->vao);
  glBindBuffer(GL_ARRAY_BUFFER, r->vbo);
  glBufferData(GL_ARRAY_BUFFER, 1, NULL, GL_STREAM_DRAW);
  glEnableVertexAttribArray(0);
  glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(RingGPU),
                        (void *)offsetof(RingGPU, x));
  glEnableVertexAttribArray(1);
  glVertexAttribPointer(1, 1, GL_FLOAT, GL_FALSE, sizeof(RingGPU),
                        (void *)offsetof(RingGPU, outer));
  glEnableVertexAttribArray(2);
  glVertexAttribPointer(2, 1, GL_FLOAT, GL_FALSE, sizeof(RingGPU),
                        (void *)offsetof(RingGPU, inner_ratio));
  glEnableVertexAttribArray(3);
  glVertexAttribPointer(3, 3, GL_FLOAT, GL_FALSE, sizeof(RingGPU),
                        (void *)offsetof(RingGPU, nx));
  glVertexAttribDivisor(0, 1);
  glVertexAttribDivisor(1, 1);
  glVertexAttribDivisor(2, 1);
  glVertexAttribDivisor(3, 1);
  glBindVertexArray(0);

  glUseProgram(r->prog);
  glUniform1i(r->uRingTex, 3);
  glUniform1i(r->uRadialStrip, 0);
  glUseProgram(0);
  return true;
}

void rings_destroy(RingRenderer *r) {
  if (!r)
    return;

  if (r->cache) {
    for (size_t i = 0; i < r->cache_count; i++) {
      if (r->cache[i].tex)
        glDeleteTextures(1, &r->cache[i].tex);
    }
  }
  free(r->cache);
  free(r->cpu);
  free(r->tex);
  free(r->radial);

  if (r->vbo)
    glDeleteBuffers(1, &r->vbo);
  if (r->vao)
    glDeleteVertexArrays(1, &r->vao);
  if (r->prog)
    glDeleteProgram(r->prog);

  memset(r, 0, sizeof(*r));
}

static void rings_instances_clear(RingRenderer *r) { r->count = 0; }

static void rings_instances_push(RingRenderer *r, const RingGPU *gpu, GLuint tex,
                                 bool radial_strip) {
  if (r->count + 1 > r->cap) {
    size_t next = r->cap ? (r->cap * 2) : 16;
    while (next < r->count + 1)
      next *= 2;
    RingGPU *ng = (RingGPU *)realloc(r->cpu, next * sizeof(RingGPU));
    GLuint *nt = (GLuint *)realloc(r->tex, next * sizeof(GLuint));
    uint8_t *nm = (uint8_t *)realloc(r->radial, next * sizeof(uint8_t));
    if (!ng || !nt || !nm) {
      free(ng);
      free(nt);
      free(nm);
      return;
    }
    r->cpu = ng;
    r->tex = nt;
    r->radial = nm;
    r->cap = next;
  }
  r->cpu[r->count] = *gpu;
  r->tex[r->count] = tex;
  r->radial[r->count] = radial_strip ? 1u : 0u;
  r->count++;
}

void rings_build_and_draw(RingRenderer *r, const Mat4 *world_to_clip,
                          const Sim *sim, const char *assets_dir,
                          bool realistic_scale, double render_radius_scale) {
  if (!r || !r->prog || !world_to_clip || !sim || !assets_dir)
    return;

  rings_instances_clear(r);

  for (size_t i = 0; i < sim->count; i++) {
    const Body *b = &sim->bodies[i];
    bool radial_strip = false;
    GLuint ring_tex = ring_texture_for_body(assets_dir, &r->cache,
                                            &r->cache_count, &r->cache_cap, b,
                                            &radial_strip);
    if (!ring_tex)
      continue;

    const double pr =
        body_visual_radius_world(b, realistic_scale, render_radius_scale);
    const double outer = pr * 2.5;
    const double inner = pr * 1.2;

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

    RingGPU gpu = {.x = (float)b->x,
                   .y = (float)b->y,
                   .z = (float)b->z,
                   .outer = (float)outer,
                   .inner_ratio = (float)(inner / outer),
                   .nx = nx,
                   .ny = ny,
                   .nz = nz};
    rings_instances_push(r, &gpu, ring_tex, radial_strip);
  }

  if (!r->count)
    return;

  glUseProgram(r->prog);
  glUniformMatrix4fv(r->uWorldToClip, 1, GL_FALSE, world_to_clip->m);
  glActiveTexture(GL_TEXTURE3);
  glBindVertexArray(r->vao);
  glBindBuffer(GL_ARRAY_BUFFER, r->vbo);

  glEnable(GL_BLEND);
  glDepthMask(GL_FALSE);
  glDepthFunc(GL_LEQUAL);

  for (size_t i = 0; i < r->count; i++) {
    glBindTexture(GL_TEXTURE_2D, r->tex[i]);
    glUniform1i(r->uRadialStrip, r->radial[i] ? 1 : 0);
    glBufferData(GL_ARRAY_BUFFER, (GLsizeiptr)sizeof(RingGPU), &r->cpu[i],
                 GL_STREAM_DRAW);
    glDrawArraysInstanced(GL_TRIANGLE_STRIP, 0, 4, 1);
  }

  glDepthFunc(GL_LESS);
  glDepthMask(GL_TRUE);
  glBindVertexArray(0);
  glUseProgram(0);
}
