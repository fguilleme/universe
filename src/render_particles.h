#pragma once

#include <GL/glew.h>

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

#include "mat4.h"

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
  GLuint prog;
  GLint uWorldToClip;
  GLint uCamRight;
  GLint uCamUp;
  GLint uTex;
  GLint uTime;
  GLuint vao;
  GLuint vbo;
} ParticleRenderer;

bool particles_init(ParticleRenderer *p);
void particles_destroy(ParticleRenderer *p);

void particles_draw(ParticleRenderer *p, const Mat4 *world_to_clip, float time,
                    float cam_right_x, float cam_right_y, float cam_right_z,
                    float cam_up_x, float cam_up_y, float cam_up_z,
                    GLuint tex_array, const ParticleGPU *particles,
                    size_t particle_count);

