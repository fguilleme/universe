#pragma once

#include <GL/glew.h>

#include <stdbool.h>
#include <stddef.h>

#include "camera.h"
#include "mat4.h"
#include "sim.h"

typedef struct {
  float x;
  float y;
  float z;
  float r;
  float g;
  float b;
  float a;
} LineVertex;

typedef struct {
  GLuint prog;
  GLint uWorldToClip;
  GLuint vao;
  GLuint vbo;
} LineRenderer;

bool lines_init(LineRenderer *r);
void lines_destroy(LineRenderer *r);

void lines_draw(LineRenderer *r, const Mat4 *world_to_clip, const LineVertex *verts,
                size_t count);

