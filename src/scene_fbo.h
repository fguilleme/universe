#pragma once

#include <GL/glew.h>

#include <stdbool.h>

typedef struct {
  GLuint fbo;
  GLuint color_tex;
  GLuint depth_rb;
  int w;
  int h;
  int depth_bits;
} SceneFbo;

void scene_fbo_destroy(SceneFbo *s);

// Ensures an offscreen (color+depth) framebuffer of size w*h.
// Returns false on failure and leaves `s` cleared.
bool scene_fbo_ensure(SceneFbo *s, int w, int h);

