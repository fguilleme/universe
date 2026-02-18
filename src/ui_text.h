#pragma once

#include <GL/glew.h>

#include <stddef.h>

typedef struct {
  float x;
  float y;
  float u;
  float v;
} TextVertex;

GLuint create_font_texture(void);

void text_append(TextVertex **verts, size_t *count, size_t *cap, float x,
                 float y, float scale, const char *text);

