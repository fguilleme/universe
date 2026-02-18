#include "ui_text.h"

#include <stdint.h>
#include <stdlib.h>

#include "font8x8_basic.h"

GLuint create_font_texture(void) {
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

void text_append(TextVertex **verts, size_t *count, size_t *cap, float x,
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

    (*verts)[(*count)++] = (TextVertex){.x = x0, .y = y0, .u = u0, .v = v0};
    (*verts)[(*count)++] = (TextVertex){.x = x1, .y = y0, .u = u1, .v = v0};
    (*verts)[(*count)++] = (TextVertex){.x = x1, .y = y1, .u = u1, .v = v1};
    (*verts)[(*count)++] = (TextVertex){.x = x0, .y = y0, .u = u0, .v = v0};
    (*verts)[(*count)++] = (TextVertex){.x = x1, .y = y1, .u = u1, .v = v1};
    (*verts)[(*count)++] = (TextVertex){.x = x0, .y = y1, .u = u0, .v = v1};

    x += gw;
  }
}

