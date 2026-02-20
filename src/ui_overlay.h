#pragma once

#include <GL/glew.h>

#include <stdbool.h>

#include "camera.h"
#include "mat4.h"
#include "ui_text.h"

typedef struct {
  GLuint prog;
  GLint uViewport;
  GLint uOffsetPx;
  GLint uFont;
  GLint uColor;
  GLuint vao;
  GLuint vbo;
  GLuint font_tex;
} UiTextRenderer;

bool ui_text_renderer_init(UiTextRenderer *r);
void ui_text_renderer_destroy(UiTextRenderer *r);

// Draws a batch of text vertices in screen space.
void ui_text_renderer_draw(const UiTextRenderer *tr, const TextVertex *verts,
                           size_t count, int dw, int dh, float offx, float offy,
                           float r, float g, float b, float a);

typedef struct {
  // Solar system initial conditions are for 2000-01-01 00:00.
  double solar_epoch_jd;
} UiOverlayConfig;

typedef struct {
  int preset;
  double sim_time;
  bool show_help;
  bool show_labels; // unused here but kept for future
  bool follow_selected;
  bool realistic_scale;
  Camera cam;
  int default_depth_bits;
  bool use_offscreen_scene_fbo;
  int scene_depth_bits;
} UiOverlayState;

void ui_overlay_draw(const UiTextRenderer *tr, const UiOverlayConfig *cfg,
                     const UiOverlayState *st, int dw, int dh);
