#pragma once

#include <GL/glew.h>

#include <stdbool.h>

#include "camera.h"

typedef struct {
  GLuint prog;
  GLint uEnv;
  GLint uHasEnv;
  GLint uInvViewRot;
  GLint uAspect;
  GLint uExposure;
  GLuint vao;
  float exposure;
} BgRenderer;

// Compiles and links the background shader program.
bool bg_init(BgRenderer *bg);

void bg_destroy(BgRenderer *bg);

// Draws the background, using `env_tex` if non-zero.
void bg_draw(const BgRenderer *bg, const Camera *cam, int dw, int dh,
             GLuint env_tex);

