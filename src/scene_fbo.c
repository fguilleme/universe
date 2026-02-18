#include "scene_fbo.h"

#include <stdio.h>

void scene_fbo_destroy(SceneFbo *s) {
  if (!s)
    return;
  if (s->depth_rb)
    glDeleteRenderbuffers(1, &s->depth_rb);
  if (s->color_tex)
    glDeleteTextures(1, &s->color_tex);
  if (s->fbo)
    glDeleteFramebuffers(1, &s->fbo);
  s->fbo = 0;
  s->color_tex = 0;
  s->depth_rb = 0;
  s->w = 0;
  s->h = 0;
  s->depth_bits = 0;
}

bool scene_fbo_ensure(SceneFbo *s, int w, int h) {
  if (!s)
    return false;
  if (w <= 0 || h <= 0)
    return false;
  if (s->fbo && s->w == w && s->h == h)
    return true;

  scene_fbo_destroy(s);

  glGenFramebuffers(1, &s->fbo);
  glBindFramebuffer(GL_FRAMEBUFFER, s->fbo);

  glGenTextures(1, &s->color_tex);
  glBindTexture(GL_TEXTURE_2D, s->color_tex);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
  glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, w, h, 0, GL_RGBA, GL_UNSIGNED_BYTE,
               NULL);
  glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D,
                         s->color_tex, 0);

  glGenRenderbuffers(1, &s->depth_rb);
  glBindRenderbuffer(GL_RENDERBUFFER, s->depth_rb);
  glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH_COMPONENT24, w, h);
  GLint depth_size = 0;
  glGetRenderbufferParameteriv(GL_RENDERBUFFER, GL_RENDERBUFFER_DEPTH_SIZE,
                               &depth_size);
  s->depth_bits = (int)depth_size;
  glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT,
                            GL_RENDERBUFFER, s->depth_rb);

  const GLenum status = glCheckFramebufferStatus(GL_FRAMEBUFFER);
  glBindFramebuffer(GL_FRAMEBUFFER, 0);
  if (status != GL_FRAMEBUFFER_COMPLETE) {
    fprintf(stderr, "Scene FBO incomplete: 0x%x\n", (unsigned)status);
    scene_fbo_destroy(s);
    return false;
  }

  s->w = w;
  s->h = h;
  return true;
}

