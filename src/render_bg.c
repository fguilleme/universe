#include "render_bg.h"

#include <math.h>
#include <string.h>

#include "glutil.h"
#include "shaders.h"

bool bg_init(BgRenderer *bg) {
  if (!bg)
    return false;
  memset(bg, 0, sizeof(*bg));
  bg->exposure = 1.2f;

  GLuint bvs = glutil_compile_shader(GL_VERTEX_SHADER, kBgVS);
  GLuint bfs = glutil_compile_shader(GL_FRAGMENT_SHADER, kBgFS);
  if (!bvs || !bfs)
    return false;
  bg->prog = glutil_link_program(bvs, bfs);
  glDeleteShader(bvs);
  glDeleteShader(bfs);
  if (!bg->prog)
    return false;

  bg->uEnv = glGetUniformLocation(bg->prog, "uEnv");
  bg->uHasEnv = glGetUniformLocation(bg->prog, "uHasEnv");
  bg->uInvViewRot = glGetUniformLocation(bg->prog, "uInvViewRot");
  bg->uAspect = glGetUniformLocation(bg->prog, "uAspect");
  bg->uExposure = glGetUniformLocation(bg->prog, "uExposure");

  glGenVertexArrays(1, &bg->vao);
  glBindVertexArray(0);

  glUseProgram(bg->prog);
  glUniform1i(bg->uEnv, 1);
  glUniform1f(bg->uExposure, bg->exposure);
  glUseProgram(0);
  return true;
}

void bg_destroy(BgRenderer *bg) {
  if (!bg)
    return;
  if (bg->vao)
    glDeleteVertexArrays(1, &bg->vao);
  if (bg->prog)
    glDeleteProgram(bg->prog);
  memset(bg, 0, sizeof(*bg));
}

void bg_draw(const BgRenderer *bg, const Camera *cam, int dw, int dh,
             GLuint env_tex) {
  if (!bg || !bg->prog || !cam)
    return;

  const float cy = cosf((float)-cam->yaw);
  const float sy = sinf((float)-cam->yaw);
  const float cp = cosf((float)-cam->pitch);
  const float sp = sinf((float)-cam->pitch);
  // uInvViewRot = Rz(-yaw) * Rx(-pitch)
  float invR[9];
  invR[0] = cy;
  invR[1] = sy;
  invR[2] = 0.0f;
  invR[3] = -sy * cp;
  invR[4] = cy * cp;
  invR[5] = sp;
  invR[6] = sy * sp;
  invR[7] = -cy * sp;
  invR[8] = cp;

  glBindVertexArray(bg->vao);
  glDisable(GL_DEPTH_TEST);
  glDepthMask(GL_FALSE);
  glDisable(GL_BLEND);
  glUseProgram(bg->prog);
  glUniformMatrix3fv(bg->uInvViewRot, 1, GL_FALSE, invR);
  glUniform1f(bg->uAspect, (dh > 0) ? ((float)dw / (float)dh) : 1.0f);
  glUniform1i(bg->uHasEnv, env_tex ? 1 : 0);
  glActiveTexture(GL_TEXTURE1);
  glBindTexture(GL_TEXTURE_2D, env_tex);
  glDrawArrays(GL_TRIANGLES, 0, 3);
  glUseProgram(0);
  glBindVertexArray(0);
  glDepthMask(GL_TRUE);
  glEnable(GL_DEPTH_TEST);
  glEnable(GL_BLEND);
}

