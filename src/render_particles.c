#include "render_particles.h"

#include <stddef.h>
#include <string.h>

#include "glutil.h"
#include "shaders.h"

bool particles_init(ParticleRenderer *p) {
  if (!p)
    return false;
  memset(p, 0, sizeof(*p));

  GLuint vs = glutil_compile_shader(GL_VERTEX_SHADER, kParticlesVS);
  GLuint fs = glutil_compile_shader(GL_FRAGMENT_SHADER, kParticlesFS);
  if (!vs || !fs)
    return false;
  p->prog = glutil_link_program(vs, fs);
  glDeleteShader(vs);
  glDeleteShader(fs);
  if (!p->prog)
    return false;

  p->uWorldToClip = glGetUniformLocation(p->prog, "uWorldToClip");
  p->uCamRight = glGetUniformLocation(p->prog, "uCamRight");
  p->uCamUp = glGetUniformLocation(p->prog, "uCamUp");
  p->uTex = glGetUniformLocation(p->prog, "uTex");
  p->uTime = glGetUniformLocation(p->prog, "uTime");

  glGenVertexArrays(1, &p->vao);
  glGenBuffers(1, &p->vbo);
  glBindVertexArray(p->vao);
  glBindBuffer(GL_ARRAY_BUFFER, p->vbo);
  glBufferData(GL_ARRAY_BUFFER, 1, NULL, GL_STREAM_DRAW);
  glEnableVertexAttribArray(0);
  glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(ParticleGPU),
                        (void *)offsetof(ParticleGPU, x));
  glEnableVertexAttribArray(1);
  glVertexAttribPointer(1, 1, GL_FLOAT, GL_FALSE, sizeof(ParticleGPU),
                        (void *)offsetof(ParticleGPU, radius));
  glEnableVertexAttribArray(2);
  glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, sizeof(ParticleGPU),
                        (void *)offsetof(ParticleGPU, r));
  glEnableVertexAttribArray(3);
  glVertexAttribIPointer(3, 1, GL_UNSIGNED_INT, sizeof(ParticleGPU),
                         (void *)offsetof(ParticleGPU, layer));
  glEnableVertexAttribArray(4);
  glVertexAttribPointer(4, 3, GL_FLOAT, GL_FALSE, sizeof(ParticleGPU),
                        (void *)offsetof(ParticleGPU, tilt));
  glVertexAttribDivisor(0, 1);
  glVertexAttribDivisor(1, 1);
  glVertexAttribDivisor(2, 1);
  glVertexAttribDivisor(3, 1);
  glVertexAttribDivisor(4, 1);
  glBindVertexArray(0);

  glUseProgram(p->prog);
  glUniform1i(p->uTex, 0);
  glUniform1f(p->uTime, 0.0f);
  glUseProgram(0);
  return true;
}

void particles_destroy(ParticleRenderer *p) {
  if (!p)
    return;
  if (p->vbo)
    glDeleteBuffers(1, &p->vbo);
  if (p->vao)
    glDeleteVertexArrays(1, &p->vao);
  if (p->prog)
    glDeleteProgram(p->prog);
  memset(p, 0, sizeof(*p));
}

void particles_draw(ParticleRenderer *p, const Mat4 *world_to_clip, float time,
                    float cam_right_x, float cam_right_y, float cam_right_z,
                    float cam_up_x, float cam_up_y, float cam_up_z,
                    GLuint tex_array, const ParticleGPU *particles,
                    size_t particle_count) {
  if (!p || !p->prog || !world_to_clip || !particles || !particle_count)
    return;

  glUseProgram(p->prog);
  glActiveTexture(GL_TEXTURE0);
  glBindTexture(GL_TEXTURE_2D_ARRAY, tex_array);
  glUniform1f(p->uTime, time);
  glUniformMatrix4fv(p->uWorldToClip, 1, GL_FALSE, world_to_clip->m);
  glUniform3f(p->uCamRight, cam_right_x, cam_right_y, cam_right_z);
  glUniform3f(p->uCamUp, cam_up_x, cam_up_y, cam_up_z);
  glBindVertexArray(p->vao);
  glBindBuffer(GL_ARRAY_BUFFER, p->vbo);
  glBufferData(GL_ARRAY_BUFFER,
               (GLsizeiptr)(particle_count * sizeof(ParticleGPU)), particles,
               GL_STREAM_DRAW);

  // Depth pre-pass so nearer bodies properly occlude farther ones.
  glDisable(GL_BLEND);
  glColorMask(GL_FALSE, GL_FALSE, GL_FALSE, GL_FALSE);
  glDepthMask(GL_TRUE);
  glDepthFunc(GL_LESS);
  glDrawArraysInstanced(GL_TRIANGLE_STRIP, 0, 4, (GLsizei)particle_count);

  // Color pass: draw only the visible fragments.
  glEnable(GL_BLEND);
  glColorMask(GL_TRUE, GL_TRUE, GL_TRUE, GL_TRUE);
  glDepthMask(GL_FALSE);
  glDepthFunc(GL_EQUAL);
  glDrawArraysInstanced(GL_TRIANGLE_STRIP, 0, 4, (GLsizei)particle_count);

  // Restore defaults for other draws.
  glDepthFunc(GL_LESS);
  glDepthMask(GL_TRUE);
  glBindVertexArray(0);
  glUseProgram(0);
}

