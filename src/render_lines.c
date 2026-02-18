#include "render_lines.h"

#include <stddef.h>
#include <string.h>

#include "glutil.h"
#include "shaders.h"

bool lines_init(LineRenderer *r) {
  if (!r)
    return false;
  memset(r, 0, sizeof(*r));

  GLuint lvs = glutil_compile_shader(GL_VERTEX_SHADER, kLinesVS);
  GLuint lfs = glutil_compile_shader(GL_FRAGMENT_SHADER, kLinesFS);
  if (!lvs || !lfs)
    return false;
  r->prog = glutil_link_program(lvs, lfs);
  glDeleteShader(lvs);
  glDeleteShader(lfs);
  if (!r->prog)
    return false;
  r->uWorldToClip = glGetUniformLocation(r->prog, "uWorldToClip");

  glGenVertexArrays(1, &r->vao);
  glGenBuffers(1, &r->vbo);
  glBindVertexArray(r->vao);
  glBindBuffer(GL_ARRAY_BUFFER, r->vbo);
  glBufferData(GL_ARRAY_BUFFER, 1, NULL, GL_STREAM_DRAW);
  glEnableVertexAttribArray(0);
  glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(LineVertex),
                        (void *)offsetof(LineVertex, x));
  glEnableVertexAttribArray(1);
  glVertexAttribPointer(1, 4, GL_FLOAT, GL_FALSE, sizeof(LineVertex),
                        (void *)offsetof(LineVertex, r));
  glBindVertexArray(0);
  return true;
}

void lines_destroy(LineRenderer *r) {
  if (!r)
    return;
  if (r->vbo)
    glDeleteBuffers(1, &r->vbo);
  if (r->vao)
    glDeleteVertexArrays(1, &r->vao);
  if (r->prog)
    glDeleteProgram(r->prog);
  memset(r, 0, sizeof(*r));
}

void lines_draw(LineRenderer *r, const Mat4 *world_to_clip,
                const LineVertex *verts, size_t count) {
  if (!r || !r->prog || !world_to_clip || !verts || !count)
    return;
  glUseProgram(r->prog);
  glUniformMatrix4fv(r->uWorldToClip, 1, GL_FALSE, world_to_clip->m);
  glBindVertexArray(r->vao);
  glBindBuffer(GL_ARRAY_BUFFER, r->vbo);
  glBufferData(GL_ARRAY_BUFFER, (GLsizeiptr)(count * sizeof(LineVertex)), verts,
               GL_STREAM_DRAW);
  glDrawArrays(GL_LINES, 0, (GLsizei)count);
  glBindVertexArray(0);
  glUseProgram(0);
}

