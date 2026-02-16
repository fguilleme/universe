#include "glutil.h"

#include <stdio.h>
#include <stdlib.h>

GLuint glutil_compile_shader(GLenum type, const char *src) {
  GLuint sh = glCreateShader(type);
  glShaderSource(sh, 1, &src, NULL);
  glCompileShader(sh);

  GLint ok = 0;
  glGetShaderiv(sh, GL_COMPILE_STATUS, &ok);
  if (ok) return sh;

  GLint log_len = 0;
  glGetShaderiv(sh, GL_INFO_LOG_LENGTH, &log_len);
  char *log = (char *)calloc((size_t)log_len + 1, 1);
  if (log) glGetShaderInfoLog(sh, log_len, NULL, log);

  fprintf(stderr, "Shader compile failed (%u):\n%s\n", (unsigned)type, log ? log : "(no log)");
  free(log);
  glDeleteShader(sh);
  return 0;
}

GLuint glutil_link_program(GLuint vs, GLuint fs) {
  GLuint prog = glCreateProgram();
  glAttachShader(prog, vs);
  glAttachShader(prog, fs);
  glLinkProgram(prog);

  GLint ok = 0;
  glGetProgramiv(prog, GL_LINK_STATUS, &ok);
  if (ok) return prog;

  GLint log_len = 0;
  glGetProgramiv(prog, GL_INFO_LOG_LENGTH, &log_len);
  char *log = (char *)calloc((size_t)log_len + 1, 1);
  if (log) glGetProgramInfoLog(prog, log_len, NULL, log);

  fprintf(stderr, "Program link failed:\n%s\n", log ? log : "(no log)");
  free(log);
  glDeleteProgram(prog);
  return 0;
}

