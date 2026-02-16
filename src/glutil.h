#pragma once

#include <GL/glew.h>

GLuint glutil_compile_shader(GLenum type, const char *src);
GLuint glutil_link_program(GLuint vs, GLuint fs);

