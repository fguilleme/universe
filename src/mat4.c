#include "mat4.h"

#include <math.h>

Mat4 mat4_identity(void) {
  Mat4 out = {.m = {0}};
  out.m[0] = 1.0f;
  out.m[5] = 1.0f;
  out.m[10] = 1.0f;
  out.m[15] = 1.0f;
  return out;
}

Mat4 mat4_mul(Mat4 a, Mat4 b) {
  Mat4 out = {.m = {0}};
  for (int c = 0; c < 4; c++) {
    for (int r = 0; r < 4; r++) {
      out.m[c * 4 + r] =
          a.m[0 * 4 + r] * b.m[c * 4 + 0] + a.m[1 * 4 + r] * b.m[c * 4 + 1] +
          a.m[2 * 4 + r] * b.m[c * 4 + 2] + a.m[3 * 4 + r] * b.m[c * 4 + 3];
    }
  }
  return out;
}

Mat4 mat4_translate(float x, float y, float z) {
  Mat4 out = mat4_identity();
  out.m[12] = x;
  out.m[13] = y;
  out.m[14] = z;
  return out;
}

Mat4 mat4_rot_x(float a) {
  Mat4 out = mat4_identity();
  const float c = cosf(a);
  const float s = sinf(a);
  out.m[5] = c;
  out.m[6] = s;
  out.m[9] = -s;
  out.m[10] = c;
  return out;
}

Mat4 mat4_rot_z(float a) {
  Mat4 out = mat4_identity();
  const float c = cosf(a);
  const float s = sinf(a);
  out.m[0] = c;
  out.m[1] = s;
  out.m[4] = -s;
  out.m[5] = c;
  return out;
}

Mat4 mat4_ortho(float left, float right, float bottom, float top, float near_z,
                float far_z) {
  Mat4 out = {.m = {0}};
  out.m[0] = 2.0f / (right - left);
  out.m[5] = 2.0f / (top - bottom);
  out.m[10] = -2.0f / (far_z - near_z);
  out.m[12] = -(right + left) / (right - left);
  out.m[13] = -(top + bottom) / (top - bottom);
  out.m[14] = -(far_z + near_z) / (far_z - near_z);
  out.m[15] = 1.0f;
  return out;
}

Mat4 mat4_inset_world_to_clip(int w, int h, double cx, double cy,
                              double half_world) {
  const double aspect = (h > 0) ? ((double)w / (double)h) : 1.0;
  const float half_h = (float)half_world;
  const float half_w = (float)(half_world * aspect);
  Mat4 proj = mat4_ortho(-half_w, half_w, -half_h, half_h, -1.0f, 1.0f);
  Mat4 view = mat4_translate((float)-cx, (float)-cy, 0.0f);
  return mat4_mul(proj, view);
}

void mat4_mul_vec4(const Mat4 *m, float x, float y, float z, float w,
                   float out4[4]) {
  out4[0] = m->m[0] * x + m->m[4] * y + m->m[8] * z + m->m[12] * w;
  out4[1] = m->m[1] * x + m->m[5] * y + m->m[9] * z + m->m[13] * w;
  out4[2] = m->m[2] * x + m->m[6] * y + m->m[10] * z + m->m[14] * w;
  out4[3] = m->m[3] * x + m->m[7] * y + m->m[11] * z + m->m[15] * w;
}

