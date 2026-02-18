#pragma once

typedef struct {
  float m[16];
} Mat4;

Mat4 mat4_identity(void);
Mat4 mat4_mul(Mat4 a, Mat4 b);
Mat4 mat4_translate(float x, float y, float z);
Mat4 mat4_rot_x(float a);
Mat4 mat4_rot_z(float a);
Mat4 mat4_ortho(float left, float right, float bottom, float top, float near_z,
                float far_z);
Mat4 mat4_inset_world_to_clip(int w, int h, double cx, double cy,
                              double half_world);

void mat4_mul_vec4(const Mat4 *m, float x, float y, float z, float w,
                   float out4[4]);

