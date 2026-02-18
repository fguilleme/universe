#include "camera.h"

#include <math.h>

#include "mat4.h"

bool world_to_screen_px(const Mat4 *world_to_clip, int w, int h, float x,
                        float y, float z, float *sx, float *sy) {
  float clip[4];
  mat4_mul_vec4(world_to_clip, x, y, z, 1.0f, clip);
  if (fabsf(clip[3]) < 1e-8f)
    return false;
  const float ndc_x = clip[0] / clip[3];
  const float ndc_y = clip[1] / clip[3];
  *sx = (ndc_x * 0.5f + 0.5f) * (float)w;
  *sy = (0.5f - ndc_y * 0.5f) * (float)h;
  return true;
}

void window_to_drawable(SDL_Window *win, int wx, int wy, int *dx, int *dy) {
  int ww = 1, wh = 1;
  int dw = 1, dh = 1;
  SDL_GetWindowSize(win, &ww, &wh);
  SDL_GL_GetDrawableSize(win, &dw, &dh);
  const float sx = (ww > 0) ? ((float)dw / (float)ww) : 1.0f;
  const float sy = (wh > 0) ? ((float)dh / (float)wh) : 1.0f;
  *dx = (int)lrintf((float)wx * sx);
  *dy = (int)lrintf((float)wy * sy);
}

Mat4 camera_world_to_clip(const Camera *cam, int w, int h) {
  const double aspect = (h > 0) ? ((double)w / (double)h) : 1.0;
  const float half_h = (float)(cam->zoom_world_h * 0.5);
  const float half_w = (float)(cam->zoom_world_h * aspect * 0.5);
  // Keep the depth range reasonably tight for better depth precision,
  // otherwise small Z differences (e.g. planet behind the Sun) quantize away.
  const float depth_range = (float)fmax(1000.0, cam->zoom_world_h * 100.0);
  Mat4 proj =
      mat4_ortho(-half_w, half_w, -half_h, half_h, -depth_range, depth_range);

  Mat4 t = mat4_translate((float)-cam->cx, (float)-cam->cy, (float)-cam->cz);
  Mat4 rz = mat4_rot_z((float)cam->yaw);
  Mat4 rx = mat4_rot_x((float)cam->pitch);
  Mat4 view = mat4_mul(rx, mat4_mul(rz, t));
  return mat4_mul(proj, view);
}

void view_to_world_dir(const Camera *cam, float vx, float vy, float vz,
                       float *wx, float *wy, float *wz) {
  // Inverse of view rotation: world = Rz(-yaw)*Rx(-pitch)*view
  const float cp = cosf((float)-cam->pitch);
  const float sp = sinf((float)-cam->pitch);
  float x1 = vx;
  float y1 = cp * vy + -sp * vz;
  float z1 = sp * vy + cp * vz;

  const float cy = cosf((float)-cam->yaw);
  const float sy = sinf((float)-cam->yaw);
  *wx = cy * x1 + -sy * y1;
  *wy = sy * x1 + cy * y1;
  *wz = z1;
}

void screen_to_world_on_z0(const Camera *cam, int w, int h, int sx, int sy,
                           double *wx, double *wy, double *wz) {
  const double aspect = (h > 0) ? ((double)w / (double)h) : 1.0;
  const double x_ndc = 2.0 * ((double)sx / (double)w) - 1.0;
  const double y_ndc = 1.0 - 2.0 * ((double)sy / (double)h);
  const float half_h = (float)(cam->zoom_world_h * 0.5);
  const float half_w = (float)(cam->zoom_world_h * aspect * 0.5);

  // View-space position on the z=0 slice.
  const float vx = (float)(x_ndc * half_w);
  const float vy = (float)(y_ndc * half_h);

  // World-space ray: origin at that view point, direction is +viewZ.
  float ox, oy, oz;
  view_to_world_dir(cam, vx, vy, 0.0f, &ox, &oy, &oz);
  ox += (float)cam->cx;
  oy += (float)cam->cy;
  oz += (float)cam->cz;

  float dx, dy, dz;
  view_to_world_dir(cam, 0.0f, 0.0f, 1.0f, &dx, &dy, &dz);
  if (fabsf(dz) < 1e-8f) {
    *wx = ox;
    *wy = oy;
    *wz = 0.0;
    return;
  }
  const float t = -oz / dz;
  *wx = ox + dx * t;
  *wy = oy + dy * t;
  *wz = 0.0;
}

