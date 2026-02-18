#pragma once

#include <stdbool.h>

#include <SDL2/SDL.h>

#include "mat4.h"

typedef struct {
  double cx;
  double cy;
  double cz;
  double yaw;
  double pitch;
  double zoom_world_h;
} Camera;

bool world_to_screen_px(const Mat4 *world_to_clip, int w, int h, float x,
                        float y, float z, float *sx, float *sy);

void window_to_drawable(SDL_Window *win, int wx, int wy, int *dx, int *dy);

Mat4 camera_world_to_clip(const Camera *cam, int w, int h);

void view_to_world_dir(const Camera *cam, float vx, float vy, float vz,
                       float *wx, float *wy, float *wz);

void screen_to_world_on_z0(const Camera *cam, int w, int h, int sx, int sy,
                           double *wx, double *wy, double *wz);
