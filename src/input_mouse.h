#pragma once

#include <SDL2/SDL.h>

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

#include "camera.h"
#include "sim.h"

typedef struct {
  bool panning;
  int pan_start_x;
  int pan_start_y;
  double pan_start_cx;
  double pan_start_cy;
  double pan_start_cz;

  bool rotate_down;
  bool rotating_cam;
  int rotate_start_x;
  int rotate_start_y;
  double rotate_start_yaw;
  double rotate_start_pitch;

  bool spawning;
  int spawn_start_sx;
  int spawn_start_sy;
  double spawn_start_wx;
  double spawn_start_wy;
  double spawn_start_wz;
  double spawn_vel_scale;
} MouseState;

void mouse_state_init(MouseState *ms);

void mouse_handle_wheel(const SDL_MouseWheelEvent *wheel, Camera *cam);

void mouse_handle_button_down(MouseState *ms, SDL_Window *win,
                              const SDL_MouseButtonEvent *btn, Camera *cam,
                              Sim *sim);

void mouse_handle_button_up(MouseState *ms, SDL_Window *win,
                            const SDL_MouseButtonEvent *btn, Camera *cam,
                            Sim *sim, bool realistic_scale,
                            double render_radius_scale, double spawn_mass,
                            uint32_t *selected_id, bool *follow_selected);

void mouse_handle_motion(MouseState *ms, SDL_Window *win,
                         const SDL_MouseMotionEvent *motion, Camera *cam);

