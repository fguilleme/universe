#include "input_mouse.h"

#include <math.h>
#include <stdlib.h>

#include "util.h"

void mouse_state_init(MouseState *ms) {
  if (!ms)
    return;
  *ms = (MouseState){0};
  ms->spawn_vel_scale = 1.0;
}

void mouse_handle_wheel(const SDL_MouseWheelEvent *wheel, Camera *cam) {
  if (!wheel || !cam)
    return;
  const double z = (wheel->y > 0) ? (1.0 / 1.15) : 1.15;
  cam->zoom_world_h = clampd(cam->zoom_world_h * z, 1e-3, 1e9);
}

void mouse_handle_button_down(MouseState *ms, SDL_Window *win,
                              const SDL_MouseButtonEvent *btn, Camera *cam,
                              Sim *sim) {
  (void)sim;
  if (!ms || !win || !btn || !cam)
    return;

  if (btn->button == SDL_BUTTON_MIDDLE) {
    int dx = 0, dy = 0;
    window_to_drawable(win, btn->x, btn->y, &dx, &dy);
    ms->panning = true;
    ms->pan_start_x = dx;
    ms->pan_start_y = dy;
    ms->pan_start_cx = cam->cx;
    ms->pan_start_cy = cam->cy;
    ms->pan_start_cz = cam->cz;
    return;
  }

  // LMB: rotate/select (selection on button-up if not rotated)
  if (btn->button == SDL_BUTTON_LEFT) {
    int dx = 0, dy = 0;
    window_to_drawable(win, btn->x, btn->y, &dx, &dy);
    ms->rotate_down = true;
    ms->rotating_cam = false;
    ms->rotate_start_x = dx;
    ms->rotate_start_y = dy;
    ms->rotate_start_yaw = cam->yaw;
    ms->rotate_start_pitch = cam->pitch;
    return;
  }

  // RMB: spawn
  if (btn->button == SDL_BUTTON_RIGHT) {
    int dx = 0, dy = 0;
    window_to_drawable(win, btn->x, btn->y, &dx, &dy);
    ms->spawning = true;
    ms->spawn_start_sx = dx;
    ms->spawn_start_sy = dy;
    int dw = 1, dh = 1;
    SDL_GL_GetDrawableSize(win, &dw, &dh);
    screen_to_world_on_z0(cam, dw, dh, ms->spawn_start_sx, ms->spawn_start_sy,
                          &ms->spawn_start_wx, &ms->spawn_start_wy,
                          &ms->spawn_start_wz);

    const uint16_t mod = SDL_GetModState();
    ms->spawn_vel_scale = (mod & KMOD_SHIFT) ? 5.0 : 1.0;
  }
}

static uint32_t pick_body(const Sim *sim, const Camera *cam, int dw, int dh,
                          int mdx, int mdy, bool realistic_scale,
                          double render_radius_scale) {
  Mat4 world_to_clip = camera_world_to_clip(cam, dw, dh);
  double best_d2_px = 0.0;
  uint32_t best_id = 0;
  for (size_t i = 0; i < sim->count; i++) {
    const Body *b = &sim->bodies[i];
    float clip[4];
    mat4_mul_vec4(&world_to_clip, (float)b->x, (float)b->y, (float)b->z, 1.0f,
                  clip);
    if (fabsf(clip[3]) < 1e-8f)
      continue;
    const float ndc_x = clip[0] / clip[3];
    const float ndc_y = clip[1] / clip[3];
    const float sx = (ndc_x * 0.5f + 0.5f) * (float)dw;
    const float sy = (0.5f - ndc_y * 0.5f) * (float)dh;

    const float dx = sx - (float)mdx;
    const float dy = sy - (float)mdy;
    const double d2 = (double)dx * (double)dx + (double)dy * (double)dy;
    const double radius_px_raw =
        body_visual_radius_world(b, realistic_scale, render_radius_scale) *
        ((double)dh / cam->zoom_world_h);
    const double pick_r_px = clampd(radius_px_raw * 1.4, 10.0, 60.0);
    if (d2 > pick_r_px * pick_r_px)
      continue;
    if (!best_id || d2 < best_d2_px) {
      best_d2_px = d2;
      best_id = b->id;
    }
  }
  if (best_id && best_d2_px <= 60.0 * 60.0)
    return best_id;
  return 0;
}

void mouse_handle_button_up(MouseState *ms, SDL_Window *win,
                            const SDL_MouseButtonEvent *btn, Camera *cam,
                            Sim *sim, bool realistic_scale,
                            double render_radius_scale, double spawn_mass,
                            uint32_t *selected_id, bool *follow_selected) {
  if (!ms || !win || !btn || !cam || !sim)
    return;

  if (btn->button == SDL_BUTTON_MIDDLE) {
    ms->panning = false;
    return;
  }

  if (btn->button == SDL_BUTTON_RIGHT) {
    if (ms->spawning) {
      int dw = 1, dh = 1;
      SDL_GL_GetDrawableSize(win, &dw, &dh);
      double end_wx = 0.0, end_wy = 0.0, end_wz = 0.0;
      int dx = 0, dy = 0;
      window_to_drawable(win, btn->x, btn->y, &dx, &dy);
      screen_to_world_on_z0(cam, dw, dh, dx, dy, &end_wx, &end_wy, &end_wz);

      Body b = {0};
      b.x = ms->spawn_start_wx;
      b.y = ms->spawn_start_wy;
      b.z = ms->spawn_start_wz;
      const uint16_t mod = SDL_GetModState();
      if (mod & KMOD_CTRL) {
        // Resolve center: selected body or heaviest body.
        const Body *center_body = NULL;
        if (selected_id && *selected_id) {
          for (size_t i = 0; i < sim->count; i++) {
            if (sim->bodies[i].id == *selected_id) {
              center_body = &sim->bodies[i];
              break;
            }
          }
        }
        if (!center_body) {
          for (size_t i = 0; i < sim->count; i++) {
            const Body *bb = &sim->bodies[i];
            if (!center_body || bb->mass > center_body->mass)
              center_body = bb;
          }
        }

        if (center_body) {
          const double rx = b.x - center_body->x;
          const double ry = b.y - center_body->y;
          const double rz = b.z - center_body->z;
          const double r = sqrt(rx * rx + ry * ry + rz * rz);
          if (r > 1e-6) {
            const double invr = 1.0 / r;
            double tx = -ry * invr;
            double ty = rx * invr;
            if (mod & KMOD_ALT) {
              tx = -tx;
              ty = -ty;
            }
            const double v_circ = sqrt(sim->G * center_body->mass / r);
            b.vx = center_body->vx + tx * v_circ;
            b.vy = center_body->vy + ty * v_circ;
            b.vz = center_body->vz;
          }
        }
      } else {
        b.vx = (end_wx - ms->spawn_start_wx) * ms->spawn_vel_scale;
        b.vy = (end_wy - ms->spawn_start_wy) * ms->spawn_vel_scale;
        b.vz = (end_wz - ms->spawn_start_wz) * ms->spawn_vel_scale;
      }
      b.mass = spawn_mass;
      (void)sim_add_body(sim, b);
    }
    ms->spawning = false;
    return;
  }

  if (btn->button == SDL_BUTTON_LEFT) {
    if (ms->rotate_down && !ms->rotating_cam) {
      int dw = 1, dh = 1;
      SDL_GL_GetDrawableSize(win, &dw, &dh);
      int mdx = 0, mdy = 0;
      window_to_drawable(win, btn->x, btn->y, &mdx, &mdy);
      const uint32_t best = pick_body(sim, cam, dw, dh, mdx, mdy,
                                      realistic_scale, render_radius_scale);
      if (selected_id)
        *selected_id = best;
      if (!best && follow_selected)
        *follow_selected = false;
    }
    ms->rotate_down = false;
    ms->rotating_cam = false;
  }
}

void mouse_handle_motion(MouseState *ms, SDL_Window *win,
                         const SDL_MouseMotionEvent *motion, Camera *cam) {
  if (!ms || !win || !motion || !cam)
    return;

  if (ms->panning) {
    int dw = 1, dh = 1;
    SDL_GL_GetDrawableSize(win, &dw, &dh);
    int mdx = 0, mdy = 0;
    window_to_drawable(win, motion->x, motion->y, &mdx, &mdy);
    const int dx = mdx - ms->pan_start_x;
    const int dy = mdy - ms->pan_start_y;
    const double aspect = (dh > 0) ? ((double)dw / (double)dh) : 1.0;
    const double world_per_px_y = cam->zoom_world_h / (double)dh;
    const double world_per_px_x = cam->zoom_world_h * aspect / (double)dw;

    const float vx = (float)(-(double)dx * world_per_px_x);
    const float vy = (float)((double)dy * world_per_px_y);
    float wx, wy, wz;
    view_to_world_dir(cam, vx, vy, 0.0f, &wx, &wy, &wz);
    cam->cx = ms->pan_start_cx + (double)wx;
    cam->cy = ms->pan_start_cy + (double)wy;
    cam->cz = ms->pan_start_cz + (double)wz;
  }

  if (ms->rotate_down) {
    int mdx = 0, mdy = 0;
    window_to_drawable(win, motion->x, motion->y, &mdx, &mdy);
    const int dx = mdx - ms->rotate_start_x;
    const int dy = mdy - ms->rotate_start_y;
    if (!ms->rotating_cam && (abs(dx) + abs(dy) > 3))
      ms->rotating_cam = true;
    if (ms->rotating_cam) {
      const double sens = 0.005;
      cam->yaw = ms->rotate_start_yaw + (double)dx * sens;
      cam->pitch = clampd(ms->rotate_start_pitch + (double)dy * sens, -1.55,
                          1.55);
    }
  }
}
