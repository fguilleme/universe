#pragma once

#include <stdbool.h>

#include "camera.h"
#include "mat4.h"
#include "sim.h"
#include "ui_overlay.h"
#include "ui_text.h"

void ui_draw_selected_hud(const UiTextRenderer *tr, TextVertex **verts,
                          size_t *cap, int dw, int dh,
                          const Mat4 *world_to_clip, const Camera *cam,
                          const Body *sel, bool realistic_scale,
                          double render_radius_scale);

void ui_draw_body_labels(const UiTextRenderer *tr, TextVertex **verts,
                         size_t *cap, int dw, int dh,
                         const Mat4 *world_to_clip, const Camera *cam,
                         const Sim *sim, bool realistic_scale,
                         double render_radius_scale);

