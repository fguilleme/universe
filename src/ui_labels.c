#include "ui_labels.h"

#include <stdio.h>
#include <string.h>

#include "util.h"

static void ui_text_draw(const UiTextRenderer *tr, TextVertex *verts,
                         size_t count, int dw, int dh, float offx, float offy,
                         float r, float g, float b, float a) {
  if (!count)
    return;
  glBindVertexArray(tr->vao);
  glBindBuffer(GL_ARRAY_BUFFER, tr->vbo);
  glBufferData(GL_ARRAY_BUFFER, (GLsizeiptr)(count * sizeof(TextVertex)), verts,
               GL_STREAM_DRAW);

  glUseProgram(tr->prog);
  glUniform2f(tr->uViewport, (float)dw, (float)dh);
  glActiveTexture(GL_TEXTURE2);
  glBindTexture(GL_TEXTURE_2D, tr->font_tex);
  glUniform2f(tr->uOffsetPx, offx, offy);
  glUniform4f(tr->uColor, r, g, b, a);
  glDrawArrays(GL_TRIANGLES, 0, (GLsizei)count);
  glUseProgram(0);
  glBindVertexArray(0);
}

void ui_draw_selected_hud(const UiTextRenderer *tr, TextVertex **verts,
                          size_t *cap, int dw, int dh,
                          const Mat4 *world_to_clip, const Camera *cam,
                          const Body *sel, bool realistic_scale,
                          double render_radius_scale) {
  if (!tr || !verts || !cap || !world_to_clip || !cam || !sel)
    return;

  glDisable(GL_DEPTH_TEST);
  glDepthMask(GL_FALSE);

  // Selected: HUD text.
  {
    char hud[256];
    snprintf(hud, sizeof(hud), "Selected: %s", sel->name);
    size_t count = 0;
    text_append(verts, &count, cap, 12.0f, 12.0f, 1.0f, hud);
    ui_text_draw(tr, *verts, count, dw, dh, 1.0f, 1.0f, 0, 0, 0, 0.55f);
    ui_text_draw(tr, *verts, count, dw, dh, 0.0f, 0.0f, 1, 1, 1, 0.9f);
  }

  // Label above the selected body.
  {
    float sx = 0.0f, sy = 0.0f;
    if (world_to_screen_px(world_to_clip, dw, dh, (float)sel->x, (float)sel->y,
                           (float)sel->z, &sx, &sy)) {
      const float scale = 1.5f;
      float radius_px =
          (float)(body_visual_radius_world(sel, realistic_scale,
                                           render_radius_scale) *
                  (double)dh / cam->zoom_world_h);
      radius_px = (float)clampd(radius_px, 6.0, 70.0);
      const float text_w = (float)strlen(sel->name) * 8.0f * scale;
      const float x =
          (float)clampd((double)(sx - text_w * 0.5f), 4.0,
                        (double)dw - 4.0 - (double)text_w);
      const float y =
          (float)clampd((double)(sy - radius_px - 18.0f * scale), 4.0,
                        (double)dh - 20.0);

      size_t count = 0;
      text_append(verts, &count, cap, x, y, scale, sel->name);
      ui_text_draw(tr, *verts, count, dw, dh, 1.0f, 1.0f, 0, 0, 0, 0.6f);
      ui_text_draw(tr, *verts, count, dw, dh, 0.0f, 0.0f, 1, 1, 1, 0.95f);
    }
  }

  glDepthMask(GL_TRUE);
  glEnable(GL_DEPTH_TEST);
}

void ui_draw_body_labels(const UiTextRenderer *tr, TextVertex **verts,
                         size_t *cap, int dw, int dh,
                         const Mat4 *world_to_clip, const Camera *cam,
                         const Sim *sim, bool realistic_scale,
                         double render_radius_scale) {
  if (!tr || !verts || !cap || !world_to_clip || !cam || !sim)
    return;

  glDisable(GL_DEPTH_TEST);
  glDepthMask(GL_FALSE);

  size_t count = 0;
  for (size_t i = 0; i < sim->count; i++) {
    const Body *b = &sim->bodies[i];
    float sx = 0.0f, sy = 0.0f;
    if (!world_to_screen_px(world_to_clip, dw, dh, (float)b->x, (float)b->y,
                            (float)b->z, &sx, &sy))
      continue;

    const float scale = 1.0f;
    float radius_px =
        (float)(body_visual_radius_world(b, realistic_scale,
                                         render_radius_scale) *
                (double)dh / cam->zoom_world_h);
    radius_px = (float)clampd(radius_px, 6.0, 60.0);
    const float text_w = (float)strlen(b->name) * 8.0f * scale;
    const float x =
        (float)clampd((double)(sx - text_w * 0.5f), 4.0,
                      (double)dw - 4.0 - (double)text_w);
    const float y =
        (float)clampd((double)(sy - radius_px - 16.0f * scale), 4.0,
                      (double)dh - 20.0);
    text_append(verts, &count, cap, x, y, scale, b->name);
  }

  ui_text_draw(tr, *verts, count, dw, dh, 1.0f, 1.0f, 0, 0, 0, 0.55f);
  ui_text_draw(tr, *verts, count, dw, dh, 0.0f, 0.0f, 1, 1, 1, 0.9f);

  glDepthMask(GL_TRUE);
  glEnable(GL_DEPTH_TEST);
}

