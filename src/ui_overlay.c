#include "ui_overlay.h"

#include <stdio.h>
#include <string.h>

#include "glutil.h"
#include "shaders.h"
#include "util.h"

void ui_text_renderer_draw(const UiTextRenderer *tr, const TextVertex *verts,
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

bool ui_text_renderer_init(UiTextRenderer *r) {
  if (!r)
    return false;
  memset(r, 0, sizeof(*r));

  GLuint vs = glutil_compile_shader(GL_VERTEX_SHADER, kTextVS);
  GLuint fs = glutil_compile_shader(GL_FRAGMENT_SHADER, kTextFS);
  if (!vs || !fs)
    return false;
  r->prog = glutil_link_program(vs, fs);
  glDeleteShader(vs);
  glDeleteShader(fs);
  if (!r->prog)
    return false;

  r->uViewport = glGetUniformLocation(r->prog, "uViewport");
  r->uOffsetPx = glGetUniformLocation(r->prog, "uOffsetPx");
  r->uFont = glGetUniformLocation(r->prog, "uFont");
  r->uColor = glGetUniformLocation(r->prog, "uColor");

  glGenVertexArrays(1, &r->vao);
  glGenBuffers(1, &r->vbo);
  glBindVertexArray(r->vao);
  glBindBuffer(GL_ARRAY_BUFFER, r->vbo);
  glBufferData(GL_ARRAY_BUFFER, 1, NULL, GL_STREAM_DRAW);
  glEnableVertexAttribArray(0);
  glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, sizeof(TextVertex),
                        (void *)offsetof(TextVertex, x));
  glEnableVertexAttribArray(1);
  glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, sizeof(TextVertex),
                        (void *)offsetof(TextVertex, u));
  glBindVertexArray(0);

  r->font_tex = create_font_texture();
  if (!r->font_tex)
    return false;

  glUseProgram(r->prog);
  glUniform1i(r->uFont, 2);
  glUseProgram(0);
  return true;
}

void ui_text_renderer_destroy(UiTextRenderer *r) {
  if (!r)
    return;
  if (r->font_tex)
    glDeleteTextures(1, &r->font_tex);
  if (r->vbo)
    glDeleteBuffers(1, &r->vbo);
  if (r->vao)
    glDeleteVertexArrays(1, &r->vao);
  if (r->prog)
    glDeleteProgram(r->prog);
  memset(r, 0, sizeof(*r));
}

void ui_overlay_draw(const UiTextRenderer *tr, const UiOverlayConfig *cfg,
                     const UiOverlayState *st, int dw, int dh) {
  if (!tr || !tr->prog || !cfg || !st)
    return;

  glDisable(GL_DEPTH_TEST);
  glDepthMask(GL_FALSE);

  TextVertex *verts = NULL;
  size_t count = 0;
  size_t cap = 0;

  // Date/time.
  if (st->preset == 3) {
    const double jd = cfg->solar_epoch_jd + st->sim_time * 365.25;
    int yy = 0, mo = 0, dd = 0, hh = 0, mm = 0, ss = 0;
    jd_to_gregorian_utc(jd, &yy, &mo, &dd, &hh, &mm, &ss);
    char dt[64];
    snprintf(dt, sizeof(dt), "%04d-%02d-%02d %02d:%02d:%02d", yy, mo, dd, hh,
             mm, ss);
    const float scale = 3.0f;
    const float text_w = (float)strlen(dt) * 8.0f * scale;
    const float x = (float)clampd(((double)dw - (double)text_w) * 0.5, 4.0,
                                  (double)dw - 4.0 - text_w);
    const float y = 10.0f;
    text_append(&verts, &count, &cap, x, y, scale, dt);
    ui_text_renderer_draw(tr, verts, count, dw, dh, 2.0f, 2.0f, 0, 0, 0,
                          0.65f);
    ui_text_renderer_draw(tr, verts, count, dw, dh, 0.0f, 0.0f, 1, 1, 1,
                          0.95f);
    count = 0;
  }

  // Help.
  if (st->show_help) {
    const float scale = 2.0f;
    const float x0 = 12.0f;
    float y0 = 34.0f;
    const char *lines[] = {
        "F1: toggle help",
        "F2: toggle body labels",
        "F3: realistic scaling",
        "",
        "RMB drag: spawn body (Shift faster)",
        "  Ctrl: circular orbit around selection",
        "  Alt: reverse orbit direction",
        "LMB click: select body",
        "LMB drag: rotate camera",
        "MMB drag / WASD: pan",
        "Wheel: zoom",
        "Z: zoom to selected",
        "F: follow selected",
        "T: trails   V: velocity vectors",
        "PgUp/PgDn: trail length",
        "Space: pause   .: step",
        "- / =: time scale",
        "G / H: gravity",
        "1: two-body   2: galaxy   3: solar system",
        "R: reset solar system   C: clear",
        "Esc: quit",
        NULL,
    };
    for (int i = 0; lines[i]; i++) {
      text_append(&verts, &count, &cap, x0, y0, scale, lines[i]);
      y0 += 10.0f * scale;
    }
    ui_text_renderer_draw(tr, verts, count, dw, dh, 1.0f, 1.0f, 0, 0, 0,
                          0.55f);
    ui_text_renderer_draw(tr, verts, count, dw, dh, 0.0f, 0.0f, 1, 1, 1,
                          0.88f);
    count = 0;
  }

  // Camera info overlay.
  {
    char cam_txt[192];
    snprintf(cam_txt, sizeof(cam_txt),
             "cam: (%.3f %.3f %.3f)  yaw=%.3f  pitch=%.3f  zoom=%.3f",
             st->cam.cx, st->cam.cy, st->cam.cz, st->cam.yaw, st->cam.pitch,
             st->cam.zoom_world_h);
    const float scale = 2.5f;
    const float text_w = (float)strlen(cam_txt) * 8.0f * scale;
    const float x = (float)clampd(((double)dw - (double)text_w) * 0.5, 4.0,
                                  (double)dw - 4.0 - text_w);
    const float y = (float)dh - (8.0f * scale) - 8.0f - 22.0f;
    text_append(&verts, &count, &cap, x, y, scale, cam_txt);
    ui_text_renderer_draw(tr, verts, count, dw, dh, 2.0f, 2.0f, 0, 0, 0,
                          0.55f);
    ui_text_renderer_draw(tr, verts, count, dw, dh, 0.0f, 0.0f, 1, 1, 1,
                          0.85f);
    count = 0;
  }

  // Depth debug.
  {
    char dbg[128];
    const int scene_bits =
        st->use_offscreen_scene_fbo ? st->scene_depth_bits : st->default_depth_bits;
    snprintf(dbg, sizeof(dbg), "Depth: default=%d scene=%d%s",
             st->default_depth_bits, scene_bits,
             st->use_offscreen_scene_fbo ? " (FBO)" : "");
    const float scale = 2.0f;
    const float x = 12.0f;
    const float y = (float)dh - (8.0f * scale) - 6.0f;
    text_append(&verts, &count, &cap, x, y, scale, dbg);
    ui_text_renderer_draw(tr, verts, count, dw, dh, 1.0f, 1.0f, 0, 0, 0,
                          0.55f);
    ui_text_renderer_draw(tr, verts, count, dw, dh, 0.0f, 0.0f, 1, 1, 1,
                          0.75f);
    count = 0;
  }

  // Scaling mode.
  {
    const char *mode = st->realistic_scale ? "Scale: realistic (F3)"
                                           : "Scale: exaggerated (F3)";
    const float scale = 2.0f;
    const float text_w = (float)strlen(mode) * 8.0f * scale;
    const float x = (float)clampd(((double)dw - (double)text_w) * 0.5, 4.0,
                                  (double)dw - 4.0 - text_w);
    const float y = (float)dh - (8.0f * scale) - 6.0f;
    text_append(&verts, &count, &cap, x, y, scale, mode);
    ui_text_renderer_draw(tr, verts, count, dw, dh, 1.5f, 1.5f, 0, 0, 0,
                          0.55f);
    ui_text_renderer_draw(tr, verts, count, dw, dh, 0.0f, 0.0f, 1, 1, 1,
                          0.80f);
    count = 0;
  }

  free(verts);
  glDepthMask(GL_TRUE);
  glEnable(GL_DEPTH_TEST);
}
