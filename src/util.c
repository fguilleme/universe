#include "util.h"

#include <math.h>
#include <stdio.h>
#include <string.h>

double clampd(double x, double lo, double hi) {
  if (x < lo)
    return lo;
  if (x > hi)
    return hi;
  return x;
}

uint32_t hash_u32(uint32_t x) {
  // Cheap integer hash.
  x ^= x >> 16;
  x *= 0x7feb352du;
  x ^= x >> 15;
  x *= 0x846ca68bu;
  x ^= x >> 16;
  return x;
}

void hsv_to_rgb(float h, float s, float v, float *r, float *g, float *b) {
  const float hh = h - floorf(h);
  const float c = v * s;
  const float x = c * (1.0f - fabsf(fmodf(hh * 6.0f, 2.0f) - 1.0f));
  const float m = v - c;
  float rr = 0.0f, gg = 0.0f, bb = 0.0f;
  const float k = hh * 6.0f;
  if (k < 1.0f) {
    rr = c;
    gg = x;
  } else if (k < 2.0f) {
    rr = x;
    gg = c;
  } else if (k < 3.0f) {
    gg = c;
    bb = x;
  } else if (k < 4.0f) {
    gg = x;
    bb = c;
  } else if (k < 5.0f) {
    rr = x;
    bb = c;
  } else {
    rr = c;
    bb = x;
  }
  *r = rr + m;
  *g = gg + m;
  *b = bb + m;
}

void trail_color_from_id(uint32_t id, float *r, float *g, float *b) {
  const uint32_t h = hash_u32(id ? id : 1u);
  const float hue = (float)(h & 0xFFFFu) / 65535.0f;
  hsv_to_rgb(hue, 0.75f, 0.95f, r, g, b);
}

void str_to_lower_ascii(char *dst, size_t cap, const char *src) {
  if (!cap)
    return;
  size_t i = 0;
  for (; src[i] && i + 1 < cap; i++) {
    char c = src[i];
    if (c >= 'A' && c <= 'Z')
      c = (char)(c - 'A' + 'a');
    dst[i] = c;
  }
  dst[i] = 0;
}

void jd_to_gregorian_utc(double jd, int *year, int *month, int *day, int *hour,
                         int *min, int *sec) {
  // Astronomical Julian Day -> proleptic Gregorian calendar.
  // This is sufficient for an on-screen simulation clock.
  const double jd0 = jd + 0.5;
  long z = (long)floor(jd0);
  double f = jd0 - (double)z;

  long a = z;
  if (z >= 2299161) {
    const long alpha = (long)floor(((double)z - 1867216.25) / 36524.25);
    a = z + 1 + alpha - alpha / 4;
  }
  const long b = a + 1524;
  const long c = (long)floor(((double)b - 122.1) / 365.25);
  const long d = (long)floor(365.25 * (double)c);
  const long e = (long)floor(((double)b - (double)d) / 30.6001);

  const double day_f = (double)b - (double)d - floor(30.6001 * (double)e) + f;
  const int day_i = (int)floor(day_f);
  const double frac = day_f - (double)day_i;

  int m = (e < 14) ? (int)(e - 1) : (int)(e - 13);
  int y = (m > 2) ? (int)(c - 4716) : (int)(c - 4715);

  const double h_f = frac * 24.0;
  int h = (int)floor(h_f);
  const double min_f = (h_f - (double)h) * 60.0;
  int mi = (int)floor(min_f);
  int s = (int)lrint((min_f - (double)mi) * 60.0);
  if (s >= 60) {
    s = 0;
    mi++;
  }
  if (mi >= 60) {
    mi = 0;
    h++;
  }
  if (h >= 24) {
    // Rare (rounding) overflow; keep it simple.
    h = 0;
  }

  *year = y;
  *month = m;
  *day = day_i;
  *hour = h;
  *min = mi;
  *sec = s;
}

double body_visual_radius_world(const Body *b, bool realistic_scale,
                                double exaggerated_scale) {
  const double base = realistic_scale ? b->radius : b->render_radius;
  const double scale = realistic_scale ? 1.0 : exaggerated_scale;
  return base * scale;
}
