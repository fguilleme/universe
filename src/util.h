#pragma once

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

#include "sim.h"

double clampd(double x, double lo, double hi);

uint32_t hash_u32(uint32_t x);
void hsv_to_rgb(float h, float s, float v, float *r, float *g, float *b);
void trail_color_from_id(uint32_t id, float *r, float *g, float *b);

void str_to_lower_ascii(char *dst, size_t cap, const char *src);

void jd_to_gregorian_utc(double jd, int *year, int *month, int *day, int *hour,
                         int *min, int *sec);

double body_visual_radius_world(const Body *b, bool realistic_scale,
                                double exaggerated_scale);
