#pragma once

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

typedef struct {
  uint32_t id;
  char name[32];
  double x;
  double y;
  double z;
  double vx;
  double vy;
  double vz;
  double ax;
  double ay;
  double az;
  double mass;
  // Physical radius used for collisions.
  double radius;
  // Visual radius used for rendering.
  double render_radius;
  // Texture layer index for rendering (interpretation is renderer-specific).
  uint32_t tex_layer;
  // Visual spin parameters (do not affect physics).
  float tilt_rad;
  float spin_rate_rad_per_time;
  float spin_phase_rad;
  float r;
  float g;
  float b;
  bool alive;
} Body;

typedef struct {
  Body *bodies;
  size_t count;
  size_t capacity;

  uint32_t next_id;

  int num_threads;
  double *scratch_ax;
  double *scratch_ay;
  double *scratch_az;
  size_t scratch_capacity;
  int scratch_threads;

  // Optional Metal compute backend (macOS only).
  void *metal;
  bool use_metal;
  float *metal_in4;
  float *metal_out4;
  size_t metal_capacity;

  double G;
  double softening;
  double density;
  bool merge_on_collision;

  bool have_accel;
} Sim;

void sim_init(Sim *sim);
void sim_destroy(Sim *sim);

// Sets number of CPU threads used for gravity computation (>=1).
void sim_set_num_threads(Sim *sim, int threads);

void sim_reset(Sim *sim);
void sim_step(Sim *sim, double dt);

bool sim_add_body(Sim *sim, Body body);
void sim_set_capacity(Sim *sim, size_t new_capacity);

// Enables Metal compute acceleration on supported macOS systems.
void sim_enable_metal(Sim *sim, bool enabled);
