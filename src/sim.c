#include "sim.h"

#include "metal_accel.h"

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <pthread.h>
#include <stdlib.h>
#include <string.h>

static double clampd(double x, double lo, double hi) {
  if (x < lo) return lo;
  if (x > hi) return hi;
  return x;
}

static double radius_from_mass_density(double mass, double density) {
  // 2D "radius" used for collisions/rendering; treat bodies as spheres with
  // volume = mass/density.
  if (mass <= 0.0 || density <= 0.0) return 1.0;
  const double volume = mass / density;
  const double r = cbrt((3.0 * volume) / (4.0 * M_PI));
  return fmax(r, 1e-6);
}

void sim_init(Sim *sim) {
  memset(sim, 0, sizeof(*sim));
  sim->G = 1.0;
  sim->softening = 0.01;
  sim->density = 1.0;
  sim->merge_on_collision = true;
  sim->next_id = 1;
  sim->num_threads = 1;
  sim->use_metal = false;
  sim->metal = (void *)metal_accel_create();
}

void sim_destroy(Sim *sim) {
  free(sim->bodies);
  free(sim->scratch_ax);
  free(sim->scratch_ay);
  free(sim->scratch_az);
  free(sim->metal_in4);
  free(sim->metal_out4);
  if (sim->metal) metal_accel_destroy((MetalAccel *)sim->metal);
  memset(sim, 0, sizeof(*sim));
}

void sim_enable_metal(Sim *sim, bool enabled) {
  sim->use_metal = enabled && (sim->metal != NULL);
}

static void sim_ensure_scratch(Sim *sim, size_t capacity) {
  if (sim->num_threads <= 1) return;
  if ((int)sim->scratch_threads != sim->num_threads || sim->scratch_capacity < capacity) {
    const size_t n = capacity * (size_t)sim->num_threads;
    double *ax = (double *)malloc(n * sizeof(double));
    double *ay = (double *)malloc(n * sizeof(double));
    double *az = (double *)malloc(n * sizeof(double));
    if (!ax || !ay || !az) {
      free(ax);
      free(ay);
      free(az);
      // If allocation fails, fall back to single-threaded.
      free(sim->scratch_ax);
      free(sim->scratch_ay);
      free(sim->scratch_az);
      sim->scratch_ax = sim->scratch_ay = sim->scratch_az = NULL;
      sim->scratch_capacity = 0;
      sim->scratch_threads = 0;
      sim->num_threads = 1;
      return;
    }
    free(sim->scratch_ax);
    free(sim->scratch_ay);
    free(sim->scratch_az);
    sim->scratch_ax = ax;
    sim->scratch_ay = ay;
    sim->scratch_az = az;
    sim->scratch_capacity = capacity;
    sim->scratch_threads = sim->num_threads;
  }
}

void sim_set_num_threads(Sim *sim, int threads) {
  if (threads < 1) threads = 1;
  sim->num_threads = threads;
  sim_ensure_scratch(sim, sim->capacity);
}

void sim_set_capacity(Sim *sim, size_t new_capacity) {
  if (new_capacity <= sim->capacity) return;
  Body *nb = (Body *)realloc(sim->bodies, new_capacity * sizeof(Body));
  if (!nb) return;
  sim->bodies = nb;
  sim->capacity = new_capacity;
  sim_ensure_scratch(sim, new_capacity);

  if (new_capacity > sim->metal_capacity) {
    float *in4 = (float *)realloc(sim->metal_in4, new_capacity * 4 * sizeof(float));
    if (in4) sim->metal_in4 = in4;
    float *out4 = (float *)realloc(sim->metal_out4, new_capacity * 4 * sizeof(float));
    if (out4) sim->metal_out4 = out4;
    if (sim->metal_in4 && sim->metal_out4) {
      sim->metal_capacity = new_capacity;
    }
  }
}

bool sim_add_body(Sim *sim, Body body) {
  if (sim->count == sim->capacity) {
    size_t next = sim->capacity ? sim->capacity * 2 : 256;
    sim_set_capacity(sim, next);
    if (sim->count == sim->capacity) return false;
  }

  body.alive = true;
  body.id = sim->next_id++;
  if (body.name[0] == '\0') {
    snprintf(body.name, sizeof(body.name), "Body %u", (unsigned)body.id);
  }
  if (!(body.mass > 0.0)) body.mass = 1.0;
  if (!(body.radius > 0.0)) body.radius = radius_from_mass_density(body.mass, sim->density);
  if (!(body.render_radius > 0.0)) body.render_radius = body.radius;
  // Defaults for visuals.
  // (tilt/spin are left as provided; default is 0)
  if (body.r == 0 && body.g == 0 && body.b == 0) {
    // Color by mass (roughly).
    const double t = clampd(log10(body.mass + 1.0) / 6.0, 0.0, 1.0);
    body.r = (float)(0.2 + 0.8 * t);
    body.g = (float)(0.4 + 0.4 * (1.0 - t));
    body.b = (float)(1.0 - 0.7 * t);
  }
  sim->bodies[sim->count++] = body;
  sim->have_accel = false;
  return true;
}

void sim_reset(Sim *sim) {
  sim->count = 0;
  sim->have_accel = false;
  sim->next_id = 1;
}

typedef struct {
  Sim *sim;
  size_t i0;
  size_t i1;
  double *ax;
  double *ay;
  double *az;
} AccelJob;

static void *compute_accel_worker(void *arg) {
  AccelJob *job = (AccelJob *)arg;
  Sim *sim = job->sim;
  const size_t n = sim->count;
  const double eps2 = sim->softening * sim->softening;

  for (size_t i = job->i0; i < job->i1; i++) {
    const Body *bi = &sim->bodies[i];
    if (!bi->alive) continue;
    for (size_t j = i + 1; j < n; j++) {
      const Body *bj = &sim->bodies[j];
      if (!bj->alive) continue;

      const double dx = bj->x - bi->x;
      const double dy = bj->y - bi->y;
      const double dz = bj->z - bi->z;
      const double r2 = dx * dx + dy * dy + dz * dz + eps2;
      const double inv_r = 1.0 / sqrt(r2);
      const double inv_r3 = inv_r * inv_r * inv_r;
      const double s = sim->G * inv_r3;

      const double aix = s * bj->mass * dx;
      const double aiy = s * bj->mass * dy;
      const double aiz = s * bj->mass * dz;
      const double ajx = -s * bi->mass * dx;
      const double ajy = -s * bi->mass * dy;
      const double ajz = -s * bi->mass * dz;

      job->ax[i] += aix;
      job->ay[i] += aiy;
      job->az[i] += aiz;
      job->ax[j] += ajx;
      job->ay[j] += ajy;
      job->az[j] += ajz;
    }
  }
  return NULL;
}

static void compute_accel(Sim *sim) {
  for (size_t i = 0; i < sim->count; i++) {
    sim->bodies[i].ax = 0.0;
    sim->bodies[i].ay = 0.0;
    sim->bodies[i].az = 0.0;
  }

  // Metal compute path (macOS). Uses float math for speed.
  if (sim->use_metal && sim->metal && sim->count >= 512 && sim->metal_capacity >= sim->count) {
    for (size_t i = 0; i < sim->count; i++) {
      sim->metal_in4[i * 4 + 0] = (float)sim->bodies[i].x;
      sim->metal_in4[i * 4 + 1] = (float)sim->bodies[i].y;
      sim->metal_in4[i * 4 + 2] = (float)sim->bodies[i].z;
      sim->metal_in4[i * 4 + 3] = (float)sim->bodies[i].mass;
    }
    const bool ok = metal_accel_compute((MetalAccel *)sim->metal, sim->metal_in4, sim->count, (float)sim->G,
                                        (float)sim->softening, sim->metal_out4);
    if (ok) {
      for (size_t i = 0; i < sim->count; i++) {
        sim->bodies[i].ax = (double)sim->metal_out4[i * 4 + 0];
        sim->bodies[i].ay = (double)sim->metal_out4[i * 4 + 1];
        sim->bodies[i].az = (double)sim->metal_out4[i * 4 + 2];
      }
      return;
    }
  }

  // For small N or single-thread mode, use the simple loop.
  if (sim->num_threads <= 1 || sim->count < 512) {
    const double eps2 = sim->softening * sim->softening;
    for (size_t i = 0; i < sim->count; i++) {
      Body *bi = &sim->bodies[i];
      if (!bi->alive) continue;
      for (size_t j = i + 1; j < sim->count; j++) {
        Body *bj = &sim->bodies[j];
        if (!bj->alive) continue;

        const double dx = bj->x - bi->x;
        const double dy = bj->y - bi->y;
        const double dz = bj->z - bi->z;
        const double r2 = dx * dx + dy * dy + dz * dz + eps2;
        const double inv_r = 1.0 / sqrt(r2);
        const double inv_r3 = inv_r * inv_r * inv_r;
        const double s = sim->G * inv_r3;

        const double aix = s * bj->mass * dx;
        const double aiy = s * bj->mass * dy;
        const double aiz = s * bj->mass * dz;
        const double ajx = -s * bi->mass * dx;
        const double ajy = -s * bi->mass * dy;
        const double ajz = -s * bi->mass * dz;

        bi->ax += aix;
        bi->ay += aiy;
        bi->az += aiz;
        bj->ax += ajx;
        bj->ay += ajy;
        bj->az += ajz;
      }
    }
    return;
  }

  sim_ensure_scratch(sim, sim->capacity);
  if (!sim->scratch_ax || !sim->scratch_ay || !sim->scratch_az) {
    sim->num_threads = 1;
    compute_accel(sim);
    return;
  }

  const int tcount = sim->num_threads;
  const size_t n = sim->count;
  for (int t = 0; t < tcount; t++) {
    memset(sim->scratch_ax + (size_t)t * sim->capacity, 0, n * sizeof(double));
    memset(sim->scratch_ay + (size_t)t * sim->capacity, 0, n * sizeof(double));
    memset(sim->scratch_az + (size_t)t * sim->capacity, 0, n * sizeof(double));
  }

  pthread_t *threads = (pthread_t *)malloc((size_t)tcount * sizeof(pthread_t));
  AccelJob *jobs = (AccelJob *)malloc((size_t)tcount * sizeof(AccelJob));
  if (!threads || !jobs) {
    free(threads);
    free(jobs);
    sim->num_threads = 1;
    compute_accel(sim);
    return;
  }

  for (int t = 0; t < tcount; t++) {
    const size_t i0 = (n * (size_t)t) / (size_t)tcount;
    const size_t i1 = (n * (size_t)(t + 1)) / (size_t)tcount;
    jobs[t] = (AccelJob){.sim = sim,
                         .i0 = i0,
                         .i1 = i1,
                         .ax = sim->scratch_ax + (size_t)t * sim->capacity,
                         .ay = sim->scratch_ay + (size_t)t * sim->capacity,
                         .az = sim->scratch_az + (size_t)t * sim->capacity};
    pthread_create(&threads[t], NULL, compute_accel_worker, &jobs[t]);
  }

  for (int t = 0; t < tcount; t++) {
    pthread_join(threads[t], NULL);
  }

  for (size_t i = 0; i < n; i++) {
    double ax = 0.0, ay = 0.0, az = 0.0;
    for (int t = 0; t < tcount; t++) {
      ax += sim->scratch_ax[(size_t)t * sim->capacity + i];
      ay += sim->scratch_ay[(size_t)t * sim->capacity + i];
      az += sim->scratch_az[(size_t)t * sim->capacity + i];
    }
    sim->bodies[i].ax = ax;
    sim->bodies[i].ay = ay;
    sim->bodies[i].az = az;
  }

  free(threads);
  free(jobs);
}

static void compact_dead(Sim *sim) {
  size_t out = 0;
  for (size_t i = 0; i < sim->count; i++) {
    if (!sim->bodies[i].alive) continue;
    if (out != i) sim->bodies[out] = sim->bodies[i];
    out++;
  }
  sim->count = out;
}

static void handle_collisions(Sim *sim) {
  if (!sim->merge_on_collision) return;

  for (size_t i = 0; i < sim->count; i++) {
    Body *a = &sim->bodies[i];
    if (!a->alive) continue;
    for (size_t j = i + 1; j < sim->count; j++) {
      Body *b = &sim->bodies[j];
      if (!b->alive) continue;

      const double dx = b->x - a->x;
      const double dy = b->y - a->y;
      const double dz = b->z - a->z;
      const double dist2 = dx * dx + dy * dy + dz * dz;
      const double r = a->radius + b->radius;
      if (dist2 > r * r) continue;

      // Merge: conserve momentum and mass.
      Body *big = a;
      Body *small = b;
      if (b->mass > a->mass) {
        big = b;
        small = a;
      }

      const double m1 = big->mass;
      const double m2 = small->mass;
      const double m = m1 + m2;
      if (!(m > 0.0)) continue;

      const double cx = (big->x * m1 + small->x * m2) / m;
      const double cy = (big->y * m1 + small->y * m2) / m;
      const double cz = (big->z * m1 + small->z * m2) / m;
      const double cvx = (big->vx * m1 + small->vx * m2) / m;
      const double cvy = (big->vy * m1 + small->vy * m2) / m;
      const double cvz = (big->vz * m1 + small->vz * m2) / m;

      big->x = cx;
      big->y = cy;
      big->z = cz;
      big->vx = cvx;
      big->vy = cvy;
      big->vz = cvz;
      big->mass = m;
      big->radius = radius_from_mass_density(m, sim->density);
      big->render_radius = radius_from_mass_density(m, sim->density);
      // Blend colors by mass.
      big->r = (float)((big->r * m1 + small->r * m2) / m);
      big->g = (float)((big->g * m1 + small->g * m2) / m);
      big->b = (float)((big->b * m1 + small->b * m2) / m);

      small->alive = false;
    }
  }
  compact_dead(sim);
}

void sim_step(Sim *sim, double dt) {
  if (sim->count == 0) return;
  if (!sim->have_accel) {
    compute_accel(sim);
    sim->have_accel = true;
  }

  // Velocity Verlet (symplectic-ish) integration.
  const double half_dt = 0.5 * dt;
  const double dt2_half = 0.5 * dt * dt;
  for (size_t i = 0; i < sim->count; i++) {
    Body *b = &sim->bodies[i];
    b->x += b->vx * dt + b->ax * dt2_half;
    b->y += b->vy * dt + b->ay * dt2_half;
    b->z += b->vz * dt + b->az * dt2_half;
    b->vx += b->ax * half_dt;
    b->vy += b->ay * half_dt;
    b->vz += b->az * half_dt;
  }

  handle_collisions(sim);
  compute_accel(sim);

  for (size_t i = 0; i < sim->count; i++) {
    Body *b = &sim->bodies[i];
    b->vx += b->ax * half_dt;
    b->vy += b->ay * half_dt;
    b->vz += b->az * half_dt;
  }
}
