#pragma once

#include <stdbool.h>
#include <stddef.h>

typedef struct MetalAccel MetalAccel;

// Creates a Metal compute context. Returns NULL if unavailable.
MetalAccel *metal_accel_create(void);
void metal_accel_destroy(MetalAccel *ma);

// Computes accelerations for N bodies.
// Inputs are in simulation units (AU, years, solar masses).
// Returns false if the compute path failed.
bool metal_accel_compute(MetalAccel *ma, const float *pos_mass4, size_t count, float G, float softening,
                         float *out_accel4);

