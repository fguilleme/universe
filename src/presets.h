#pragma once

#include <stdbool.h>

#include "camera.h"
#include "sim.h"

typedef enum {
  PRESET_TWO_BODY = 1,
  PRESET_GALAXY = 2,
  PRESET_SOLAR = 3,
} PresetId;

void preset_seed_two_body_orbit(Sim *sim);
void preset_seed_disk_galaxy(Sim *sim, int n, double disk_radius,
                             double central_mass);
void preset_seed_solar_system(Sim *sim);

// Loads preset configuration from `data/presets.yaml` (or UNIVERSE_PRESETS_YAML)
// and applies it to the simulation + camera.
bool preset_apply_from_yaml(const char *preset_name, Sim *sim, Camera *cam,
                            double *time_scale, uint32_t *selected_id,
                            bool *follow_selected);
