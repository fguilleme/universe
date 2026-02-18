#pragma once

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

