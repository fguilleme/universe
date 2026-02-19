#pragma once

#include <stdbool.h>

#include "sim.h"

// Dumps the current simulation bodies to a YAML file at `path`.
// Returns true on success.
bool dump_bodies_yaml(const char *path, const Sim *sim, double sim_time);

