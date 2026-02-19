#include "dump_yaml.h"

#include <inttypes.h>
#include <stdio.h>
#include <string.h>

static void yaml_write_quoted_string(FILE *f, const char *s) {
  fputc('\'', f);
  if (s) {
    for (const char *p = s; *p; p++) {
      if (*p == '\'') {
        fputc('\'', f);
        fputc('\'', f);
      } else {
        fputc(*p, f);
      }
    }
  }
  fputc('\'', f);
}

bool dump_bodies_yaml(const char *path, const Sim *sim, double sim_time) {
  if (!path || !path[0] || !sim)
    return false;

  FILE *f = fopen(path, "wb");
  if (!f)
    return false;

  fprintf(f, "version: 1\n");
  fprintf(f, "sim_time: %.17g\n", sim_time);
  fprintf(f, "G: %.17g\n", sim->G);
  fprintf(f, "softening: %.17g\n", sim->softening);
  fprintf(f, "density: %.17g\n", sim->density);
  fprintf(f, "merge_on_collision: %s\n", sim->merge_on_collision ? "true" : "false");
  fprintf(f, "count: %zu\n", sim->count);
  fprintf(f, "bodies:\n");

  for (size_t i = 0; i < sim->count; i++) {
    const Body *b = &sim->bodies[i];
    fprintf(f, "  - id: %" PRIu32 "\n", b->id);
    fprintf(f, "    name: ");
    yaml_write_quoted_string(f, b->name);
    fputc('\n', f);
    fprintf(f, "    alive: %s\n", b->alive ? "true" : "false");
    fprintf(f, "    mass: %.17g\n", b->mass);

    fprintf(f, "    pos: [%.17g, %.17g, %.17g]\n", b->x, b->y, b->z);
    fprintf(f, "    vel: [%.17g, %.17g, %.17g]\n", b->vx, b->vy, b->vz);
    fprintf(f, "    accel: [%.17g, %.17g, %.17g]\n", b->ax, b->ay, b->az);

    fprintf(f, "    radius: %.17g\n", b->radius);
    fprintf(f, "    render_radius: %.17g\n", b->render_radius);
    fprintf(f, "    tex_layer: %" PRIu32 "\n", b->tex_layer);

    fprintf(f, "    tilt_rad: %.9g\n", (double)b->tilt_rad);
    fprintf(f, "    spin_rate_rad_per_time: %.9g\n", (double)b->spin_rate_rad_per_time);
    fprintf(f, "    spin_phase_rad: %.9g\n", (double)b->spin_phase_rad);

    fprintf(f, "    color: [%.9g, %.9g, %.9g]\n", (double)b->r, (double)b->g,
            (double)b->b);
  }

  fclose(f);
  return true;
}

