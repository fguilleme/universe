#include "presets.h"

#include <math.h>
#include <stdint.h>
#include <string.h>

#include "textures.h"

static uint32_t xorshift32(uint32_t *state) {
  uint32_t x = *state;
  x ^= x << 13;
  x ^= x >> 17;
  x ^= x << 5;
  *state = x;
  return x;
}

static double frand01(uint32_t *state) {
  return (double)(xorshift32(state) & 0xFFFFFF) / (double)0x1000000;
}

static double au_from_km(double km) { return km / 149597870.7; }

static double au_per_year_from_km_per_s(double km_s) {
  return km_s * 31557600.0 / 149597870.7;
}

static double visual_radius_from_physical(double r_au, double extra_scale) {
  // Exaggerate sizes while keeping relative ordering.
  // Use a low exponent to keep smaller bodies visible without making the Sun
  // unrealistically large compared to orbital distances.
  return pow(fmax(r_au, 1e-12), 0.30) * 0.60 * extra_scale;
}

void preset_seed_two_body_orbit(Sim *sim) {
  sim_reset(sim);
  sim->G = 1.0;
  sim->softening = 0.02;
  sim->density = 2.0;
  sim->merge_on_collision = true;

  const double sun_mass = 1.0e6;
  const double planet_mass = 10.0;
  const double r = 30.0;
  const double v = sqrt(sim->G * sun_mass / r);

  (void)sim_add_body(sim, (Body){.name = "Primary",
                                 .x = 0.0,
                                 .y = 0.0,
                                 .vx = 0.0,
                                 .vy = 0.0,
                                 .mass = sun_mass,
                                 .r = 1.0f,
                                 .g = 0.9f,
                                 .b = 0.6f});
  (void)sim_add_body(sim, (Body){.name = "Secondary",
                                 .x = r,
                                 .y = 0.0,
                                 .vx = 0.0,
                                 .vy = v,
                                 .mass = planet_mass,
                                 .r = 0.5f,
                                 .g = 0.8f,
                                 .b = 1.0f});
}

void preset_seed_disk_galaxy(Sim *sim, int n, double disk_radius,
                             double central_mass) {
  sim_reset(sim);
  sim->G = 1.0;
  sim->softening = 0.05;
  sim->density = 1.5;
  sim->merge_on_collision = true;

  (void)sim_add_body(sim, (Body){.name = "Center",
                                 .x = 0.0,
                                 .y = 0.0,
                                 .vx = 0.0,
                                 .vy = 0.0,
                                 .mass = central_mass,
                                 .r = 1.0f,
                                 .g = 0.85f,
                                 .b = 0.5f});

  uint32_t rng = 0x12345678u;
  for (int i = 0; i < n; i++) {
    const double u = frand01(&rng);
    const double v = frand01(&rng);
    const double r = disk_radius * sqrt(u);
    const double a = v * 2.0 * M_PI;
    const double x = r * cos(a);
    const double y = r * sin(a);

    // Circular-ish orbit around center with slight noise.
    const double speed = sqrt(sim->G * central_mass / fmax(r, 1e-3));
    const double tx = -sin(a);
    const double ty = cos(a);
    const double noise = (frand01(&rng) - 0.5) * 0.08;

    Body b = {0};
    b.x = x;
    b.y = y;
    b.vx = tx * speed * (1.0 + noise);
    b.vy = ty * speed * (1.0 + noise);
    b.mass = 1.0 + 15.0 * frand01(&rng);
    b.r = (float)(0.4 + 0.6 * frand01(&rng));
    b.g = (float)(0.4 + 0.6 * frand01(&rng));
    b.b = 1.0f;
    (void)sim_add_body(sim, b);
  }
}

void preset_seed_solar_system(Sim *sim) {
  // Accurate-ish 3D initial conditions from JPL Horizons (DE441),
  // barycentric ecliptic-of-J2000 frame at 2000-01-01 00:00 TDB.
  // Units in this simulation:
  // - distance: AU
  // - mass: solar mass
  // - time: years
  // In these units, G = 4*pi^2.
  sim_reset(sim);
  sim->G = 4.0 * M_PI * M_PI;
  sim->softening = 1e-6;
  sim->merge_on_collision = false;

  typedef struct {
    const char *name;
    double mass_solar;
    double radius_km;
    uint32_t layer;
    double tilt_deg;
    double rotation_period_days;
    double x_km, y_km, z_km;
    double vx_km_s, vy_km_s, vz_km_s;
  } BodyInit;

  const BodyInit bodies[] = {
      {"Sun", 1.0, 695700.0, TEX_SUN, 7.25, 25.38, -1.068108951496322E+06,
       -4.177210908491462E+05, 3.086887010002915E+04, 9.305302656256911E-03,
       -1.283177282717393E-02, -1.631700118015769E-04},
      {"Mercury", 1.660e-7, 2439.7, TEX_MERCURY, 0.034, 58.646,
       -2.212073002393702E+07, -6.682435921338345E+07, -3.461577076477692E+06,
       3.666229234452722E+01, -1.230266984222893E+01, -4.368336206255391E+00},
      {"Venus", 2.447e-6, 6051.8, TEX_VENUS, 177.36, -243.025,
       -1.075068040813442E+08, -1.053381041265540E+07,
       6.024630234834461E+06, 2.949586892816910E+00,
       -3.520915385418408E+01, -5.945882137071832E-01},
      {"Earth", 3.003e-6, 6371.0, TEX_EARTH, 23.439, 0.99726968,
       -2.521092392536846E+07, 1.449177150228508E+08,
       -2.394650324223201E+04, -2.983862964912880E+01,
       -5.218619426624726E+00, 1.193210383927141E-03},
      {"Mars", 3.227e-7, 3389.5, TEX_MARS, 25.19, 1.02595675,
       -2.155838897834890E+08, -3.539753056080124E+07,
       4.532312700768938E+06, 3.886195690590875E+00,
       -2.315847231020980E+01, -5.777378581341476E-01},
      {"Jupiter", 9.545e-4, 69911.0, TEX_JUPITER, 3.13, 0.41354,
       5.990175449652451E+08, 4.394132557364160E+08,
       -1.520962039552424E+07, -7.882322868602207E+00,
       1.134831390779115E+01, 1.246425392362529E-01},
      {"Saturn", 2.858e-4, 58232.0, TEX_SATURN, 26.73, 0.44401,
       9.587267807903839E+08, -9.825068848392826E+08,
       -1.695064654623866E+07, 5.908189281781105E+00,
       6.296823371892003E+00, -3.704231008630100E-01},
      {"Uranus", 4.366e-5, 25362.0, TEX_URANUS, 97.77, -0.71833,
       2.158708519561221E+09, 2.058822953381595E+09,
       -1.785026690391928E+07, -4.813504855508626E+00,
       4.665693498302723E+00, 8.548558001407977E-02},
      {"Neptune", 5.151e-5, 24622.0, TEX_NEPTUNE, 28.32, 0.67125,
       2.510985726269384E+09, -3.732999840511448E+09,
       1.873087787685826E+07, 4.445846366126391E+00,
       3.067772238252208E+00, -1.638828425643102E-01},
  };

  for (size_t i = 0; i < sizeof(bodies) / sizeof(bodies[0]); i++) {
    const BodyInit *in = &bodies[i];
    Body b = {0};
    strncpy(b.name, in->name, sizeof(b.name) - 1);
    b.name[sizeof(b.name) - 1] = 0;

    b.x = au_from_km(in->x_km);
    b.y = au_from_km(in->y_km);
    b.z = au_from_km(in->z_km);
    b.vx = au_per_year_from_km_per_s(in->vx_km_s);
    b.vy = au_per_year_from_km_per_s(in->vy_km_s);
    b.vz = au_per_year_from_km_per_s(in->vz_km_s);
    b.mass = in->mass_solar;

    const double r_au = au_from_km(in->radius_km);
    b.radius = r_au;
    b.render_radius = visual_radius_from_physical(r_au, 1.0);
    b.tex_layer = in->layer;
    b.tilt_rad = (float)(in->tilt_deg * (M_PI / 180.0));
    {
      const double pd = in->rotation_period_days;
      if (pd != 0.0) {
        const double sign = (pd >= 0.0) ? 1.0 : -1.0;
        const double pd_abs = fabs(pd);
        b.spin_rate_rad_per_time =
            (float)(sign * (2.0 * M_PI * 365.25) / pd_abs);
      }
    }
    b.spin_phase_rad = 0.0f;
    b.r = 1.0f;
    b.g = 1.0f;
    b.b = 1.0f;
    (void)sim_add_body(sim, b);
  }
}

