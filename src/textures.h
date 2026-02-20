#pragma once

#include <GL/glew.h>

enum {
  TEX_GENERIC = 0,
  TEX_SUN = 1,
  TEX_MERCURY = 2,
  TEX_VENUS = 3,
  TEX_EARTH = 4,
  TEX_MARS = 5,
  TEX_JUPITER = 6,
  TEX_SATURN = 7,
  TEX_URANUS = 8,
  TEX_NEPTUNE = 9,
  TEX_MOON = 10,
  TEX_LAYER_COUNT = 11,
};

// Creates a GL_TEXTURE_2D_ARRAY with TEX_LAYER_COUNT layers.
// Loads files from `dir` when present; otherwise uses procedural fallbacks.
// Expected filenames (recommended from solarsystemscope packs):
//   2k_sun.jpg
//   2k_mercury.jpg
//   2k_venus_surface.jpg
//   2k_earth_daymap.jpg
//   2k_moon.jpg
//   2k_mars.jpg
//   2k_jupiter.jpg
//   2k_saturn.jpg
//   2k_uranus.jpg
//   2k_neptune.jpg
GLuint textures_create_planet_array(const char *dir);

// Tries to load an RGBA image into a GL_TEXTURE_2D from `dir`.
// Returns 0 if none of the candidate files can be loaded.
// Texture is created with linear filtering and clamp-to-edge.
GLuint textures_try_load_rgba2d(const char *dir, const char *label, const char *const *candidates, int *out_w, int *out_h);

// Loads a Radiance .hdr equirectangular environment map into a GL_TEXTURE_2D.
// Returns 0 if the file can't be loaded.
GLuint textures_load_hdr_equirect(const char *path);
