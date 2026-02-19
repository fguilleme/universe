# Universe Sandbox (minimal clone)

Standalone C / SDL2 / OpenGL gravity sandbox.

This is a small, self-contained N-body gravity simulator (Newtonian gravity) rendered with OpenGL point sprites.

It includes a 3D "solar system" preset (AU / years units) seeded from JPL Horizons barycentric state vectors (J2000 ecliptic frame), with:

- Non-circular orbits (eccentricity)
- Out-of-plane motion (inclination)
- Per-planet axial tilt + spin (texture rotation on a spherical impostor)

## Optional real textures

By default, the app uses procedural planet textures so it stays standalone.
If you want real textures, download them from https://www.solarsystemscope.com/textures/ and place them in `assets/`.

Note: asset paths are resolved relative to the current working directory; if you run the binary from another directory, set `UNIVERSE_ASSETS_DIR`.

Recommended filenames:

- `assets/2k_sun.jpg`
- `assets/2k_mercury.jpg`
- `assets/2k_venus_surface.jpg`
- `assets/2k_earth_daymap.jpg`
- `assets/2k_mars.jpg`
- `assets/2k_jupiter.jpg`
- `assets/2k_saturn.jpg`
- `assets/2k_uranus.jpg`
- `assets/2k_neptune.jpg`

The loader will resample any size to 512x512 internally.

Licensing: please ensure your use complies with the provider’s terms.

## Optional HDR background

Download an equirectangular `.hdr` and place it at:

- `assets/HDR_blue_nebulae-1.hdr`

Example: https://www.spacespheremaps.com/wp-content/uploads/HDR_blue_nebulae-1.hdr

If missing, the app renders a procedural starfield background.

### Helper scripts

- `scripts/check_assets.sh`: prints which expected files are present
- `scripts/fetch_hdr.sh`: downloads the example HDR into `assets/`

You can override the asset directory with:

```sh
UNIVERSE_ASSETS_DIR=/full/path/to/assets ./build/universe
```

## Build

Dependencies:

- SDL2
- GLEW
- OpenGL
- CMake >= 3.16

```sh
cmake -S . -B build
cmake --build build -j
./build/universe
```

Performance tips:

- CPU threads (gravity): `UNIVERSE_THREADS=8 ./build/universe`
- macOS Metal compute (gravity): `UNIVERSE_METAL=1 ./build/universe`

## Presets YAML

Preset definitions live in `data/presets.yaml`.

- Override path: `UNIVERSE_PRESETS_YAML=/path/to/presets.yaml ./build/universe`
- Each entry under `presets:` can be either:
  - `kind: bodies` with a `bodies:` list (like the solar system)
  - `kind: two_body` / `kind: disk_galaxy` with generator parameters

## Controls

- `F1`: toggle on-screen help
- `F2`: toggle body labels
- `LMB drag`: spawn a body; drag direction sets initial velocity (hold Shift for faster)
  - Hold `Ctrl` to spawn a circular orbit around selected (hold `Alt` to reverse)
- `RMB drag`: rotate camera
- `MMB drag` or `W/A/S/D`: pan camera
- `Mouse wheel`: zoom
- `Space`: pause/resume
- `R`: reset to solar system
- `1`: two-body orbit demo
- `2`: disk galaxy demo
- `3`: solar system demo
- `C`: clear all bodies
- `M`: toggle merge-on-collision
- `F`: follow selected body
- `T`: toggle trails
- `V`: toggle velocity vectors
- `[` / `]`: decrease/increase spawn mass
- `-` / `=`: slow down / speed up time
- `G` / `H`: decrease/increase gravity constant
- `Esc`: quit
