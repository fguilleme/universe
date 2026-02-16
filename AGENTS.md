# Repository Guidelines

## Project Structure & Module Organization

- `src/`: Core application code.
  - `main.c`: SDL2 window/input loop, camera, rendering, UI overlays.
  - `sim.c` / `sim.h`: N-body gravity simulation (currently 3D) + integration.
  - `glutil.c` / `glutil.h`: Shader compile/link helpers.
  - `textures.c` / `textures.h`: Planet texture array + optional HDR loader.
- `assets/`: Optional runtime assets (not committed). See `assets/README.md`.
- `third_party/`: Vendored single-header deps (e.g. `stb_image.h`, bitmap font).
- `scripts/`: Small helper scripts (e.g. checking/downloading assets).
- `build/`: Out-of-tree build directory (generated).

## Build, Test, and Development Commands

- `cmake -S . -B build`: Configure the project.
- `cmake --build build -j`: Build the `universe` binary.
- `./build/universe`: Run from repo root (recommended so `assets/` resolves).
- `UNIVERSE_ASSETS_DIR=/path/to/assets ./build/universe`: Override asset lookup.
- `bash scripts/check_assets.sh`: Show which expected textures/HDR are present.

## Coding Style & Naming Conventions

- Language: C11.
- Indentation: 2 spaces; keep functions small and readable.
- Naming: `snake_case` for functions and files; `CamelCase` for `typedef struct`.
- Avoid introducing new heavy dependencies; prefer small single-file libraries.

## Testing Guidelines

- No automated test suite yet.
- When changing physics/rendering, do a quick manual sanity pass:
  - Run solar system preset (`3`) and verify bodies move and UI overlays render.
  - Verify selection/picking and help overlays (`F1`, `F2`) still work.

## Commit & Pull Request Guidelines

- Keep commits focused and descriptive (e.g. “Render HDR background” vs “fix”).
- PRs should include:
  - What changed and why, with controls/screenshots for visual changes.
  - Notes on any new assets and their licenses/attribution (don’t commit them).

## Security & Configuration Tips

- Treat downloaded assets as untrusted input; the loaders should fail gracefully.
- Large files belong in `assets/` and should remain gitignored.

