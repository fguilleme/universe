#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
ASSETS_DIR="${UNIVERSE_ASSETS_DIR:-$ROOT_DIR/assets}"

echo "Assets dir: $ASSETS_DIR"

req=(
  2k_sun.jpg
  2k_mercury.jpg
  2k_venus_surface.jpg
  2k_earth_daymap.jpg
  2k_mars.jpg
  2k_jupiter.jpg
  2k_saturn.jpg
  2k_uranus.jpg
  2k_neptune.jpg
  HDR_blue_nebulae-1.hdr
)

missing=0
for f in "${req[@]}"; do
  if [[ -f "$ASSETS_DIR/$f" ]]; then
    echo "OK:      $f"
  else
    echo "MISSING: $f"
    missing=1
  fi
done

exit $missing

