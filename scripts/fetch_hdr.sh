#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
ASSETS_DIR="${UNIVERSE_ASSETS_DIR:-$ROOT_DIR/assets}"

mkdir -p "$ASSETS_DIR"

HDR_URL="https://www.spacespheremaps.com/wp-content/uploads/HDR_blue_nebulae-1.hdr"
OUT="$ASSETS_DIR/HDR_blue_nebulae-1.hdr"

echo "Downloading: $HDR_URL"
echo "To:         $OUT"

curl -L --fail --retry 3 -o "$OUT" "$HDR_URL"

echo "Done."

