#!/usr/bin/env python3

from __future__ import annotations

import math
import re
import sys


def parse_bodies(path: str) -> dict[str, dict[str, float]]:
  txt = open(path, "r", encoding="utf-8", errors="replace").read().splitlines()
  bodies: dict[str, dict[str, float]] = {}
  in_bodies = False
  cur: str | None = None
  for line in txt:
    if line.strip() == "bodies:":
      in_bodies = True
      continue
    if not in_bodies:
      continue
    m = re.match(r"\s*-\s*name:\s*(.+)$", line)
    if m:
      cur = m.group(1).strip()
      bodies[cur] = {}
      continue
    if not cur:
      continue
    m = re.match(r"\s*(x_km|y_km|z_km|vx_km_s|vy_km_s|vz_km_s):\s*([^#]+)", line)
    if m:
      k = m.group(1)
      v = float(m.group(2).strip())
      bodies[cur][k] = v
  return bodies


def main() -> int:
  path = sys.argv[1] if len(sys.argv) > 1 else "data/presets.yaml"
  bodies = parse_bodies(path)
  for name in ("Earth", "Moon"):
    if name not in bodies:
      print(f"Missing body '{name}' in {path}")
      return 2
    for k in ("x_km", "y_km", "z_km"):
      if k not in bodies[name]:
        print(f"Missing {name}.{k} in {path}")
        return 2

  ex, ey, ez = bodies["Earth"]["x_km"], bodies["Earth"]["y_km"], bodies["Earth"]["z_km"]
  mx, my, mz = bodies["Moon"]["x_km"], bodies["Moon"]["y_km"], bodies["Moon"]["z_km"]
  dx, dy, dz = mx - ex, my - ey, mz - ez
  d = math.sqrt(dx * dx + dy * dy + dz * dz)
  print(f"Earth-Moon distance: {d:,.3f} km")
  if d > 1_000_000:
    print("WARNING: distance is very large; likely wrong Horizons center/frame/epoch")
  return 0


if __name__ == "__main__":
  raise SystemExit(main())

