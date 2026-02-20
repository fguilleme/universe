#!/usr/bin/env python3

"""Convert a JPL Horizons 'Vectors' output to a YAML body entry.

Usage:
  python scripts/horizons_vectors_to_yaml.py --name Moon --tex_layer TEX_MOON \
    --mass_solar 3.694e-8 --radius_km 1737.4 --tilt_deg 6.68 \
    --rotation_period_days 27.321661 < horizons.txt

Expected Horizons settings (recommended):
  - Ephemeris Type: VECTORS
  - Center: Solar System Barycenter (500@0)
  - Reference frame: Ecliptic of J2000
  - Units: KM-S
  - Output: full precision

The script looks for a single state vector block and extracts X,Y,Z,VX,VY,VZ.
"""

from __future__ import annotations

import argparse
import re
import subprocess
import sys


_NUM = r"[-+]?\d+(?:\.\d+)?(?:[Ee][-+]?\d+)?"


def _slice_soe_block(text: str) -> str:
  # Horizons text output usually wraps the data between $$SOE / $$EOE.
  i = text.find("$$SOE")
  if i == -1:
    return text
  j = text.find("$$EOE", i)
  if j == -1:
    return text[i:]
  return text[i:j]


def _split_columns(line: str) -> list[str]:
  # Split either CSV header or space-aligned header.
  if "," in line:
    cols = [c.strip() for c in line.split(",")]
  else:
    cols = [c.strip() for c in re.split(r"\s+", line.strip())]
  cols = [c for c in cols if c]
  return cols


def _find_table_state(text: str) -> dict[str, float] | None:
  lines = [ln.strip() for ln in text.splitlines() if ln.strip()]
  # Look for a header that mentions our desired components.
  want_syms = {"X", "Y", "Z", "VX", "VY", "VZ"}
  for i, line in enumerate(lines):
    up = line.upper()
    if not all(sym in up for sym in want_syms):
      continue

    cols = [c.upper() for c in _split_columns(line)]
    idx = {name: j for j, name in enumerate(cols)}
    if not all(sym in idx for sym in want_syms):
      continue

    # Next line with enough numbers is assumed to be the data row.
    for j in range(i + 1, min(i + 10, len(lines))):
      nums = re.findall(_NUM, lines[j])
      if len(nums) < len(cols):
        continue
      vals = [float(x) for x in nums[: len(cols)]]
      return {
          "x_km": vals[idx["X"]],
          "y_km": vals[idx["Y"]],
          "z_km": vals[idx["Z"]],
          "vx_km_s": vals[idx["VX"]],
          "vy_km_s": vals[idx["VY"]],
          "vz_km_s": vals[idx["VZ"]],
      }
  return None


def _find_state(text: str) -> dict[str, float]:
  text = _slice_soe_block(text)

  # Robust extraction: find X/Y/Z/VX/VY/VZ assignments anywhere in the block.
  # Horizons formats vary, e.g.:
  #   X = ... Y = ... Z = ...
  #   VX= ... VY= ... VZ= ...
  # or with units:
  #   X (km) = ...
  want = {
      "x_km": "X",
      "y_km": "Y",
      "z_km": "Z",
      "vx_km_s": "VX",
      "vy_km_s": "VY",
      "vz_km_s": "VZ",
  }
  out: dict[str, float] = {}
  for key, sym in want.items():
    m = re.search(
        rf"\b{sym}\b[^=]*=\s*({_NUM})",
        text,
        flags=re.IGNORECASE,
    )
    if not m:
      table = _find_table_state(text)
      if table is not None:
        return table
      raise ValueError(
          f"Could not find {sym}=... in input. "
          "Clipboard likely contains a table output; ensure you copied the VECTORS block (between $$SOE/$$EOE), "
          "or copy the header+data row that includes X,Y,Z,VX,VY,VZ."
      )
    out[key] = float(m.group(1))
  return out


def main() -> int:
  ap = argparse.ArgumentParser()
  ap.add_argument("--name", required=True)
  ap.add_argument("--tex_layer", required=True)
  ap.add_argument("--mass_solar", required=True, type=float)
  ap.add_argument("--radius_km", required=True, type=float)
  ap.add_argument("--tilt_deg", required=True, type=float)
  ap.add_argument("--rotation_period_days", required=True, type=float)
  ap.add_argument(
      "--clipboard",
      action="store_true",
      help="Read input from macOS clipboard via pbpaste (avoids waiting for EOF)",
  )
  ap.add_argument(
      "input",
      nargs="?",
      default="-",
      help="Input file path (default: stdin)",
  )
  args = ap.parse_args()

  if args.clipboard:
    try:
      text = subprocess.check_output(["pbpaste"], text=True)
    except Exception as e:
      raise SystemExit(f"Failed to read clipboard via pbpaste: {e}")
  else:
    if args.input != "-":
      with open(args.input, "r", encoding="utf-8", errors="replace") as f:
        text = f.read()
    else:
      if sys.stdin.isatty():
        print(
            "Paste the Horizons VECTORS output, then press Ctrl-D to finish...",
            file=sys.stderr,
        )
      text = sys.stdin.read()
  state = _find_state(text)

  # Emit YAML snippet compatible with data/presets.yaml.
  print(f"      - name: {args.name}")
  print(f"        mass_solar: {args.mass_solar:.12g}")
  print(f"        radius_km: {args.radius_km:.6g}")
  print(f"        tex_layer: {args.tex_layer}")
  print(f"        tilt_deg: {args.tilt_deg:.6g}")
  print(f"        rotation_period_days: {args.rotation_period_days:.9g}")
  print(f"        x_km: {state['x_km']:.15E}")
  print(f"        y_km: {state['y_km']:.15E}")
  print(f"        z_km: {state['z_km']:.15E}")
  print(f"        vx_km_s: {state['vx_km_s']:.15E}")
  print(f"        vy_km_s: {state['vy_km_s']:.15E}")
  print(f"        vz_km_s: {state['vz_km_s']:.15E}")
  return 0


if __name__ == "__main__":
  raise SystemExit(main())
