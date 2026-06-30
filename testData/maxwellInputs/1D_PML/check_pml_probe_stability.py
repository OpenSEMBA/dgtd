#!/usr/bin/env python3
"""PASS if max |Ey| and max |Hz| <= threshold at each uniform PML probe."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

THRESHOLD = 3.0
COLUMNS = ("time", "Ex", "Ey", "Ez", "Hx", "Hy", "Hz")


def parse_probe_dat(path: Path) -> list[dict[str, float]]:
    rows: list[dict[str, float]] = []
    in_data = False
    with path.open() as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith("Time (s)"):
                in_data = True
                continue
            if not in_data:
                continue
            parts = line.split()
            if len(parts) < 7:
                continue
            vals = dict(zip(COLUMNS, (float(x) for x in parts[:7])))
            rows.append(vals)
    return rows


def probe_stats(path: Path) -> tuple[float, float]:
    rows = parse_probe_dat(path)
    if not rows:
        raise ValueError(f"no data rows in {path}")
    max_ey = max(abs(r["Ey"]) for r in rows)
    max_hz = max(abs(r["Hz"]) for r in rows)
    return max_ey, max_hz


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "export_dir",
        type=Path,
        help="Exports/<mode>/<case>/PointProbes (or parent containing PointProbes)",
    )
    parser.add_argument(
        "--threshold",
        type=float,
        default=THRESHOLD,
        help=f"max |Ey|, |Hz| limit (default {THRESHOLD})",
    )
    parser.add_argument(
        "--probes",
        type=int,
        nargs="*",
        default=[0, 1, 2],
        help="PointProbe IDs (default 0 1 2)",
    )
    args = parser.parse_args()

    base = args.export_dir
    if (base / "PointProbes").is_dir():
        base = base / "PointProbes"

    ok = True
    print(f"threshold={args.threshold}")
    for pid in args.probes:
        path = base / f"PointProbe{pid}.dat"
        if not path.is_file():
            print(f"FAIL probe {pid}: missing {path}")
            ok = False
            continue
        max_ey, max_hz = probe_stats(path)
        passed = max_ey <= args.threshold and max_hz <= args.threshold
        status = "PASS" if passed else "FAIL"
        print(
            f"{status} PointProbe{pid}: max|Ey|={max_ey:.6e} max|Hz|={max_hz:.6e}"
        )
        ok = ok and passed

    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
