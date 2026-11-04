#!/usr/bin/env python3
"""Compare a 2D coated-cylinder SGBC run with Figure 7 of the SGBC paper.

The reference columns were read off the published plot of normalized
monostatic scattering width, sigma_2D / (2a), in dB versus ka. The solid
reference and SGBC G2 lie on top of each other in that figure, so the
test matches the DGTD backscatter to the digitized solid curve over the
energetic band ka = 0.5 to 4.09. Digitization is about 0.4 dB, wider than
the 0.03 dB G2-to-reference difference reported in the text.
"""

import argparse
import csv
import math
import sys
from pathlib import Path

C0 = 299792458.0
RADIUS_M = 1.0
KA_MIN = 0.5
KA_MAX = 4.09
MEAN_DB = 1.0
MAX_DB = 2.5

HERE = Path(__file__).resolve().parent
REFERENCE = HERE / "fig7_digitized.csv"


def load_reference(path):
    rows = []
    with path.open() as handle:
        for rec in csv.DictReader(handle):
            rows.append((float(rec["ka"]), float(rec["sim_ref_db"])))
    if len(rows) < 10:
        raise SystemExit(f"reference {path} has too few samples")
    return rows


def load_dgtd(path):
    samples = []
    with path.open() as handle:
        header = handle.readline()
        if "rcs" not in header:
            raise SystemExit(f"unexpected RCS header in {path}: {header.strip()}")
        for line in handle:
            parts = line.split()
            if len(parts) < 4:
                continue
            freq = float(parts[2])
            rcs = float(parts[3])
            ka = 2.0 * math.pi * freq / C0 * RADIUS_M
            if rcs <= 0.0:
                continue
            db = 10.0 * math.log10(rcs / (2.0 * RADIUS_M))
            samples.append((ka, db))
    if not samples:
        raise SystemExit(f"no RCS samples in {path}")
    return samples


def interp(reference, ka):
    if ka <= reference[0][0]:
        return reference[0][1]
    if ka >= reference[-1][0]:
        return reference[-1][1]
    for i in range(1, len(reference)):
        k1, v1 = reference[i]
        k0, v0 = reference[i - 1]
        if ka <= k1:
            span = k1 - k0
            if span == 0.0:
                return v1
            w = (ka - k0) / span
            return v0 + w * (v1 - v0)
    return reference[-1][1]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("rcs_dat", type=Path)
    args = parser.parse_args()
    reference = load_reference(REFERENCE)
    samples = load_dgtd(args.rcs_dat)
    errors = []
    for ka, db in samples:
        if not (KA_MIN <= ka <= KA_MAX):
            continue
        errors.append((abs(db - interp(reference, ka)), ka, db, interp(reference, ka)))
    if len(errors) < 5:
        raise SystemExit(
            f"only {len(errors)} samples inside ka [{KA_MIN}, {KA_MAX}]"
        )
    mean = sum(item[0] for item in errors) / len(errors)
    worst = max(errors)
    print(
        f"ka {KA_MIN}-{KA_MAX}: n={len(errors)} "
        f"mean |err|={mean:.3f} dB  max={worst[0]:.3f} dB at ka={worst[1]:.3f} "
        f"(DGTD {worst[2]:.2f}, figure {worst[3]:.2f})"
    )
    if mean > MEAN_DB or worst[0] > MAX_DB:
        print(
            f"FAIL mean limit {MEAN_DB} dB, point limit {MAX_DB} dB",
            file=sys.stderr,
        )
        return 1
    print("PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
