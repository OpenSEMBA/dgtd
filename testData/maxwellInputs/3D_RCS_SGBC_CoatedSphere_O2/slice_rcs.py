#!/usr/bin/env python3
"""Slice per-angle NTFF CSVs into the monostatic and bistatic comparison files."""

import argparse
import csv
import math
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

HEADER_RE = re.compile(r"theta_rad=([0-9eE.+-]+)\s+phi_rad=([0-9eE.+-]+)")
PI = math.pi


def load_angles(rcs_dir: Path):
    rows = []
    for path in sorted(rcs_dir.glob("rcsComplex_Th_*_Phi_*_dgtd.csv")):
        text = path.read_text().splitlines()
        match = None
        for line in text:
            match = HEADER_RE.search(line)
            if match:
                break
        if match is None:
            raise SystemExit(f"missing theta/phi header in {path}")
        theta = float(match.group(1))
        phi = float(match.group(2))
        data = []
        for rec in csv.DictReader(line for line in text if not line.startswith("#")):
            data.append({k: float(v) for k, v in rec.items()})
        rows.append((theta, phi, data))
    if not rows:
        raise SystemExit(f"no rcsComplex CSVs in {rcs_dir}")
    return rows


def nearest_freq(data, target_hz):
    return min(data, key=lambda rec: abs(rec["frequency_Hz"] - target_hz))


def write_csv(path: Path, fieldnames, records):
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        for rec in records:
            writer.writerow(rec)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--rcs-dir", type=Path, required=True)
    args = parser.parse_args()
    rcs_dir = args.rcs_dir
    summary = rcs_dir.parent / "summary"
    summary.mkdir(parents=True, exist_ok=True)

    angles = load_angles(rcs_dir)
    mono = min(angles, key=lambda item: abs(item[0] - PI) + abs(item[1]))
    eplane = [item for item in angles if abs(item[1]) < 1e-8]
    hplane = [item for item in angles if abs(item[1] - PI / 2) < 1e-8]
    if not eplane or not hplane:
        raise SystemExit("missing E-plane (phi=0) or H-plane (phi=pi/2) files")

    complex_fields = [
        "rcs_m2",
        "E_theta_re",
        "E_theta_im",
        "E_phi_re",
        "E_phi_im",
    ]
    mono_fields = ["frequency_Hz", *complex_fields]
    mono_rows = [{k: rec[k] for k in mono_fields} for rec in mono[2]]
    mono_rows.sort(key=lambda rec: rec["frequency_Hz"])
    write_csv(summary / "monostatic_rcs.csv", mono_fields, mono_rows)

    bistatic_fields = ["theta_deg", *complex_fields]
    cuts = []
    for freq_hz, tag in ((200e6, "200MHz"), (400e6, "400MHz")):
        for plane, name in ((eplane, "Eplane"), (hplane, "Hplane")):
            records = []
            for theta, _phi, data in sorted(plane, key=lambda item: item[0]):
                rec = nearest_freq(data, freq_hz)
                if abs(rec["frequency_Hz"] - freq_hz) > 1.0:
                    raise SystemExit(
                        f"no sample at {freq_hz} Hz (nearest {rec['frequency_Hz']})"
                    )
                row = {"theta_deg": theta * 180.0 / PI}
                row.update({k: rec[k] for k in complex_fields})
                records.append(row)
            out = summary / f"bistatic_{tag}_{name}.csv"
            write_csv(out, bistatic_fields, records)
            cuts.append((freq_hz, name, records))

    fig, ax = plt.subplots()
    ax.plot(
        [rec["frequency_Hz"] / 1e6 for rec in mono_rows],
        [rec["rcs_m2"] for rec in mono_rows],
    )
    ax.set_xlabel("Frequency (MHz)")
    ax.set_ylabel("RCS (m$^2$)")
    ax.set_title("Monostatic RCS, theta = 180 deg")
    fig.tight_layout()
    fig.savefig(summary / "monostatic_rcs.png", dpi=120)
    plt.close(fig)

    for freq_hz, tag in ((200e6, "200MHz"), (400e6, "400MHz")):
        fig, ax = plt.subplots()
        for name, style in (("Eplane", "-"), ("Hplane", "--")):
            records = next(item[2] for item in cuts if item[0] == freq_hz and item[1] == name)
            ax.plot(
                [rec["theta_deg"] for rec in records],
                [rec["rcs_m2"] for rec in records],
                style,
                label="E-plane" if name == "Eplane" else "H-plane",
            )
        ax.set_xlabel("Theta (deg)")
        ax.set_ylabel("RCS (m$^2$)")
        ax.set_title(f"Bistatic RCS at {tag}")
        ax.legend()
        fig.tight_layout()
        fig.savefig(summary / f"bistatic_{tag}.png", dpi=120)
        plt.close(fig)

    print(f"wrote {summary}")


if __name__ == "__main__":
    main()
