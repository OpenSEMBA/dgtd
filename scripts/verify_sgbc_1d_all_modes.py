#!/usr/bin/env python3
"""Automated, non-interactive SGBC 1D verification across all supported run modes.

Runs the automated 1D SGBC scikit-rf GoogleTests in `cases_tests`:
  - CasesTest.1D_SGBC_Multilayer_PEC_Backed_RAM_vs_Skrf
  - CasesTest.1D_SGBC_Resonant_Passthrough_Slab_vs_Skrf
  - CasesTest.1D_SGBC_Through_Passthrough_vs_Skrf
  - CasesTest.1D_SGBC_Through_PEC_Backed_vs_Skrf
  - CasesTest.1D_SGBC_Multiple_Passthrough_vs_Skrf
  - CasesTest.1D_SGBC_Multiple_PEC_Backed_vs_Skrf

across all four supported execution modes:
  1. Serial CPU (1 rank, -d cpu): all 6 tests (including the long 24k/120k-step paper slabs)
  2. MPI-2 CPU  (2 ranks, -d cpu): multi-rank tests (`1D_SGBC_Through*` and `1D_SGBC_Multiple*`)
  3. MPI-4 CPU  (4 ranks, -d cpu): multi-rank tests (`1D_SGBC_Through*` and `1D_SGBC_Multiple*`)
  4. Single-rank CUDA GPU (1 rank, -d cuda): GPU tests (`1D_SGBC_Through*` and `1D_SGBC_Multiple*`)

Exits 0 if all modes pass, or 1 if any mode fails.
"""

from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
CANDIDATE_BINS = (
    ROOT / "build/gnu-release-cuda-sm120/bin/cases_tests",
    ROOT / "build/gnu-release-cuda/bin/cases_tests",
    ROOT / "build/gnu-release-mpi/bin/cases_tests",
)


def find_cases_tests_bin() -> Path:
    for candidate in CANDIDATE_BINS:
        if candidate.is_file():
            return candidate
    raise FileNotFoundError(
        "Could not find cases_tests binary in build/gnu-release-{cuda-sm120,cuda,mpi}/bin/"
    )


def main() -> int:
    cases_bin = find_cases_tests_bin()
    has_cuda = "cuda" in str(cases_bin)
    all_filter = "--gtest_filter=*1D_SGBC*"
    parallel_filter = "--gtest_filter=*1D_SGBC_Through*:*1D_SGBC_Multiple*"

    modes = [
        ("Serial CPU (1 rank)", [str(cases_bin), "-d", "cpu", all_filter]),
        ("MPI-2 CPU (2 ranks)", ["mpirun", "-np", "2", str(cases_bin), "-d", "cpu", parallel_filter]),
        ("MPI-4 CPU (4 ranks)", ["mpirun", "-np", "4", str(cases_bin), "-d", "cpu", parallel_filter]),
    ]
    if has_cuda:
        modes.append(
            ("CUDA GPU (1 rank)", [str(cases_bin), "-d", "cuda", parallel_filter])
        )

    failures = []
    for label, cmd in modes:
        print(f"=== Running {label}: {' '.join(cmd)} ===", flush=True)
        proc = subprocess.run(cmd, cwd=ROOT)
        if proc.returncode == 0:
            print(f"[PASS] {label}\n", flush=True)
        else:
            print(f"[FAIL] {label} (exit code {proc.returncode})\n", flush=True)
            failures.append(label)

    if failures:
        print(f"SGBC_VERIFICATION: FAIL ({len(failures)} mode(s) failed: {', '.join(failures)})")
        return 1

    print(f"SGBC_VERIFICATION: PASS (all {len(modes)} modes passed)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
