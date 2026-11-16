#!/usr/bin/env bash
# Rebuild the linear G1 mesh. Elementary tags are set equal to physical tags
# so fixGmshMesh / MFEM attributes match the JSON material tags (as on G2).
set -euo pipefail
cd "$(dirname "$0")"
gmsh -2 2D_Dipole_PML_G1.geo -format msh2 -o 2D_Dipole_PML_G1.msh
python3 - <<'PY'
from pathlib import Path
p = Path("2D_Dipole_PML_G1.msh")
t = p.read_text()
pre, rest = t.split("$Elements\n", 1)
body, post = rest.split("$EndElements", 1)
lines = body.strip().splitlines()
out = [lines[0]]
for line in lines[1:]:
    parts = line.split()
    parts[4] = parts[3]
    out.append(" ".join(parts))
p.write_text(pre + "$Elements\n" + "\n".join(out) + "\n$EndElements" + post)
print("wrote", p, "with elementary := physical")
PY
