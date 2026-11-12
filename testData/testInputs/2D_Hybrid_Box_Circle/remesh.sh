#!/usr/bin/env bash
# Rebuild hybrid mesh and set elementary := physical (fixGmshMesh / MFEM attrs).
# Also flips any remaining inverted triangles (negative signed area).
set -euo pipefail
cd "$(dirname "$0")"
# Elementary surfaces are tagged 1..14 outer->inner; gmsh -2 meshes in tag order.
gmsh -2 2D_Hybrid_Box_Circle.geo -format msh2 -o 2D_Hybrid_Box_Circle.msh
python3 - <<'PY'
from pathlib import Path
p = Path("2D_Hybrid_Box_Circle.msh")
t = p.read_text()
pre, rest = t.split("$Elements\n", 1)
body, post = rest.split("$EndElements", 1)
lines = body.strip().splitlines()

# Parse nodes for orientation check
nodes = {}
nb = t.split("$Nodes\n")[1].split("$EndNodes")[0].strip().splitlines()
n = int(nb[0])
for L in nb[1 : n + 1]:
    i, x, y, z = L.split()[:4]
    nodes[int(i)] = (float(x), float(y))

out = [lines[0]]
pos = neg = zero = flipped = 0
for line in lines[1:]:
    parts = line.split()
    et = int(parts[1])
    ntags = int(parts[2])
    if ntags >= 2:
        parts[4] = parts[3]  # elementary := physical
    if et == 2:  # triangle
        ids = [int(x) for x in parts[3 + ntags : 6 + ntags]]
        a, b, c = nodes[ids[0]], nodes[ids[1]], nodes[ids[2]]
        area = 0.5 * ((b[0] - a[0]) * (c[1] - a[1]) - (c[0] - a[0]) * (b[1] - a[1]))
        if area < -1e-18:
            # swap last two nodes → flip orientation
            parts[3 + ntags + 1], parts[3 + ntags + 2] = (
                parts[3 + ntags + 2],
                parts[3 + ntags + 1],
            )
            flipped += 1
            neg += 1
        elif area > 1e-18:
            pos += 1
        else:
            zero += 1
    out.append(" ".join(parts))
p.write_text(pre + "$Elements\n" + "\n".join(out) + "\n$EndElements" + post)
print(
    f"wrote {p}: elementary:=physical; tri orient before flip "
    f"pos={pos} neg={neg} zero={zero}; flipped={flipped}"
)
if zero:
    raise SystemExit(f"ERROR: {zero} zero-area triangles")
PY
