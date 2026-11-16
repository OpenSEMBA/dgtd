// Hybrid 2D mesh: unstructured core (circle -> TFSF), structured quads outside.
//
//   circle r=0.5
//   TFSF        |x|,|y|=1     Physical Curve 9
//   NTFF        |x|,|y|=3     Physical Curve 11
//   PML         |x|,|y|=4->5   Physical Surfaces 1–4
//   Outer SMA   |x|,|y|=5     Physical Curves 4,6,13,15
//
// Face convention: surfaces are *meshed outer -> inner* (PML, then 3->4, then
// 1->3, then TFSF core, then disk) so shared-face Elem1/Elem2 orientation is
// established from the truncated exterior inward.
//
// Topology is still built inside-out (shared edges). Only the Mesh Surface
// calls run outside-in.
//
// Tangential node count is identical on every square (Transfinite requirement).
// Rebuild: ./remesh.sh

Mesh.Algorithm = 6;
Mesh.ElementOrder = 1;
Mesh.SecondOrderLinear = 1;

r = 0.5; a1 = 1; a3 = 3; a4 = 4; a5 = 5;
n_arc = 16;
n_tang = 80;
n_r13 = 40; n_r34 = 20; n_r45 = 24;

// ---------------------------------------------------------------------------
// Geometry (points / curves). Rings created inside-out for shared edges.
// ---------------------------------------------------------------------------
Point(1) = {0, 0, 0};
Point(2) = { r,  0, 0};
Point(3) = { 0,  r, 0};
Point(4) = {-r,  0, 0};
Point(5) = { 0, -r, 0};
Circle(1) = {2, 1, 3};
Circle(2) = {3, 1, 4};
Circle(3) = {4, 1, 5};
Circle(4) = {5, 1, 2};
Curve Loop(1) = {1, 2, 3, 4};
Transfinite Curve{1:4} = n_arc;

Macro Corners
  Point(p+0) = {A, -A, 0};
  Point(p+1) = {A,  A, 0};
  Point(p+2) = {-A, A, 0};
  Point(p+3) = {-A,-A, 0};
Return
A=a1; p=10; Call Corners;
A=a3; p=20; Call Corners;
A=a4; p=30; Call Corners;
A=a5; p=40; Call Corners;

// TFSF square edges
Line(10) = {10, 11};
Line(11) = {11, 12};
Line(12) = {12, 13};
Line(13) = {13, 10};
Transfinite Curve{10:13} = n_tang;
Curve Loop(2) = {10, 11, 12, 13};

// Ring line macro: outer sides + radials (no surfaces yet)
Macro RingLines
  Line(eO+0) = {pO+0, pO+1};
  Line(eO+1) = {pO+1, pO+2};
  Line(eO+2) = {pO+2, pO+3};
  Line(eO+3) = {pO+3, pO+0};
  Line(eR+0) = {pI+0, pO+0};
  Line(eR+1) = {pI+1, pO+1};
  Line(eR+2) = {pI+2, pO+2};
  Line(eR+3) = {pI+3, pO+3};
  Transfinite Curve{eO:eO+3} = n_tang;
  Transfinite Curve{eR:eR+3} = nR;
  // CCW loops: out on SE radial, along outer CCW, back on NE radial, reverse inner.
  Curve Loop(c+0) = {eR+0, eO+0, -(eR+1), -(eI+0)};
  Curve Loop(c+1) = {eR+1, eO+1, -(eR+2), -(eI+1)};
  Curve Loop(c+2) = {eR+2, eO+2, -(eR+3), -(eI+2)};
  Curve Loop(c+3) = {eR+3, eO+3, -(eR+0), -(eI+3)};
Return

// 1->3 vacuum lines
pI=10; pO=20; eI=10; eO=20; eR=30; c=20; nR=n_r13; Call RingLines;
// 3->4 vacuum lines
pI=20; pO=30; eI=20; eO=40; eR=50; c=40; nR=n_r34; Call RingLines;
// 4->5 PML lines
pI=30; pO=40; eI=40; eO=60; eR=70; c=60; nR=n_r45; Call RingLines;

MeshSize{ PointsOf{ Curve{1:4, 10:13}; } } = 2*a1/n_tang;

// ---------------------------------------------------------------------------
// Surfaces numbered and defined OUTER -> INNER (low tags = exterior)
//   1–4   PML +X,+Y,-X,-Y   (loops c=60)
//   5–8   vacuum ring 3->4   (loops c=40)
//   9–12  vacuum ring 1->3   (loops c=20)
//   13    TFSF square \ disk
//   14    disk
// ---------------------------------------------------------------------------
Macro QuadSurf
  // CCW corner order (p0→p3→p2→p1) so Transfinite tris have positive area.
  Plane Surface(s) = {cl};
  Transfinite Surface{s} = {p0, p3, p2, p1};
Return

// PML 4->5: loops 60..63, points 30..33 -> 40..43
s=1; cl=60; p0=30; p1=31; p2=41; p3=40; Call QuadSurf; // +X
s=2; cl=61; p0=31; p1=32; p2=42; p3=41; Call QuadSurf; // +Y
s=3; cl=62; p0=32; p1=33; p2=43; p3=42; Call QuadSurf; // -X
s=4; cl=63; p0=33; p1=30; p2=40; p3=43; Call QuadSurf; // -Y

// Vacuum 3->4
s=5; cl=40; p0=20; p1=21; p2=31; p3=30; Call QuadSurf;
s=6; cl=41; p0=21; p1=22; p2=32; p3=31; Call QuadSurf;
s=7; cl=42; p0=22; p1=23; p2=33; p3=32; Call QuadSurf;
s=8; cl=43; p0=23; p1=20; p2=30; p3=33; Call QuadSurf;

// Vacuum 1->3
s=9;  cl=20; p0=10; p1=11; p2=21; p3=20; Call QuadSurf;
s=10; cl=21; p0=11; p1=12; p2=22; p3=21; Call QuadSurf;
s=11; cl=22; p0=12; p1=13; p2=23; p3=22; Call QuadSurf;
s=12; cl=23; p0=13; p1=10; p2=20; p3=23; Call QuadSurf;

Plane Surface(13) = {2, 1}; // TFSF square \ disk
Plane Surface(14) = {1};    // disk

// Physical groups (material / BC tags)
Physical Surface(4) = {1};              // +X PML
Physical Surface(1) = {2};              // +Y PML
Physical Surface(3) = {3};              // -X PML
Physical Surface(2) = {4};              // -Y PML
Physical Surface(9) = {5, 6, 7, 8, 9, 10, 11, 12, 13};
Physical Surface(10) = {14};            // disk

Physical Curve(1) = {1, 2, 3, 4};       // circle
Physical Curve(9) = {10, 11, 12, 13};   // TFSF
Physical Curve(11) = {20, 21, 22, 23};  // NTFF
Physical Curve(15) = {60};              // SMA +X
Physical Curve(4) = {61};               // SMA +Y
Physical Curve(13) = {62};              // SMA -X
Physical Curve(6) = {63};               // SMA -Y

// Mesh order: elementary surfaces are tagged 1..14 outer->inner.
// Gmsh's default 2D pass meshes by increasing surface tag, so faces are
// established from the PML inward to the disk (Elem1/Elem2 / TFSF faces).
