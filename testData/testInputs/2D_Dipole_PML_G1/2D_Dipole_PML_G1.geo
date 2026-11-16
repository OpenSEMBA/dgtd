// 2D Dipole + box PML — geometry order 1 (linear triangles everywhere).
// Vacuum |x|,|y|<=3 (tags 9,10), PML out to |x|,|y|=4 (tags 1–8).
// Dipole circle r=0.5 (curve tag 1, disk tag 10). Outer SMA tags match G2 JSON.
// Rebuild with ./remesh.sh (sets elementary := physical for fixGmshMesh).
SetFactory("OpenCASCADE");

r = 0.5;
a = 3.0;
b = 4.0;

Disk(1) = {0, 0, 0, r};
Rectangle(2) = {-a, -a, 0, 2*a, 2*a};           // vacuum square
Rectangle(3) = {-a, a, 0, 2*a, b-a};             // +Y slab
Rectangle(4) = {-a, -b, 0, 2*a, b-a};            // -Y slab
Rectangle(5) = {-b, -a, 0, b-a, 2*a};            // -X slab
Rectangle(6) = {a, -a, 0, b-a, 2*a};             // +X slab
Rectangle(7) = {a, a, 0, b-a, b-a};              // +X+Y corner
Rectangle(8) = {a, -b, 0, b-a, b-a};             // +X-Y corner
Rectangle(9) = {-b, -b, 0, b-a, b-a};            // -X-Y corner
Rectangle(10) = {-b, a, 0, b-a, b-a};            // -X+Y corner

BooleanFragments{ Surface{1:10}; Delete; }{}

// Map surfaces by centre (fragment renumbers).
s_list() = Surface{:};
disk = 0; vac = {};
pml1 = {}; pml2 = {}; pml3 = {}; pml4 = {};
pml5 = {}; pml6 = {}; pml7 = {}; pml8 = {};

For i In {0:#s_list()-1}
  s = s_list(i);
  bb() = BoundingBox Surface{s};
  cx = 0.5*(bb(0)+bb(3));
  cy = 0.5*(bb(1)+bb(4));
  dx = bb(3)-bb(0); dy = bb(4)-bb(1);
  If (dx < 1.2 && dy < 1.2 && Sqrt(cx*cx + cy*cy) < 0.3)
    disk = s;
  ElseIf (Abs(cx) < a-1e-3 && Abs(cy) < a-1e-3)
    vac() += {s};
  ElseIf (Abs(cx) < a-1e-3 && cy > 0)
    pml1() += {s};
  ElseIf (Abs(cx) < a-1e-3 && cy < 0)
    pml2() += {s};
  ElseIf (Abs(cy) < a-1e-3 && cx < 0)
    pml3() += {s};
  ElseIf (Abs(cy) < a-1e-3 && cx > 0)
    pml4() += {s};
  ElseIf (cx > 0 && cy > 0)
    pml5() += {s};
  ElseIf (cx > 0 && cy < 0)
    pml6() += {s};
  ElseIf (cx < 0 && cy < 0)
    pml7() += {s};
  ElseIf (cx < 0 && cy > 0)
    pml8() += {s};
  EndIf
EndFor

If (disk == 0)
  Error("Failed to identify dipole disk surface");
EndIf

Physical Surface(1) = pml1();
Physical Surface(2) = pml2();
Physical Surface(3) = pml3();
Physical Surface(4) = pml4();
Physical Surface(5) = pml5();
Physical Surface(6) = pml6();
Physical Surface(7) = pml7();
Physical Surface(8) = pml8();
Physical Surface(9) = vac();
Physical Surface(10) = {disk};

// Boundary curves
c_list() = Unique(Abs(Boundary{ Surface{s_list()}; }));
dipole_c = {};
sma4 = {}; sma6 = {}; sma13 = {}; sma15 = {};
sma18 = {}; sma19 = {}; sma20 = {}; sma21 = {};
sma22 = {}; sma23 = {}; sma24 = {}; sma25 = {};

For i In {0:#c_list()-1}
  c = c_list(i);
  bb() = BoundingBox Curve{c};
  cx = 0.5*(bb(0)+bb(3));
  cy = 0.5*(bb(1)+bb(4));
  // Full circle (or arc) around the origin: bbox ≈ [-r,r]^2
  If (Abs(cx) < 1e-3 && Abs(cy) < 1e-3 &&
      Abs(bb(0) + r) < 0.05 && Abs(bb(3) - r) < 0.05 &&
      Abs(bb(1) + r) < 0.05 && Abs(bb(4) - r) < 0.05)
    dipole_c() += {c};
  ElseIf (Abs(bb(1)-b) < 1e-6 && Abs(bb(4)-b) < 1e-6)
    If (bb(0) >= a-1e-6)
      sma19() += {c};
    ElseIf (bb(3) <= -a+1e-6)
      sma24() += {c};
    Else
      sma4() += {c};
    EndIf
  ElseIf (Abs(bb(1)+b) < 1e-6 && Abs(bb(4)+b) < 1e-6)
    If (bb(0) >= a-1e-6)
      sma20() += {c};
    ElseIf (bb(3) <= -a+1e-6)
      sma22() += {c};
    Else
      sma6() += {c};
    EndIf
  ElseIf (Abs(bb(0)-b) < 1e-6 && Abs(bb(3)-b) < 1e-6)
    If (bb(1) >= a-1e-6)
      sma18() += {c};
    ElseIf (bb(4) <= -a+1e-6)
      sma21() += {c};
    Else
      sma15() += {c};
    EndIf
  ElseIf (Abs(bb(0)+b) < 1e-6 && Abs(bb(3)+b) < 1e-6)
    If (bb(1) >= a-1e-6)
      sma25() += {c};
    ElseIf (bb(4) <= -a+1e-6)
      sma23() += {c};
    Else
      sma13() += {c};
    EndIf
  EndIf
EndFor

Physical Curve(1) = dipole_c();
Physical Curve(4) = sma4();
Physical Curve(6) = sma6();
Physical Curve(13) = sma13();
Physical Curve(15) = sma15();
Physical Curve(18) = sma18();
Physical Curve(19) = sma19();
Physical Curve(20) = sma20();
Physical Curve(21) = sma21();
Physical Curve(22) = sma22();
Physical Curve(23) = sma23();
Physical Curve(24) = sma24();
Physical Curve(25) = sma25();

lc_vac = 0.25;
lc_pml = 0.12;
lc_circ = 0.06;
n_circ = 48;   // finer than G2's 11 quadratic edges

MeshSize{ PointsOf{ Surface{disk}; } } = lc_circ;
MeshSize{ PointsOf{ Curve{dipole_c()}; } } = lc_circ;
MeshSize{ PointsOf{ Surface{vac()}; } } = lc_vac;
MeshSize{ PointsOf{ Surface{pml1(), pml2(), pml3(), pml4(),
                            pml5(), pml6(), pml7(), pml8()}; } } = lc_pml;
Transfinite Curve{dipole_c()} = n_circ;

Mesh.Algorithm = 6;
Mesh.ElementOrder = 1;
Mesh.SecondOrderLinear = 1;
