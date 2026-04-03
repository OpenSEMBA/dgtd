// 3D RCS SGBC Sphere in Box — Reduced size
// Sphere r=0.5m, TFSF box l=1.8m, SMA sphere r=2.5m
SetFactory("OpenCASCADE");

// --- Geometry ---
Sphere(1) = {0, 0, 0, 0.5, -Pi/2, Pi/2, 2*Pi};
Sphere(2) = {0, 0, 0, 2.5, -Pi/2, Pi/2, 2*Pi};
l = 1.8;
Box(3) = {-l/2, -l/2, -l/2, l, l, l};

Delete { Volume{1}; }
Delete { Volume{2}; }
Delete { Volume{3}; }

// --- Volumes ---
Surface Loop(4) = {2};
Surface Loop(5) = {6, 7, 3, 5, 8, 4};
Volume(1) = {4, 5};          // between SMA sphere and TFSF box

Surface Loop(6) = {6, 7, 3, 5, 8, 4};
Surface Loop(7) = {1};
Volume(2) = {6, 7};          // between TFSF box and SGBC sphere

Surface Loop(8) = {1};
Volume(3) = {8};              // inside SGBC sphere (no fields)

// --- Physical groups ---
Physical Volume("vacuum", 1) = {1, 2, 3};
Physical Surface("SMA", 2)   = {2};
Physical Surface("TFSF", 3)  = {6, 4, 8, 3, 7, 5};
Physical Surface("SGBC", 4)  = {1};

// --- Mesh sizing ---
// At fc=150 MHz, λ = 2.0 m.  Order 2 needs ~5 pts/λ → h ≈ 0.40 m max.
// Use finer sizes on TFSF and scatterer for accuracy.

// Sphere surface: h = 0.15 m (curvature + scattering resolution)
MeshSize{ PointsOf{ Surface{1}; } } = 0.15;

// TFSF box faces: h = 0.15 m (clean wave injection)
MeshSize{ PointsOf{ Surface{3, 4, 5, 6, 7, 8}; } } = 0.15;

// SMA sphere: h = 0.40 m (coarse, only absorbs outgoing scattered field)
MeshSize{ PointsOf{ Surface{2}; } } = 0.40;

// Grade from TFSF box outward to SMA
Field[1] = Distance;
Field[1].SurfacesList = {3, 4, 5, 6, 7, 8};

Field[2] = Threshold;
Field[2].InField = 1;
Field[2].SizeMin = 0.15;   // at TFSF faces
Field[2].SizeMax = 0.40;   // approaching SMA
Field[2].DistMin = 0.0;
Field[2].DistMax = 1.5;

// Grade from sphere surface inward (sphere interior is dead — SGBC+PEC)
Field[3] = Distance;
Field[3].SurfacesList = {1};

Field[4] = Threshold;
Field[4].InField = 3;
Field[4].SizeMin = 0.15;   // on sphere surface
Field[4].SizeMax = 0.35;   // sphere interior bulk
Field[4].DistMin = 0.0;
Field[4].DistMax = 0.3;

// Keep interior of TFSF box fine (scattered field region)
Field[5] = Box;
Field[5].VIn  = 0.15;
Field[5].VOut = 0.40;
Field[5].XMin = -l/2;  Field[5].XMax = l/2;
Field[5].YMin = -l/2;  Field[5].YMax = l/2;
Field[5].ZMin = -l/2;  Field[5].ZMax = l/2;
Field[5].Thickness = 0.2;

Field[6] = Min;
Field[6].FieldsList = {2, 4, 5};

Background Field = 6;

Mesh.Algorithm3D = 10;      // HXT
Mesh.OptimizeNetgen = 1;
