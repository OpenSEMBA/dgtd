// Coated PEC sphere, SGBC at r = 0.50 m (coating exterior).
// The ball r < 0.50 stays in the mesh, coarse, so the coating is an
// interior face. The PEC is the inner end of the SGBC layer, not a meshed wall.
// TFSF and NTFF are axis-aligned boxes. SMA stays a sphere.
// Order-2 curvilinear tets. JSON tags are elementary ids, filled by bind_case.py.
SetFactory("OpenCASCADE");

r_sgbc = 0.50;
a_tfsf = 0.80;
a_ntff = 1.15;
r_sma  = 2.60;

Sphere(1) = {0, 0, 0, r_sma};
Box(2) = {-a_ntff, -a_ntff, -a_ntff, 2*a_ntff, 2*a_ntff, 2*a_ntff};
Box(3) = {-a_tfsf, -a_tfsf, -a_tfsf, 2*a_tfsf, 2*a_tfsf, 2*a_tfsf};
Sphere(4) = {0, 0, 0, r_sgbc};

BooleanFragments{ Volume{1}; Delete; }{ Volume{2}; Volume{3}; Volume{4}; Delete; }

// Keep the ball. Its bounding box is the only one inside r = 0.55.
vols[] = Volume{:};
inner = -1;
outer[] = {};
For i In {0:#vols()-1}
  bb() = BoundingBox Volume{vols(i)};
  rx = (bb(3) - bb(0)) / 2;
  ry = (bb(4) - bb(1)) / 2;
  rz = (bb(5) - bb(2)) / 2;
  rmax = rx;
  If (ry > rmax)
    rmax = ry;
  EndIf
  If (rz > rmax)
    rmax = rz;
  EndIf
  If (rmax < 0.55)
    inner = vols(i);
  Else
    outer[] += {vols(i)};
  EndIf
EndFor
If (inner < 0)
  Error("inner ball volume was not found");
EndIf
Printf("inner volume %g, outer volumes %g", inner, #outer());

surfs[] = Surface{:};
Physical Surface("boundaries", 1) = {surfs[]};
vac[] = Volume{:};
Physical Volume("vacuum", 2) = {vac[]};
Printf("Coated sphere: %g volumes, %g surfaces", #vac(), #surfs());

// h = 0.08 on the coating, 0.09 through the NTFF box, 0.40 at the SMA.
// The NTFF box corner sits near r = 2, so the fine Box field covers that box
// and the Ball field only tightens the sphere.
Field[1] = Ball;
Field[1].XCenter = 0;
Field[1].YCenter = 0;
Field[1].ZCenter = 0;
Field[1].Radius = 0.55;
Field[1].Thickness = 0.15;
Field[1].VIn = 0.08;
Field[1].VOut = 1.0;

Field[2] = Box;
Field[2].XMin = -a_ntff;
Field[2].XMax = a_ntff;
Field[2].YMin = -a_ntff;
Field[2].YMax = a_ntff;
Field[2].ZMin = -a_ntff;
Field[2].ZMax = a_ntff;
Field[2].VIn = 0.09;
Field[2].VOut = 0.40;
Field[2].Thickness = 0.70;

Field[3] = Min;
Field[3].FieldsList = {1, 2};

// Coarse fill of the inner ball. The sphere surface stays at the outer size.
Field[4] = Constant;
Field[4].VIn = 0.30;

Field[5] = Restrict;
Field[5].InField = 3;
Field[5].VolumesList = {outer[]};

Field[6] = Restrict;
Field[6].InField = 4;
Field[6].VolumesList = {inner};

Field[7] = Min;
Field[7].FieldsList = {5, 6};
Background Field = 7;

Mesh.MeshSizeExtendFromBoundary = 0;
Mesh.MeshSizeFromPoints = 0;
Mesh.MeshSizeFromCurvature = 0;

Mesh.Algorithm3D = 10;
Mesh.Optimize = 1;
Mesh.OptimizeNetgen = 1;
Mesh.ElementOrder = 2;
Mesh.HighOrderOptimize = 2;
Mesh.MshFileVersion = 2.2;
Mesh.Binary = 0;
