// Fast 2D Lorentz-profile regression case.
// x in [0, 2], y in [0, 1].
// Volume 1: vacuum, x in [0, 0.5] (scattered-field side)
// Volume 2: Lorentz profile, x in [0.5, 1.0]
// Volume 3: Lorentz profile, x in [1.0, 1.5]
// Volume 4: vacuum, x in [1.5, 2.0]
// TFSF entry between volumes 1 and 2 (x = 0.5), +x, Ey/Hz.
// SMA on the left (x = 0); PEC on the right (x = 2) and on top/bottom.
SetFactory("Built-in");
Mesh.MshFileVersion = 2.2;
Mesh.CharacteristicLengthMin = 0.08;
Mesh.CharacteristicLengthMax = 0.2;

Point(1) = {0, 0, 0};
Point(2) = {0.5, 0, 0};
Point(3) = {1.0, 0, 0};
Point(4) = {1.5, 0, 0};
Point(5) = {2.0, 0, 0};
Point(6) = {0, 1, 0};
Point(7) = {0.5, 1, 0};
Point(8) = {1.0, 1, 0};
Point(9) = {1.5, 1, 0};
Point(10) = {2.0, 1, 0};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 10};
Line(6) = {10, 9};
Line(7) = {9, 8};
Line(8) = {8, 7};
Line(9) = {7, 6};
Line(10) = {6, 1};
Line(11) = {2, 7};
Line(12) = {3, 8};
Line(13) = {4, 9};

Curve Loop(1) = {1, 11, 9, 10};
Plane Surface(1) = {1};
Curve Loop(2) = {2, 12, 8, -11};
Plane Surface(2) = {2};
Curve Loop(3) = {3, 13, 7, -12};
Plane Surface(3) = {3};
Curve Loop(4) = {4, 5, 6, -13};
Plane Surface(4) = {4};

Physical Surface(1) = {1};
Physical Surface(2) = {2};
Physical Surface(3) = {3};
Physical Surface(4) = {4};

Physical Line(1) = {1};
Physical Line(2) = {2};
Physical Line(3) = {3};
Physical Line(4) = {4};
Physical Line(5) = {5};
Physical Line(6) = {6};
Physical Line(7) = {7};
Physical Line(8) = {8};
Physical Line(9) = {9};
Physical Line(10) = {10};
Physical Line(11) = {11};
Physical Line(12) = {12};
Physical Line(13) = {13};
