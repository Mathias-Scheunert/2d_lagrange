// Gmsh project created on Mon Aug 20 09:01:27 2018
SetFactory("OpenCASCADE");

// Define geometries.
Circle(1) = {0, 0, 0, 10, 0, 2*Pi};
Line Loop(1) = {1};
Plane Surface(1) = {1};

Circle(2) = {0, 0, 0, 2.5, 0, 2*Pi};
Line Loop(2) = {2};
Plane Surface(2) = {2};

Rectangle(3) = {-11, 0, 0, 22, 11, 0};

BooleanDifference{ Surface{1}; Delete;}{ Surface{3};}
BooleanDifference{ Surface{2}; Delete;}{ Surface{3}; Delete;}
BooleanDifference{ Surface{1}; Delete;}{ Surface{2};}

Rotate {{1, 0, 0}, {0, 0, 0}, -Pi/2} {
  Surface{1,2}; 
}

Point(15) = {1, 0, 1, 0.05};
Point{15} In Surface{2};

// Define physical entities.
Physical Line("surface") = {15, 26, 28};
Physical Line("subsurface") = {27};
Physical Line("interior") = {12};

Physical Surface("inner_part") = {2};
Physical Surface("outer_part") = {1};

// Add some complexity to model.

Point(16) = {-4.7, -0, 3.6, 2.0};
Point(17) = {-4.3, -0, 2.8, 2.0};
Point(18) = {-3.7, -0, 2.6, 2.0};
Point(19) = {-3.3, 0, 2.8, 2.0};
Point(20) = {-3.2, 0, 3.1, 2.0};
Point(21) = {-2.8, 0, 3.1, 2.0};
Point(22) = {-2.7, 0, 3.3, 2.0};
Point(23) = {-2.6, 0, 3.8, 2.0};
Point(24) = {-2.8, 0, 4, 2.0};
Point(25) = {-3.9, 0, 3.8, 2.0};

Line(16) = {16, 17};
Line(17) = {17, 18};
Line(18) = {18, 19};
Line(19) = {19, 20};
Line(20) = {20, 21};
Line(21) = {21, 22};
Line(22) = {22, 23};
Line(23) = {23, 24};
Line(24) = {24, 25};
Line(25) = {25, 16};

Line Loop(3) = {17, 18, 19, 20, 21, 22, 23, 24, 25, 16};
Plane Surface(3) = {3};
BooleanDifference{ Surface{1}; Delete;}{ Surface{3}; Delete;}

Physical Line("tunnel") = {17, 18, 20, 19, 22, 23, 24, 25, 16, 21};

Point(28) = {0, 0, 0, 1.0};
Point(29) = {5, 0, 0, 1.0};

Field[1] = Attractor;
Field[1].NodesList = {28, 29};

Field[2] = Threshold;
Field[2].IField = 1;
Field[2].LcMin = 0.1;
Field[2].DistMin = 0;
Field[2].DistMax = 0.2;
Background Field = 2;

Field[3] = MathEval;
Field[3].F = "0.2";

Field[4] = Restrict;
Field[4].IField = 3;
Field[4].FacesList = {2};
Field[4].EdgesList = {12, 15};

Field[5] = MathEval;
Field[5].F = "1";

Field[6] = Restrict;
Field[6].IField = 5;
Field[6].EdgesList = {26, 27, 28};
Field[6].FacesList = {1};

Field[7] = Min;
Field[7].FieldsList = {2, 4, 6};
Background Field = 7;
