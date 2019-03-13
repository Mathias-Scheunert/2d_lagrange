Point(1) = {-5000, 0, 0, 1000};
Point(2) = {6000, 0, 0, 1000};
Point(3) = {-5000, 0, 25, 1000};
Point(4) = {6000, 0, 25, 1000};
Point(5) = {-5000, 0, 5500, 2500};
Point(6) = {6000, 0, 5500, 2500};
Point(7) = {-500, 0, 0, 1};
Point(8) = {1000, 0, 0, 1};
Point(9) = {-500, 0, 25, 5};
Point(10) = {1000, 0, 25, 5};
Point(11) = {-500, 0, 250, 10};
Point(12) = {1000, 0, 250, 10};

Line(1) = {1, 3};
Line(2) = {2, 4};
Line(3) = {3, 5};
Line(4) = {4, 6};
Line(5) = {7, 9};
Line(6) = {8, 10};
Line(7) = {9, 11};
Line(8) = {10, 12};
Line(9) = {1, 7};
Line(10) = {8, 2};
Line(11) = {3, 9};
Line(12) = {10, 4};
Line(13) = {7, 8};
Line(14) = {9, 10};
Line(15) = {11, 12};
Line(16) = {5, 6};

Line Loop(1) = {1, 11, -5, -9};
Plane Surface(1) = {1};
Line Loop(2) = {5, 14, -6, -13};
Plane Surface(2) = {2};
Line Loop(3) = {6, 12, -2, -10};
Plane Surface(3) = {3};

Line Loop(4) = { 3, 16, -4, -12, 8, -15, -7, -11};
Plane Surface(4) = {4};
Line Loop(5) = { 8, -15, -7, 14};
Plane Surface(5) = {5};

Physical Line("xmin") = { 1, 3};
Physical Line("xmax") = { 2, 4};
Physical Line("ymin") = {9, 13, 10};
Physical Line("ymax") = {16};

Physical Surface("schicht_1") = { 1, 2, 3};
Physical Surface("schicht_2") = {4, 5};


