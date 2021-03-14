/*********************************************************************
 *
 *  Gmsh tutorial 1
 *
 *  Variables, elementary entities (points, lines, surfaces), physical
 *  entities (points, lines, surfaces)
 *
 *********************************************************************/

// The simplest construction in Gmsh's scripting language is the
// `affectation'. The following command defines a new variable `lc':

lc = 2e-1;

Point(1) = {-1, -1., 0., lc};
Point(2) = {0.5, -1., 0., lc};
Point(3) = {1., 0., 0., lc};
Point(4) = {0.5, 1., 0., lc};
Point(5) = {-0.5, 1., 0., lc};
Point(6) = {-1., 0., 0., lc};
Point(0) = {0., 0., 0., lc};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 1};
Line(7) = {1, 0};
Line(8) = {2, 0};
Line(9) = {3, 0};
Line(10) = {4, 0};
Line(11) = {5, 0};
Line(12) = {6, 0};

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

Line Loop(1) = {1, 8, -7};
Line Loop(2) = {2, 9, -8};
Line Loop(3) = {3, 10, -9};
Line Loop(4) = {4, 11, -10};
Line Loop(5) = {5, 12, -11};
Line Loop(6) = {6, 7, -12};

Plane Surface(1) = {1};
Plane Surface(2) = {2};
Plane Surface(3) = {3};
Plane Surface(4) = {4};
Plane Surface(5) = {5};
Plane Surface(6) = {6};

Physical Surface(1) = {1};
Physical Surface(2) = {2};
Physical Surface(3) = {3};
Physical Surface(4) = {4};
Physical Surface(5) = {5};
Physical Surface(6) = {6};
