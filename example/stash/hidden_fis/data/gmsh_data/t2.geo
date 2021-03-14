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

Point(1) = {0., 0., 0., lc};
Point(2) = {1., 0., 0., lc} ;
Point(3) = {1., 1., 0., lc} ;
Point(4) = {0., 1., 0., lc} ;

Line(1) = {1, 2} ;
Line(2) = {2, 3} ;
Line(3) = {3, 4} ;
Line(4) = {4, 1} ;
Line(5) = {3, 1} ;


Line Loop(1) = {1, 2, 5} ;
Line Loop(2) = {-5, 3, 4} ;


Plane Surface(1) = {1} ;
Plane Surface(2) = {2} ;

Physical Line(1) = {1, 2, 3, 4} ;
Physical Surface(1) = {1} ;
Physical Surface(2) = {2} ;



