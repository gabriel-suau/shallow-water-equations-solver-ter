// À adapter en fonction du raffinement de maillage souhaité
h = 0.05;

// Rectangle properties
length = 10.;
height = 1.;

// Points
Point(1) = {-0.5 * length, -0.5*height,0,h};
Point(2) = {0,-0.5*height,0,h};
Point(3) = {0,+0.5*height,0,h};
Point(4) = {-0.5 * length, +0.5*height,0,h};
Point(5) = {+0.5 * length, -0.5*height,0,h};
Point(6) = {+0.5 * length, +0.5*height,0,h};

// Lines
// Left part
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};
// Right part
Line(5) = {2,5};
Line(6) = {5,6};
Line(7) = {6,3};

// Loops
Line Loop(1) = {1,2,3,4};
Line Loop(2) = {5,6,7,-2};

// Surfaces to mesh
Plane Surface(1) = {1};
Plane Surface(2) = {2};

// Physical elements
Physical Curve(1) = {4};
Physical Curve(2) = {1,5};
Physical Curve(3) = {6};
Physical Curve(4) = {7,3};
Physical Surface(1) = {1,2};

// Uncomment for structured mesh
meshTransFinite = 50;
Transfinite Line {1,3,5,7} = 5*meshTransFinite;
Transfinite Line {2,4,6} = meshTransFinite;
Transfinite Surface {1,2};
Recombine Surface {1,2};