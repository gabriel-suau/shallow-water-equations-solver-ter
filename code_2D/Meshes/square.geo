// À adapter en fonction du raffinement de maillage souhaité
h = 0.1;

// Square properties
length = 10.;

// Points
Point(1) = {-0.5*length,-0.5*length,0,h};
Point(2) = {+0.5*length,-0.5*length,0,h};
Point(3) = {+0.5*length,+0.5*length,0,h};
Point(4) = {-0.5*length,+0.5*length,0,h};

// Lines
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};
Line Loop(1) = {1:4};
Plane Surface(1) = {1};