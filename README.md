# SWES1D/SWES2D

C++ solvers for the 1D/2D Saint-Venant equations with topography and without friction using the Finite Volumes Method. 

Work in progress...

The code will be validated with test cases and experimental data.

## Credits

All developpers are students at ENSEIRB-MATMECA, a french engineering school located in Talence.

* Robin Colombier
* Geoffrey Lebaud
* Rémi Pégouret
* Gabriel Suau
* Lucas Trautmann

## Check-list 1D
- [x] Organize the code and create a Makefile.
- [x] Implement basic 1D Finite Volumes schemes.
- [x] Verification (Dam Break problem, stationnary solutions).
- [ ] Validation (experimental data).

## Check-list 2D
- [x] Organize the code and create a Makefile.
- [x] Implement basic 2D Finite Volumes schemes.
- [x] Verification (Dam Break problem).
- [ ] Validation (experimental data).

## Main features
* Numerical fluxes : Lax-Friedrichs, Rusanov and HLL.
* Time scheme : Explicit Euler, Runge-Kutta 2 (Heun).
* Source term : Topography.
* Boundary conditions : subcritical and supercritical boundary conditions (only in the 1D code for now).
* Initial Conditions : Dam break wet/wet and wet/dry, uniform free surface height and discharge...

## Additional optional features
- [ ] Implement more robust Finite Volumes schemes (well-balanced schemes) ?
- [ ] Implement second-order schemes ?
- [ ] Implement more mesh types (make a generic cell class instead of the triangle class)
- [ ] Parallelize the 2D code ?
- [ ] Add a friction source term ?
