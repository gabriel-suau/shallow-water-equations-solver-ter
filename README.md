# TER

C++ solver for the 1D/2D Saint-Venant equations with topography and without friction using the Finite Volumes Method. 

Work in progress...

The code will be validated with test cases and experimental data.

## Credits

All developpers are students at ENSEIRB-MATMECA, a french engineering school located in Talence.

* Robin Colombier
* Théo Guichard
* Geoffrey Lebaud
* Rémi Pégouret
* Gabriel Suau
* Lucas Trautmann

## Check-list 1D
- [x] Organize the code and create a Makefile.
- [x] Implement basic 1D Finite Volumes schemes.
- [x] Verification (Dam Break problem).
- [ ] Validation (experimental data).

## Check-list 2D
- [x] Organize the code and create a Makefile.
- [ ] Implement basic 2D Finite Volumes schemes.
- [ ] Verification (Dam Break problem).
- [ ] Validation (experimental data).

## Main features
* Numerical fluxes : Lax-Friedrichs, Rusanov and HLL.
* Time scheme : Explicit Euler for both the advection term and the source term.
* Source term : Topography (not working yet).
* Boundary conditions : not working either.
* Scenarios : Resting lake, Dam Break, Sine Perturbation...

## Additional optional features
- [ ] Implement more robust Finite Volumes schemes (well-balanced schemes) ?
- [ ] Parallelize the code ?
- [ ] Add a friction source term ?
