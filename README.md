# TER

C++ solver for the 1D Saint-Venant equations with topography and without friction using the Finite Volumes Method. 

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

## Check-list
- [x] Organize the code and create a Makefile.
- [x] Implement basic 1D Finite Volumes schemes.
- [ ] Create test cases for verification.
- [x] Create experimental data files for validation.
- [ ] Document the code.

## Main features
* Numerical fluxes : Lax-Friedrichs and Rusanov.
* Time scheme : Explicit Euler for the advection term, Implicit Euler for the source term.
* Source term : Topography (not working yet).
* Boundary conditions : not implemented yet.
* Scenarios : Resting lake, Dam Break

## Additional optional features
- [ ] 2D model.
- [ ] Implement more robust Finite Volumes schemes (well-balanced schemes).
- [ ] Parallelize the code.
- [ ] Add a friction source term.
