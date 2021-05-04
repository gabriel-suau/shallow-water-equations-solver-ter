# SWES1D/SWES2D

C++ solvers for the 1D/2D Saint-Venant equations with topography and without friction using the Finite Volumes Method. 

Work in progress...

The code will be validated with test cases and experimental data.


## Installation, compilation and execution

You can get the sources from here by typing in your terminal :

```shell
git clone https://github.com/gabriel-suau/shallow-water-equations-solver-ter.git SWES
```

To compile the 1D code in release mode, you can type the following commands :

```shell
cd SWES/code_1D
make release
```

This command will produce an executable called <code>main</code>. To execute the program with the parameters written in <code>parameters.txt</code>, you can type :

```shell
./main parameters.txt
```


## Outputs

The outputs of the computation are written in the directory <code>resultsDir</code> specified in <code>parameters.txt</code>. The solution files are named following this rule :

```shell
solution_%FLUXNAME%_%SAVEITERATION%.txt
```

The code also saves the values of the parameters used for the simulation in a file named <code>params.txt</code> and the topography in the file <code>topography.txt</code>. If you simulate a built-in test case, the code will also save the exact solution at the final time of your simulation in the file <code>solution_exacte.txt</code>.

The output files are formatted in columns. Each column correspond to a quantity. For the moment there are 6 columns in the output files, containing : the abscissa, the height of the free surface, the water height, the vertically averaged water velocity, the vertically averaged water discharge and the Froude number. These files can directly be plotted using <code>gnuplot</code>. For example, let's say you ran a simulation using the HLL numerical flux. To visualise the solution at the save iteration n°10, you can type :

```shell
gnuplot
plot "solution_HLL_10.txt"
```

You may also want to visualise all your results file in the form of an animation. To do that, assuming you have 150 results files, you could type something like :

```shell
gnuplot
do for [i=0:150] {plot "solution_HLL_".i.".txt"; pause 0.02}
```


## Documentation

A documentation for the 1D code can be automatically generated with [Doxygen](https://www.doxygen.nl/index.html). Just type :

```shell
doxygen doxygen.cfg
```

This command will produce a <code>doc/</code> directory containing an HTML documentation in the <code>html/</code> directory and a PDF documentation in the <code>latex</code> directory. To read the HTML documentation using <code>firefox</code>, just type the following command (assuming you are in the root directory)

```shell
firefox doc/html/index.html
```

There is currently no Doxygen documentation for the 2D code.

## Credits

All developpers are students at ENSEIRB-MATMECA, a french engineering school located in Talence.

* Robin Colombier
* Geoffrey Lebaud
* Rémi Pégouret
* Gabriel Suau
* Lucas Trautmann


## License

This project is distributed under the [GNU-GPLv3](https://www.gnu.org/licenses/gpl-3.0.html) license. A copy of the whole license is included in the repository.


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
