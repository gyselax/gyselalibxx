## Compilation

to compile voice++:

```
git clone --recurse-submodules git@gitlab.maisondelasimulation.fr:gysela-developpers/voicexx.git
cd voicexx
mkdir build
cd build
cmake ..
make
```

to run the tests:
```
ctest --output-on-failure
```

## Dependencies

* cmake-3.15+
* a c++-17 compiler:
  * gcc-9+ should work
* Lapack
* FFTW-3
* mdspan (embedded)
* googletest (embedded)

# Architecture

## Operation classes

Classes are divided in two categories:
* "interface" abstract classes with only public pure virtual functions & a
  virtual destructor,
* implementation concrete classes.

Each interface must offer a pure virtual call operator: `operator()` and can
offer various overloads of this operator with different parameters but similar
behaviour.
Interfaces names should typically be prefixed by `I` as in `IVlasovSolver`.

Implementation classes should implement one interface each .
When relying on other operations, these classes should use the dependency
injection pattern.
Each implementation class should take `const` references to the used operations
interfaces as parameters of the constructor and store them so as to call them in
the `operator()` implementation.

### `splines` & `boundary_conditions`

Spline interpolation (mesh to spline) & evaluation (spline to mesh)

### `time`

Time iteration, including main time loop & predictor-corrector.

### `vlasov`

Vlasov solver including splitting & advections.

### `efield`

TODO field solver.


