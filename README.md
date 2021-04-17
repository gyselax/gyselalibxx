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

## Dependencies

* cmake-3.15+
* a c++-17 compiler:
  * gcc-9 should work

# Architecture

Voice++ is designed using an imperative architecture with:
* data classes with no to minimal behaviour and no polymorphism,
* operation classes with no state.

## Data classes

### `array_tools`

A simple wrapper around std::experimental::basic_mdspan to ensure better type
replacability by always using layout_stride (no layout_left or layout_right).

### `math_tools`

Simple tensor manipulation on top of `array_tools`.

### `matrix`

Various in-house sparse matrix tools.

### `geometry`

Description of the geometry (real domain).

### `mesh`

Discretization of the real domain into a discrete domain (mesh) with the
possibility to associate values to points of this mesh.

### `splines` & `boundary_conditions`

Another possible discretization based on a spline representation instead of a
mesh.

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

### `time`

Time iteration, including main time loop & predictor-corrector.

### `vlasov`

Vlasov solver including splitting & advections.

### `efield`

TODO field solver.
