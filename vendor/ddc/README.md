# The Domain DeComposition library (DDC)

## Compilation

to compile DDC:

```
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
* mdspan
* a c++-17 compiler:
  * gcc-9 should work

# Architecture

## Data classes

### `geometry`

Description of the geometry (real domain).

### `mesh`

Discretization of the real domain into a discrete domain (mesh) with the possibility to associate
values to points of this mesh.

### `splines` & `boundary_conditions`

Another possible discretization based on a spline representation instead of a mesh.
