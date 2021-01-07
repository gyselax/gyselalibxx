# GLCS Projet 2020/2021

This is an implementation of a heat equation solver:
* using the forward finite difference method on a Cartesian grid,
* parallelised using MPI,
* supporting the dependency inversion principle via dependency injection.

This project is made up of two sub-libraries each in their own directory, plus
a main code `simpleheat.cpp`.
The sub-libraries are:
* `baselib`: the base classes of the simulation, including:
  - `Distributed2DField`: a 2D array of `double`s distributed in a Cartesian
    fashion,
  - `CartesianDistribution2D`: the description of the MPI distribution used by
    Distributed2DField,
  - `Simulation`: the main simulation class,
* `heatlib`: the implementation of the finite difference heat solver,
* `simpleui`: a very simplistic user-interface taking all parameters from the
  command-line and printing the result on the console.

## Usage

Compilation:
```
mkdir build
cd build
cmake ..
make
```

Execution:
```
./simpleheat
```

## Forward finite difference method on a Cartesian grid

The solver implemented is a straight-forward forward finite difference method on
a Cartesian grid as described in the [provided document](td3-projet.pdf).

## MPI parallelization

The parallelization is mostly supported by the `Distributed2DField` class that
implements a 2D field of `double`s distributed over a Cartesian grid.
This class references a `CartesianDistribution2D` that describes the data
distribution: number of processes used in each dimension (`extents()`), position
of the local process (`coord()`), etc.
In addition, it supports a concept of ghosts where the local blocks overlap at
the frontier between processes.
Each point in the field has a specific process responsible for its computation;
it is in its local no-ghost block.
However if the point is close to the boundary of the local block, it can also be
part of the ghost block of the next process.
The ghosts are synchronized by using the `sync_ghosts()` member function.

## Inversion of control

The dependency inversion principle is supported by using dependency injection.
Whenever a class requires a service provided by another, it does not instantiate
it directly.
Instead it depends on a class with only pure virtual member functions (an
"interface") and it receives the instances implementing this interface from the
outside.
The concrete classes implementing the interfaces can either be provided from
this specific library or in another.
The `main` function acts as the injector, it constructs all concrete classes and
injects these instances in others.
