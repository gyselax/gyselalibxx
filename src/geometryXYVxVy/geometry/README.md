# Geometry X Y-Vx Vy

The geometry folder contains a helper file `geometry.hpp` which provides shortcuts to the types needed to define the specific geometry. In this case the geometry is a simulation with 1 spatial dimension (denoted X) and 1 velocity dimension (denoted Vx).

The shortcuts defined in the geometry file represent:
1.  The spatial, and velocity dimensions in real space (RDimX, RDimY, RDimVx, and RDimVy).
2.  The type of coordinates of objects represented on the dimensions (`CoordX`, `CoordY`, `CoordXY`, `CoordVx`, `CoordVy`).
3.  The type of the B-Spline bases used on the spatial and velocity dimensions (`BsplinesX`, `BsplinesY`, `BsplinesVx`, `BsplinesVy`,).
4.  The type which will describe the grid points (representing space, velocity and species) on which the simulation will evolve (`IDimX`, `IDimY`, `IDimVx`, `IDimVy`, `IDimSp`).
5.  The type of the helper class which initialises grid points in space and velocity which are compatible with the defined splines (`SplineInterpPointsX`, `SplineInterpPointsY`, `SplineInterpPointsVx`, `SplineInterpPointsVy`).
6.  The type of the objects used to build splines (`SplineXBuilder`, `SplineYBuilder`, `SplineVxBuilder`, `SplineVyBuilder`).
7.  The type which describes the index of a grid point (representing space, velocity and/or species) (e.g. `IndexX`).
8.  The type which describes a distance between grid points (e.g. `IVectX`).
9.  The type which describes the domain on which the grid points are defined (e.g. `IDomainX`).
10. The templated type of a field memory block defined on each of the domains (e.g. `FieldMemX<ElementType>`). A field memory block allocates values at each grid point in the provided domain.
11. The type of a field memory block of doubles defined on each of the domains (e.g. `DFieldMemX`).
12. The templated type of a field defined on each of the domains (e.g. `FieldX<ElementType>`). A field is similar to a reference. It provides access to a field memory block or a sub-set of a field memory block, but does not own the data itself.
13. The type of a field of doubles defined on each of the domains (e.g. `DFieldX`).
14. The templated type of a constant field defined on each of the domains (e.g. `ConstFieldX<ElementType>`).
15. The type of a constant field of doubles defined on each of the domains (e.g. `DConstFieldX`).
16. Types representing coordinates, and the grid points as well as their indices, distances and domains for the Fourier modes.
17. A class GeometryXYVxVy detailing some of the above types in a generic way which allows them to be accessed from a context where the final geometry selected is unknown.
