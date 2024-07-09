# Geometry XY

The geometry folder contains a helper file `geometry.hpp` which provides shortcuts to the types needed to define the specific geometry. In this case the geometry is a simulation with 2 spatial dimensions (denoted `X` adn `Y`).

In this geometry $`x`$ and $`y`$ are supposed **periodic** and the splines are supposed cubic. 

The shortcuts defined in the geometry file represent:
1. The continuous dimensions `RDimX` and `RDimY`. 
2. The type of coordinates of objects represented on the dimensions (`CoordX`, `CoordY`, `CoordXY`).
3. The type of the B-Spline bases used on the spatial dimensions (`BSplinesX` and `BSplinesY`). 
4.  The type which will describe the grid points on which the simulation will evolve (`IDimX`, `IDimY`).
5.  The type of the helper class which initialises grid points in space which are compatible with the defined splines (`SplineInterpPointsX`,  `SplineInterpPointsY`).
6.  The type of the objects used to build splines (`SplineXBuilder_XY`, `SplineYBuilder_XY`).
7.  The type which describes the index of a grid point (e.g. `IndexX`).
8.  The type which describes a distance between grid points (e.g. `IVectX`).
9.  The type which describes the domain on which the grid points are defined (e.g. `IDomainX`).
10. The templated type of a field defined on each of the domains (e.g. `FieldX<ElementType>`). A field defines values at each grid point in the provided domain.
11. The type of a field of doubles defined on each of the domains (e.g. `DFieldX`).
12. The templated type of a field span defined on each of the domains (e.g. `SpanX<ElementType>`). A span is similar to a reference. It provides access to a field or a sub-set of a field, but does not own the data itself.
13. The type of a span of doubles defined on each of the domains (e.g. `DSpanX`).
12. The templated type of a field view defined on each of the domains (e.g. `ViewX<ElementType>`). A view is very similar to a span. The only difference is that a view is constant. In other words, the data in the underlying field cannot be modified using a view.
13. The type of a view of doubles defined on each of the domains (e.g. `DViewX`).
14. The type of VectorField defined on the domain `IDomainXY` on the directions `NDTag<RDimX, RDimY>` (`VectorFieldXY_XY`). 
