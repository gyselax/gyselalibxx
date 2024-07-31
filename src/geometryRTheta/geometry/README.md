# Geometry RTheta

The geometry folder contains a helper file `geometry.hpp` which provides shortcuts to the types needed to define the specific geometry. In this case the geometry is a simulation with 
* 2 spatial dimensions (denoted `R` and `P`) in polar coordinates; 
* 2 spatial dimensions (denoted `X` adn `Y`) in Cartesian coordinates.

We suppose $`r`$ positive non periodic, $`\theta`$ periodic, $`x`$ and $`y`$ non periodic. We denote in the code `P` the polar coordinate $`\theta`$. 
The polar coordinates represent a **logical domain**, and the Cartesian coordinates the **physical domain**. Some mappings from the logical domain to the physical domain $`(r,\theta) \mapsto (x,y)`$ are defined in Curvilinear2DToCartesian. 

The splines applied are cubic and only defined on the logical domain. 

The polar geometry can have a degenerated point. We call *O-point* or *center point*, the points $`(r=0, \theta) \ \forall \theta`$. At these points, the loss of some properties pushes us to define another B-spline basis. This new spline is called in the code *polar spline* (see PolarSpline). Both basis are used in the code. 

## Shortcuts
The shortcuts defined in the geometry file represent:

1. The continuous dimensions `RDimR`, `RDimP`, `RDimX` and `RDimY`. 
    Continuous dimensions for the velocity on the logical and physical domains are also defined `RDimVr`, `RDimVp`, `RDimVx` and `RDimVy`.
2. The type of coordinates of objects represented on the dimensions (`CoordR`, `CoordP`, `CoordRP`, `CoordVr`, `CoordVp`, `CoordX`, `CoordY`, `CoordXY`, `CoordVx` and `CoordVy`).
3.  The type of the B-Spline bases used on the logical dimensions (`BSplinesR` and `BSplinesP`).
4.  The type which will describe the grid points on which the simulation will evolve (`IDimR`, `IDimP`). 
5.  The type of the helper class which initialises grid points in space which are compatible with the defined splines (`SplineInterpPointsR`,  `SplineInterpPointsP`).
6.  The type of the objects used to build splines and evaluate them (`SplineRPBuilder`, `SplineRPEvaluator`).
7.  The type which describes the index of a grid point (e.g. `IndexR`).
8.  The type which describes a distance between grid points (e.g. `IVectR`).
9.  The type which describes the domain on which the grid points are defined (e.g. `IDomainR`).
10. The templated type of a field defined on each of the domains (e.g. `FieldR<ElementType>`). A field defines values at each grid point in the provided domain.
11. The type of a field of doubles defined on each of the domains (e.g. `DFieldR`).
12. The templated type of a field span defined on each of the domains (e.g. `SpanR<ElementType>`). A span is similar to a reference. It provides access to a field or a sub-set of a field, but does not own the data itself.
13. The type of a span of doubles defined on each of the domains (e.g. `DSpanR`).
14. The templated type of a field view defined on each of the domains (e.g. `ViewR<ElementType>`). A view is very similar to a span. The only difference is that a view is constant. In other words, the data in the underlying field cannot be modified using a view.
15. The type of a view of doubles defined on each of the domains (e.g. `DViewR`).
16. The type of a field of doubles defined on the spline domains representing a spline on the Cartesian product base (e.g. `Spline2D`), and a spline on the polar spline base (e.g. `SplinePolar`).
17. The type of VectorField defined on the domain `IDomainRP` on the template directions `NDTag<Dim1, Dim2>` (`VectorDFieldRP<Dim1, Dim2>`). The directions used in the code are either `<RDimR, RDimP>` or `<RDimX, RDimY>`. 
