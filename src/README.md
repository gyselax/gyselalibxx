# Gyselalib++ contents

The `src/` folder contains all the code necessary to build a gyrokinetic semi-Lagrangian simulation. It is broken up into the following sub-folders:

### Operators for specific geometries: 
- [geometryRTheta](./geometryRTheta/README.md) - Code describing methods which are specific to a 2D curvilinear geometry containing a singular point.
- [geometryVparMu](./geometryVparMu/README.md) - Code describing methods which are specific to a simulation operators in velocity space (vpar,mu).
- [geometryXVx](./geometryXVx/README.md) - Code describing methods which are specific to a simulation with 1 spatial dimension and 1 velocity dimension.
- [geometryXY](./geometryXY/README.md) - Code describing methods which are specific to a simulation with 2 spatial dimensions.
- [geometryXYVxVy](./geometryXYVxVy/README.md) - Code describing methods which are specific to a simulation with 2 spatial dimensions and 2 velocity dimension.


### Operators and math tools for any geometry:
- [advection](./advection/README.md) - Code describing semi-Lagrangian advection routines.
- [collisions](./collisions/README.md) - Code describing the inter and intra species collisions.
- [interpolation](./interpolation/README.md) - Code describing interpolation methods.
- [math\_tools](./math_tools/README.md) - Code describing math tools functions.
- [pde\_solvers](./pde_solvers/README.md) - Code describing different methods for solving PDEs (e.g. Poisson's equation).
- [quadrature](./quadrature/README.md) - Code describing different quadrature methods.
- [speciesinfo](./speciesinfo/README.md) - Code used to describe the different species.
- [timestepper](./timestepper/README.md) - Code used to describe time-stepping methods (e.g. Runge-Kutta methods).

### New types 
- [data\_types](./data_types/README.md) - Code describing useful data types which are not provided by DDC.
- [matrix\_tools](./matrix_tools/README.md) - Code describing both matrix storage and the functions to solve matrix equations.
- [utils](./utils/README.md) - Code describing utility functions.

### Geometry 
- [mapping](./mapping/README.md) - Code describing tools for handling different coordinate systems.
- [multipatch](./multipatch/README.md) - Code describing multi-patch geometry.

### Code development tools
- [io](./io/README.md) - Code describing tools for inputting and outputting data in a simulation.
<!-- - [paraconfpp](./paraconfpp/README.md) - Paraconf utility functions. -->
- [mpi\_parallelisation](./mpi_parallelisation/README.md) - Code describing the MPI parallelisation.

