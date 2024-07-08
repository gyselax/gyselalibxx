# Gyselalib++ contents

The `src/` folder contains all the code necessary to build a gyrokinetic semi-Lagrangian simulation. It is broken up into the following sub-folders:

- [advection](./advection/README.md) - Code describing semi-Lagrangian advection routines.
- [collision](./collision/README.md) 
- [data\_types](./data_types/README.md)
- [geometry5D](./geometry5D/README.md) - Code describing methods which are specific to a simulation with 3 spatial dimensions and 2 velocity dimension. This will be the general geometry used for GyselaX++ development.
- [geometryRTheta](./geometryRTheta/README.md) - Code describing methods which are specific to a 2D curvilinear geometry containing a singular point.
- [geometryXVx](./geometryXVx/README.md) - Code describing methods which are specific to a simulation with 1 spatial dimension and 1 velocity dimension.
- [geometryXYVxVy](./geometryXYVxVy/README.md) - Code describing methods which are specific to a simulation with 2 spatial dimensions and 2 velocity dimension.
- [interpolation](./interpolation/README.md) - Code describing interpolation methods.
- [io](./io/README.md)
- [multipatch](./multipatch/README.md) - Code describing multipatch geometry.
<!-- - [paraconfpp](./paraconfpp/README.md) - Paraconf utility functions. -->
- [pde\_solvers](./pde_solvers/README.md) - Code describing different methods for solving PDEs (e.g. Poisson's equation).
- [quadrature](./quadrature/README.md) - Code describing different quadrature methods.
- [speciesinfo](./speciesinfo/README.md) - Code used to describe the different species.
- [timestepper](./timestepper/README.md) - Code used to describe time-stepping methods (e.g. Runge-Kutta methods).
- [utils](./utils/README.md) - Code describing utility functions.
