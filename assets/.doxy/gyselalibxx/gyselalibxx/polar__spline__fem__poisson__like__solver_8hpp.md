

# File polar\_spline\_fem\_poisson\_like\_solver.hpp



[**FileList**](files.md) **>** [**pde\_solvers**](dir_be2a347b8fed8e825bae8c199ecc63c1.md) **>** [**polar\_spline\_fem\_poisson\_like\_solver.hpp**](polar__spline__fem__poisson__like__solver_8hpp.md)

[Go to the source code of this file](polar__spline__fem__poisson__like__solver_8hpp_source.md)



* `#include <Kokkos_Profiling_ScopedRegion.hpp>`
* `#include "i_interpolation.hpp"`
* `#include "ipolar_poisson_like_solver.hpp"`
* `#include "polar_spline_fem_poisson_like_assembler.hpp"`















## Classes

| Type | Name |
| ---: | :--- |
| class | [**PolarSplineFEMPoissonLikeSolver**](classPolarSplineFEMPoissonLikeSolver.md) &lt;class [**GridR**](structGridR.md), class [**GridTheta**](structGridTheta.md), class [**PolarBSplinesRTheta**](structPolarBSplinesRTheta.md), Interpolation2D, class Mapping, MetricTensorType, class IdxRangeFull&gt;<br>_Define a polar PDE solver for a Poisson-like equation._  |
| class | [**CoeffEvaluator**](classPolarSplineFEMPoissonLikeSolver_1_1CoeffEvaluator.md) &lt;class Evaluator, class Coeff&gt;<br>_A wrapper that binds an evaluator with its batched coefficient field to present a single callable_ `double operator()(IdxBatch, CoordRTheta)` _._ |
| struct | [**QDimRMesh**](structPolarSplineFEMPoissonLikeSolver_1_1QDimRMesh.md) <br>_Tag the first dimension for the quadrature mesh._  |
| struct | [**QDimThetaMesh**](structPolarSplineFEMPoissonLikeSolver_1_1QDimThetaMesh.md) <br>_Tag the second dimension for the quadrature mesh._  |



















































------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/pde_solvers/polar_spline_fem_poisson_like_solver.hpp`

