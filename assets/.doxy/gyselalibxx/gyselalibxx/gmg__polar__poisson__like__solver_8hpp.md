

# File gmg\_polar\_poisson\_like\_solver.hpp



[**FileList**](files.md) **>** [**pde\_solvers**](dir_be2a347b8fed8e825bae8c199ecc63c1.md) **>** [**gmg\_polar\_poisson\_like\_solver.hpp**](gmg__polar__poisson__like__solver_8hpp.md)

[Go to the source code of this file](gmg__polar__poisson__like__solver_8hpp_source.md)



* `#include <cmath>`
* `#include <vector>`
* `#include <GMGPolar/gmgpolar.h>`
* `#include "ddc_alias_inline_functions.hpp"`
* `#include "i_interpolation_builder.hpp"`
* `#include "ipolar_poisson_like_solver.hpp"`













## Namespaces

| Type | Name |
| ---: | :--- |
| namespace | [**GMGPolarTools**](namespaceGMGPolarTools.md) <br> |


## Classes

| Type | Name |
| ---: | :--- |
| class | [**GMGPolarPoissonLikeSolver**](classGMGPolarPoissonLikeSolver.md) &lt;class ToPhysicalMapping, class [**GridR**](structGridR.md), class [**GridTheta**](structGridTheta.md), class InterpolatorType&gt;<br>_A Poisson-like solver using the GMGPolar multigrid library._  |
| class | [**HomogeneousDirichletBoundaryConditions**](classGMGPolarTools_1_1HomogeneousDirichletBoundaryConditions.md) <br>_Homogeneous Dirichlet boundary conditions satisfying the GMGPolar BoundaryConditions concept._  |
| class | [**MappingToDomainGeometry**](classGMGPolarTools_1_1MappingToDomainGeometry.md) &lt;class ToPhysicalMapping&gt;<br>_Wraps a gyselalibxx coordinate mapping to satisfy the GMGPolar DomainGeometry concept._  |
| class | [**PolarPoissonLikeCoefficients**](classGMGPolarTools_1_1PolarPoissonLikeCoefficients.md) &lt;class EvaluatorType, class IdxRangeCoeff, class CoordRTheta&gt;<br>_Wraps gyselalibxx interpolation-represented coefficients to satisfy the GMGPolar DensityProfileCoefficients concept._  |



















































------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/pde_solvers/gmg_polar_poisson_like_solver.hpp`

