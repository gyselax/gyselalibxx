

# File polarpoissonlikesolver.hpp



[**FileList**](files.md) **>** [**pde\_solvers**](dir_be2a347b8fed8e825bae8c199ecc63c1.md) **>** [**polarpoissonlikesolver.hpp**](polarpoissonlikesolver_8hpp.md)

[Go to the source code of this file](polarpoissonlikesolver_8hpp_source.md)



* `#include <ddc/ddc.hpp>`
* `#include "ddc_alias_inline_functions.hpp"`
* `#include "ddc_aliases.hpp"`
* `#include "gauss_legendre_integration.hpp"`
* `#include "mapping_tools.hpp"`
* `#include "math_tools.hpp"`
* `#include "matrix_batch_csr.hpp"`
* `#include "metric_tensor_evaluator.hpp"`
* `#include "polar_spline.hpp"`
* `#include "polar_spline_evaluator.hpp"`
* `#include "quadrature_coeffs_nd.hpp"`
* `#include "view.hpp"`
* `#include "volume_quadrature_nd.hpp"`















## Classes

| Type | Name |
| ---: | :--- |
| class | [**PolarSplineFEMPoissonLikeSolver**](classPolarSplineFEMPoissonLikeSolver.md) &lt;class [**GridR**](structGridR.md), class [**GridTheta**](structGridTheta.md), class [**PolarBSplinesRTheta**](structPolarBSplinesRTheta.md), class SplineRThetaEvaluatorNullBound, class IdxRangeFull&gt;<br>_Define a polar PDE solver for a Poisson-like equation._  |
| struct | [**EvalDeriv1DType**](structPolarSplineFEMPoissonLikeSolver_1_1EvalDeriv1DType.md) <br>_Object storing a value and a value of the derivative of a 1D function._  |
| struct | [**EvalDeriv2DType**](structPolarSplineFEMPoissonLikeSolver_1_1EvalDeriv2DType.md) <br>_Object storing a value and a value of the derivatives in each direction of a 2D function._  |
| struct | [**QDimRMesh**](structPolarSplineFEMPoissonLikeSolver_1_1QDimRMesh.md) <br>_Tag the first dimension for the quadrature mesh._  |
| struct | [**QDimThetaMesh**](structPolarSplineFEMPoissonLikeSolver_1_1QDimThetaMesh.md) <br>_Tag the second dimension for the quadrature mesh._  |
| struct | [**RBasisSubset**](structPolarSplineFEMPoissonLikeSolver_1_1RBasisSubset.md) <br> |
| struct | [**RCellDim**](structPolarSplineFEMPoissonLikeSolver_1_1RCellDim.md) <br> |
| struct | [**ThetaBasisSubset**](structPolarSplineFEMPoissonLikeSolver_1_1ThetaBasisSubset.md) <br> |
| struct | [**ThetaCellDim**](structPolarSplineFEMPoissonLikeSolver_1_1ThetaCellDim.md) <br> |



















































------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/pde_solvers/polarpoissonlikesolver.hpp`

