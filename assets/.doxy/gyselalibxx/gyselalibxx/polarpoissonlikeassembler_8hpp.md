

# File polarpoissonlikeassembler.hpp



[**FileList**](files.md) **>** [**pde\_solvers**](dir_be2a347b8fed8e825bae8c199ecc63c1.md) **>** [**polarpoissonlikeassembler.hpp**](polarpoissonlikeassembler_8hpp.md)

[Go to the source code of this file](polarpoissonlikeassembler_8hpp_source.md)



* `#include <ddc/ddc.hpp>`
* `#include "coord_transformation_tools.hpp"`
* `#include "ddc_alias_inline_functions.hpp"`
* `#include "ddc_aliases.hpp"`
* `#include "gauss_legendre_integration.hpp"`
* `#include "math_tools.hpp"`
* `#include "matrix_batch_csr.hpp"`
* `#include "metric_tensor_evaluator.hpp"`
* `#include "polar_spline_evaluator.hpp"`
* `#include "quadrature_coeffs_nd.hpp"`
* `#include "view.hpp"`
* `#include "volume_quadrature_nd.hpp"`













## Namespaces

| Type | Name |
| ---: | :--- |
| namespace | [**detail\_poisson**](namespacedetail__poisson.md) <br> |


## Classes

| Type | Name |
| ---: | :--- |
| class | [**PolarSplineFEMPoissonLikeAssembler**](classPolarSplineFEMPoissonLikeAssembler.md) &lt;typename [**GridR**](structGridR.md), typename [**GridTheta**](structGridTheta.md), typename [**PolarBSplinesRTheta**](structPolarBSplinesRTheta.md), typename SplineRThetaEvaluatorNullBound, typename QDimRMesh, typename QDimThetaMesh, class IdxRangeFull&gt;<br>_An operator to assemble a Poisson-like stiffness matrix using polar B-splines._  |
| struct | [**InternalBatchDim**](structPolarSplineFEMPoissonLikeAssembler_1_1InternalBatchDim.md) <br>_The tag for the batch dimension for the equation. This is public due to Cuda._  |



















































------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/pde_solvers/polarpoissonlikeassembler.hpp`

