

# Class GMGPolarTools::PolarPoissonLikeCoefficients

**template &lt;class SplineEvaluator, class [**BSplinesR**](structBSplinesR.md), class [**BSplinesTheta**](structBSplinesTheta.md)&gt;**



[**ClassList**](annotated.md) **>** [**GMGPolarTools**](namespaceGMGPolarTools.md) **>** [**PolarPoissonLikeCoefficients**](classGMGPolarTools_1_1PolarPoissonLikeCoefficients.md)



_Wraps gyselalibxx spline-represented coefficients to satisfy the GMGPolar DensityProfileCoefficients concept._ [More...](#detailed-description)

* `#include <gmg_polar_poisson_like_solver.hpp>`





































## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**PolarPoissonLikeCoefficients**](#function-polarpoissonlikecoefficients) (SplineEvaluator evaluator, DConstSplineRTheta coeff\_alpha, DConstSplineRTheta coeff\_beta) <br>_Build the class instance._  |
|  KOKKOS\_INLINE\_FUNCTION double | [**alpha**](#function-alpha) (const double & r, const double & theta) const<br>_The coefficient alpha in the Poisson-like equation._  |
|  KOKKOS\_INLINE\_FUNCTION double | [**beta**](#function-beta) (const double & r, const double & theta) const<br>_The coefficient beta in the Poisson-like equation._  |


## Public Static Functions

| Type | Name |
| ---: | :--- |
|  double | [**getAlphaJump**](#function-getalphajump) () <br>_Required for the concept, only used in custom mesh generation (refinement\_radius); not needed here._  |


























## Detailed Description




**Template parameters:**


* `SplineEvaluator` A 2D spline evaluator for ([**BSplinesR**](structBSplinesR.md), [**BSplinesTheta**](structBSplinesTheta.md)). 
* [**BSplinesR**](structBSplinesR.md) The radial B-spline type. 
* [**BSplinesTheta**](structBSplinesTheta.md) The poloidal B-spline type. 




    
## Public Functions Documentation




### function PolarPoissonLikeCoefficients 

_Build the class instance._ 
```C++
inline GMGPolarTools::PolarPoissonLikeCoefficients::PolarPoissonLikeCoefficients (
    SplineEvaluator evaluator,
    DConstSplineRTheta coeff_alpha,
    DConstSplineRTheta coeff_beta
) 
```




<hr>



### function alpha 

_The coefficient alpha in the Poisson-like equation._ 
```C++
inline KOKKOS_INLINE_FUNCTION double GMGPolarTools::PolarPoissonLikeCoefficients::alpha (
    const double & r,
    const double & theta
) const
```




<hr>



### function beta 

_The coefficient beta in the Poisson-like equation._ 
```C++
inline KOKKOS_INLINE_FUNCTION double GMGPolarTools::PolarPoissonLikeCoefficients::beta (
    const double & r,
    const double & theta
) const
```




<hr>
## Public Static Functions Documentation




### function getAlphaJump 

_Required for the concept, only used in custom mesh generation (refinement\_radius); not needed here._ 
```C++
static inline double GMGPolarTools::PolarPoissonLikeCoefficients::getAlphaJump () 
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/pde_solvers/gmg_polar_poisson_like_solver.hpp`

