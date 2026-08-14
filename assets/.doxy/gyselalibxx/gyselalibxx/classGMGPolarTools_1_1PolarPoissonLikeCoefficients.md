

# Class GMGPolarTools::PolarPoissonLikeCoefficients

**template &lt;class EvaluatorType, class IdxRangeCoeff, class CoordRTheta&gt;**



[**ClassList**](annotated.md) **>** [**GMGPolarTools**](namespaceGMGPolarTools.md) **>** [**PolarPoissonLikeCoefficients**](classGMGPolarTools_1_1PolarPoissonLikeCoefficients.md)



_Wraps gyselalibxx interpolation-represented coefficients to satisfy the GMGPolar DensityProfileCoefficients concept._ [More...](#detailed-description)

* `#include <gmg_polar_poisson_like_solver.hpp>`





































## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**PolarPoissonLikeCoefficients**](#function-polarpoissonlikecoefficients) (EvaluatorType evaluator, DConstCoeffRTheta coeff\_alpha, DConstCoeffRTheta coeff\_beta) <br>_Build the class instance._  |
|  KOKKOS\_INLINE\_FUNCTION double | [**alpha**](#function-alpha) (const double & r, const double & theta) const<br>_The coefficient alpha in the Poisson-like equation._  |
|  KOKKOS\_INLINE\_FUNCTION double | [**beta**](#function-beta) (const double & r, const double & theta) const<br>_The coefficient beta in the Poisson-like equation._  |


## Public Static Functions

| Type | Name |
| ---: | :--- |
|  double | [**getAlphaJump**](#function-getalphajump) () <br>_Required for the concept, only used in custom mesh generation (refinement\_radius); not needed here._  |


























## Detailed Description




**Template parameters:**


* `EvaluatorType` A 2D evaluator for the representation described by IdxRangeCoeff. 




    
## Public Functions Documentation




### function PolarPoissonLikeCoefficients 

_Build the class instance._ 
```C++
inline GMGPolarTools::PolarPoissonLikeCoefficients::PolarPoissonLikeCoefficients (
    EvaluatorType evaluator,
    DConstCoeffRTheta coeff_alpha,
    DConstCoeffRTheta coeff_beta
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

