

# Class GMGPolarTools::PolarPoissonLikeCoefficients



[**ClassList**](annotated.md) **>** [**GMGPolarTools**](namespaceGMGPolarTools.md) **>** [**PolarPoissonLikeCoefficients**](classGMGPolarTools_1_1PolarPoissonLikeCoefficients.md)



_Wraps gyselalibxx interpolation-represented coefficients to satisfy the GMGPolar DensityProfileCoefficients concept._ [More...](#detailed-description)

* `#include <gmg_polar_poisson_like_solver.hpp>`





































## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**PolarPoissonLikeCoefficients**](#function-polarpoissonlikecoefficients) (int nr, int ntheta) <br>_Build the class instance._  |
|  KOKKOS\_INLINE\_FUNCTION double | [**alpha**](#function-alpha) (int i\_r, int i\_theta) const<br>_The coefficient alpha in the Poisson-like equation._  |
|  KOKKOS\_INLINE\_FUNCTION double | [**beta**](#function-beta) (int i\_r, int i\_theta) const<br>_The coefficient beta in the Poisson-like equation._  |
|  void | [**update\_coefficients**](#function-update_coefficients) (KokkosConstView2D alpha, KokkosConstView2D beta) <br>_Rebuild the internal representations of α and β from grid values._  |


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
    int nr,
    int ntheta
) 
```




<hr>



### function alpha 

_The coefficient alpha in the Poisson-like equation._ 
```C++
inline KOKKOS_INLINE_FUNCTION double GMGPolarTools::PolarPoissonLikeCoefficients::alpha (
    int i_r,
    int i_theta
) const
```




<hr>



### function beta 

_The coefficient beta in the Poisson-like equation._ 
```C++
inline KOKKOS_INLINE_FUNCTION double GMGPolarTools::PolarPoissonLikeCoefficients::beta (
    int i_r,
    int i_theta
) const
```




<hr>



### function update\_coefficients 

_Rebuild the internal representations of α and β from grid values._ 
```C++
inline void GMGPolarTools::PolarPoissonLikeCoefficients::update_coefficients (
    KokkosConstView2D alpha,
    KokkosConstView2D beta
) 
```





**Parameters:**


* `alpha` Values of α at the grid interpolation points. 
* `beta` Values of β at the grid interpolation points. 




        

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

