

# Struct PolarSplineFEMPoissonLikeSolver::EvalDeriv2DType



[**ClassList**](annotated.md) **>** [**PolarSplineFEMPoissonLikeSolver**](classPolarSplineFEMPoissonLikeSolver.md) **>** [**EvalDeriv2DType**](structPolarSplineFEMPoissonLikeSolver_1_1EvalDeriv2DType.md)



_Object storing a value and a value of the derivatives in each direction of a 2D function._ 

* `#include <polarpoissonlikesolver.hpp>`





















## Public Attributes

| Type | Name |
| ---: | :--- |
|  [**DVector**](classTensor.md)&lt; [**R\_cov**](structR__cov.md), [**Theta\_cov**](structTheta__cov.md) &gt; | [**derivative**](#variable-derivative)  <br>_The gradient of the function_  _._ |
|  double | [**value**](#variable-value)  <br>_The value of the function_  _._ |












































## Public Attributes Documentation




### variable derivative 

_The gradient of the function_  _._
```C++
DVector<R_cov, Theta_cov> PolarSplineFEMPoissonLikeSolver< GridR, GridTheta, PolarBSplinesRTheta, SplineRThetaEvaluatorNullBound, IdxRangeFull >::EvalDeriv2DType::derivative;
```




<hr>



### variable value 

_The value of the function_  _._
```C++
double PolarSplineFEMPoissonLikeSolver< GridR, GridTheta, PolarBSplinesRTheta, SplineRThetaEvaluatorNullBound, IdxRangeFull >::EvalDeriv2DType::value;
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/pde_solvers/polarpoissonlikesolver.hpp`

