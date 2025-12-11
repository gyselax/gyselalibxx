

# File spline\_definitions.hpp



[**FileList**](files.md) **>** [**geometry**](dir_718520565cc7a7cfd9ba0e7c9c4c6d52.md) **>** [**spline\_definitions.hpp**](spline__definitions_8hpp.md)

[Go to the source code of this file](spline__definitions_8hpp_source.md)



* `#include <ddc/kernels/splines.hpp>`
* `#include "geometry_r_theta.hpp"`















## Classes

| Type | Name |
| ---: | :--- |
| struct | [**BSplinesR**](structBSplinesR.md) <br> |
| struct | [**BSplinesTheta**](structBSplinesTheta.md) <br> |
| struct | [**PolarBSplinesRTheta**](structPolarBSplinesRTheta.md) <br> |


## Public Types

| Type | Name |
| ---: | :--- |
| typedef DConstField&lt; IdxRangeBSRTheta &gt; | [**ConstSpline2D**](#typedef-constspline2d)  <br> |
| typedef [**VectorConstField**](classVectorField.md)&lt; double, IdxRangeBSRTheta, VectorIndexSet&lt; Dim1, Dim2 &gt; &gt; | [**ConstVectorSplineCoeffs2D**](#typedef-constvectorsplinecoeffs2d)  <br> |
| typedef Idx&lt; [**PolarBSplinesRTheta**](structPolarBSplinesRTheta.md) &gt; | [**IdxPolarBspl**](#typedef-idxpolarbspl)  <br>_Type of the index of an element of polar B-splines._  |
| typedef IdxRange&lt; [**PolarBSplinesRTheta**](structPolarBSplinesRTheta.md) &gt; | [**IdxRangeBSPolar**](#typedef-idxrangebspolar)  <br> |
| typedef IdxRange&lt; [**BSplinesR**](structBSplinesR.md) &gt; | [**IdxRangeBSR**](#typedef-idxrangebsr)  <br> |
| typedef IdxRange&lt; [**BSplinesR**](structBSplinesR.md), [**BSplinesTheta**](structBSplinesTheta.md) &gt; | [**IdxRangeBSRTheta**](#typedef-idxrangebsrtheta)  <br> |
| typedef IdxRange&lt; [**BSplinesTheta**](structBSplinesTheta.md) &gt; | [**IdxRangeBSTheta**](#typedef-idxrangebstheta)  <br> |
| typedef DFieldMem&lt; IdxRange&lt; [**PolarBSplinesRTheta**](structPolarBSplinesRTheta.md) &gt; &gt; | [**PolarSplineMemRTheta**](#typedef-polarsplinememrtheta)  <br>_Tag the polar B-splines decomposition of a function._  |
| typedef DField&lt; IdxRange&lt; [**PolarBSplinesRTheta**](structPolarBSplinesRTheta.md) &gt; &gt; | [**PolarSplineRTheta**](#typedef-polarsplinertheta)  <br> |
| typedef DField&lt; IdxRangeBSRTheta &gt; | [**Spline2D**](#typedef-spline2d)  <br> |
| typedef DFieldMem&lt; IdxRangeBSRTheta &gt; | [**Spline2DMem**](#typedef-spline2dmem)  <br> |
| typedef ddc::GrevilleInterpolationPoints&lt; [**BSplinesR**](structBSplinesR.md), SplineRBoundary, SplineRBoundary &gt; | [**SplineInterpPointsR**](#typedef-splineinterppointsr)  <br> |
| typedef ddc::GrevilleInterpolationPoints&lt; [**BSplinesTheta**](structBSplinesTheta.md), SplineThetaBoundary, SplineThetaBoundary &gt; | [**SplineInterpPointsTheta**](#typedef-splineinterppointstheta)  <br> |
| typedef ddc::SplineBuilder2D&lt; Kokkos::DefaultExecutionSpace, typename Kokkos::DefaultExecutionSpace::memory\_space, [**BSplinesR**](structBSplinesR.md), [**BSplinesTheta**](structBSplinesTheta.md), [**GridR**](structGridR.md), [**GridTheta**](structGridTheta.md), SplineRBoundary, SplineRBoundary, SplineThetaBoundary, SplineThetaBoundary, ddc::SplineSolver::LAPACK &gt; | [**SplineRThetaBuilder**](#typedef-splinerthetabuilder)  <br> |
| typedef ddc::SplineBuilder2D&lt; Kokkos::DefaultHostExecutionSpace, Kokkos::HostSpace, [**BSplinesR**](structBSplinesR.md), [**BSplinesTheta**](structBSplinesTheta.md), [**GridR**](structGridR.md), [**GridTheta**](structGridTheta.md), SplineRBoundary, SplineRBoundary, SplineThetaBoundary, SplineThetaBoundary, ddc::SplineSolver::LAPACK &gt; | [**SplineRThetaBuilder\_host**](#typedef-splinerthetabuilder_host)  <br> |
| typedef ddc::SplineEvaluator2D&lt; Kokkos::DefaultExecutionSpace, typename Kokkos::DefaultExecutionSpace::memory\_space, [**BSplinesR**](structBSplinesR.md), [**BSplinesTheta**](structBSplinesTheta.md), [**GridR**](structGridR.md), [**GridTheta**](structGridTheta.md), ddc::ConstantExtrapolationRule&lt; [**R**](structR.md), [**Theta**](structTheta.md) &gt;, ddc::ConstantExtrapolationRule&lt; [**R**](structR.md), [**Theta**](structTheta.md) &gt;, ddc::PeriodicExtrapolationRule&lt; [**Theta**](structTheta.md) &gt;, ddc::PeriodicExtrapolationRule&lt; [**Theta**](structTheta.md) &gt; &gt; | [**SplineRThetaEvaluatorConstBound**](#typedef-splinerthetaevaluatorconstbound)  <br> |
| typedef ddc::SplineEvaluator2D&lt; Kokkos::DefaultHostExecutionSpace, Kokkos::HostSpace, [**BSplinesR**](structBSplinesR.md), [**BSplinesTheta**](structBSplinesTheta.md), [**GridR**](structGridR.md), [**GridTheta**](structGridTheta.md), ddc::ConstantExtrapolationRule&lt; [**R**](structR.md), [**Theta**](structTheta.md) &gt;, ddc::ConstantExtrapolationRule&lt; [**R**](structR.md), [**Theta**](structTheta.md) &gt;, ddc::PeriodicExtrapolationRule&lt; [**Theta**](structTheta.md) &gt;, ddc::PeriodicExtrapolationRule&lt; [**Theta**](structTheta.md) &gt; &gt; | [**SplineRThetaEvaluatorConstBound\_host**](#typedef-splinerthetaevaluatorconstbound_host)  <br> |
| typedef ddc::SplineEvaluator2D&lt; Kokkos::DefaultExecutionSpace, typename Kokkos::DefaultExecutionSpace::memory\_space, [**BSplinesR**](structBSplinesR.md), [**BSplinesTheta**](structBSplinesTheta.md), [**GridR**](structGridR.md), [**GridTheta**](structGridTheta.md), ddc::NullExtrapolationRule, ddc::NullExtrapolationRule, ddc::PeriodicExtrapolationRule&lt; [**Theta**](structTheta.md) &gt;, ddc::PeriodicExtrapolationRule&lt; [**Theta**](structTheta.md) &gt; &gt; | [**SplineRThetaEvaluatorNullBound**](#typedef-splinerthetaevaluatornullbound)  <br> |
| typedef ddc::SplineEvaluator2D&lt; Kokkos::DefaultHostExecutionSpace, Kokkos::HostSpace, [**BSplinesR**](structBSplinesR.md), [**BSplinesTheta**](structBSplinesTheta.md), [**GridR**](structGridR.md), [**GridTheta**](structGridTheta.md), ddc::NullExtrapolationRule, ddc::NullExtrapolationRule, ddc::PeriodicExtrapolationRule&lt; [**Theta**](structTheta.md) &gt;, ddc::PeriodicExtrapolationRule&lt; [**Theta**](structTheta.md) &gt; &gt; | [**SplineRThetaEvaluatorNullBound\_host**](#typedef-splinerthetaevaluatornullbound_host)  <br> |
| typedef [**VectorField**](classVectorField.md)&lt; double, IdxRangeBSRTheta, VectorIndexSet&lt; Dim1, Dim2 &gt; &gt; | [**VectorSplineCoeffs2D**](#typedef-vectorsplinecoeffs2d)  <br> |
| typedef [**VectorFieldMem**](classVectorFieldMem.md)&lt; double, IdxRangeBSRTheta, VectorIndexSet&lt; Dim1, Dim2 &gt; &gt; | [**VectorSplineCoeffsMem2D**](#typedef-vectorsplinecoeffsmem2d)  <br> |




## Public Attributes

| Type | Name |
| ---: | :--- |
|  int constexpr | [**BSDegreeR**](#variable-bsdegreer)   = `3`<br> |
|  int constexpr | [**BSDegreeTheta**](#variable-bsdegreetheta)   = `3`<br> |
|  bool constexpr | [**BsplineOnUniformCellsR**](#variable-bsplineonuniformcellsr)   = `false`<br> |
|  bool constexpr | [**BsplineOnUniformCellsTheta**](#variable-bsplineonuniformcellstheta)   = `false`<br> |
|  ddc::BoundCond constexpr | [**SplineRBoundary**](#variable-splinerboundary)   = `ddc::BoundCond::GREVILLE`<br> |
|  ddc::BoundCond constexpr | [**SplineThetaBoundary**](#variable-splinethetaboundary)   = `ddc::BoundCond::PERIODIC`<br> |












































## Public Types Documentation




### typedef ConstSpline2D 

```C++
using ConstSpline2D =  DConstField<IdxRangeBSRTheta>;
```




<hr>



### typedef ConstVectorSplineCoeffs2D 

```C++
using ConstVectorSplineCoeffs2D =  VectorConstField<double, IdxRangeBSRTheta, VectorIndexSet<Dim1, Dim2> >;
```




<hr>



### typedef IdxPolarBspl 

_Type of the index of an element of polar B-splines._ 
```C++
using IdxPolarBspl =  Idx<PolarBSplinesRTheta>;
```




<hr>



### typedef IdxRangeBSPolar 

```C++
using IdxRangeBSPolar =  IdxRange<PolarBSplinesRTheta>;
```




<hr>



### typedef IdxRangeBSR 

```C++
using IdxRangeBSR =  IdxRange<BSplinesR>;
```




<hr>



### typedef IdxRangeBSRTheta 

```C++
using IdxRangeBSRTheta =  IdxRange<BSplinesR, BSplinesTheta>;
```




<hr>



### typedef IdxRangeBSTheta 

```C++
using IdxRangeBSTheta =  IdxRange<BSplinesTheta>;
```




<hr>



### typedef PolarSplineMemRTheta 

_Tag the polar B-splines decomposition of a function._ 
```C++
using PolarSplineMemRTheta =  DFieldMem<IdxRange<PolarBSplinesRTheta> >;
```



Store the polar B-splines coefficients of the function. 


        

<hr>



### typedef PolarSplineRTheta 

```C++
using PolarSplineRTheta =  DField<IdxRange<PolarBSplinesRTheta> >;
```




<hr>



### typedef Spline2D 

```C++
using Spline2D =  DField<IdxRangeBSRTheta>;
```




<hr>



### typedef Spline2DMem 

```C++
using Spline2DMem =  DFieldMem<IdxRangeBSRTheta>;
```




<hr>



### typedef SplineInterpPointsR 

```C++
using SplineInterpPointsR =  ddc::GrevilleInterpolationPoints<BSplinesR, SplineRBoundary, SplineRBoundary>;
```




<hr>



### typedef SplineInterpPointsTheta 

```C++
using SplineInterpPointsTheta =  ddc::GrevilleInterpolationPoints<BSplinesTheta, SplineThetaBoundary, SplineThetaBoundary>;
```




<hr>



### typedef SplineRThetaBuilder 

```C++
using SplineRThetaBuilder =  ddc::SplineBuilder2D< Kokkos::DefaultExecutionSpace, typename Kokkos::DefaultExecutionSpace::memory_space, BSplinesR, BSplinesTheta, GridR, GridTheta, SplineRBoundary, SplineRBoundary, SplineThetaBoundary, SplineThetaBoundary, ddc::SplineSolver::LAPACK>;
```




<hr>



### typedef SplineRThetaBuilder\_host 

```C++
using SplineRThetaBuilder_host =  ddc::SplineBuilder2D< Kokkos::DefaultHostExecutionSpace, Kokkos::HostSpace, BSplinesR, BSplinesTheta, GridR, GridTheta, SplineRBoundary, SplineRBoundary, SplineThetaBoundary, SplineThetaBoundary, ddc::SplineSolver::LAPACK>;
```




<hr>



### typedef SplineRThetaEvaluatorConstBound 

```C++
using SplineRThetaEvaluatorConstBound =  ddc::SplineEvaluator2D< Kokkos::DefaultExecutionSpace, typename Kokkos::DefaultExecutionSpace::memory_space, BSplinesR, BSplinesTheta, GridR, GridTheta, ddc::ConstantExtrapolationRule<R, Theta>, ddc::ConstantExtrapolationRule<R, Theta>, ddc::PeriodicExtrapolationRule<Theta>, ddc::PeriodicExtrapolationRule<Theta> >;
```




<hr>



### typedef SplineRThetaEvaluatorConstBound\_host 

```C++
using SplineRThetaEvaluatorConstBound_host =  ddc::SplineEvaluator2D< Kokkos::DefaultHostExecutionSpace, Kokkos::HostSpace, BSplinesR, BSplinesTheta, GridR, GridTheta, ddc::ConstantExtrapolationRule<R, Theta>, ddc::ConstantExtrapolationRule<R, Theta>, ddc::PeriodicExtrapolationRule<Theta>, ddc::PeriodicExtrapolationRule<Theta> >;
```




<hr>



### typedef SplineRThetaEvaluatorNullBound 

```C++
using SplineRThetaEvaluatorNullBound =  ddc::SplineEvaluator2D< Kokkos::DefaultExecutionSpace, typename Kokkos::DefaultExecutionSpace::memory_space, BSplinesR, BSplinesTheta, GridR, GridTheta, ddc::NullExtrapolationRule, ddc::NullExtrapolationRule, ddc::PeriodicExtrapolationRule<Theta>, ddc::PeriodicExtrapolationRule<Theta> >;
```




<hr>



### typedef SplineRThetaEvaluatorNullBound\_host 

```C++
using SplineRThetaEvaluatorNullBound_host =  ddc::SplineEvaluator2D< Kokkos::DefaultHostExecutionSpace, Kokkos::HostSpace, BSplinesR, BSplinesTheta, GridR, GridTheta, ddc::NullExtrapolationRule, ddc::NullExtrapolationRule, ddc::PeriodicExtrapolationRule<Theta>, ddc::PeriodicExtrapolationRule<Theta> >;
```




<hr>



### typedef VectorSplineCoeffs2D 

```C++
using VectorSplineCoeffs2D =  VectorField<double, IdxRangeBSRTheta, VectorIndexSet<Dim1, Dim2> >;
```




<hr>



### typedef VectorSplineCoeffsMem2D 

```C++
using VectorSplineCoeffsMem2D =  VectorFieldMem<double, IdxRangeBSRTheta, VectorIndexSet<Dim1, Dim2> >;
```




<hr>
## Public Attributes Documentation




### variable BSDegreeR 

```C++
int constexpr BSDegreeR;
```




<hr>



### variable BSDegreeTheta 

```C++
int constexpr BSDegreeTheta;
```




<hr>



### variable BsplineOnUniformCellsR 

```C++
bool constexpr BsplineOnUniformCellsR;
```




<hr>



### variable BsplineOnUniformCellsTheta 

```C++
bool constexpr BsplineOnUniformCellsTheta;
```




<hr>



### variable SplineRBoundary 

```C++
ddc::BoundCond constexpr SplineRBoundary;
```




<hr>



### variable SplineThetaBoundary 

```C++
ddc::BoundCond constexpr SplineThetaBoundary;
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/geometryRTheta/geometry/spline_definitions.hpp`

