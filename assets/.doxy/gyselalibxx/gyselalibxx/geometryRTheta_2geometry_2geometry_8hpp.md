

# File geometry.hpp



[**FileList**](files.md) **>** [**geometry**](dir_718520565cc7a7cfd9ba0e7c9c4c6d52.md) **>** [**geometry.hpp**](geometryRTheta_2geometry_2geometry_8hpp.md)

[Go to the source code of this file](geometryRTheta_2geometry_2geometry_8hpp_source.md)



* `#include <ddc/ddc.hpp>`
* `#include <ddc/kernels/splines.hpp>`
* `#include "ddc_alias_inline_functions.hpp"`
* `#include "ddc_aliases.hpp"`
* `#include "ddc_helper.hpp"`
* `#include "polar_bsplines.hpp"`
* `#include "vector_field.hpp"`
* `#include "vector_field_mem.hpp"`
* `#include "vector_index_tools.hpp"`















## Classes

| Type | Name |
| ---: | :--- |
| struct | [**BSplinesR**](structBSplinesR.md) <br> |
| struct | [**BSplinesTheta**](structBSplinesTheta.md) <br> |
| struct | [**GridR**](structGridR.md) <br> |
| struct | [**GridTheta**](structGridTheta.md) <br> |
| struct | [**PolarBSplinesRTheta**](structPolarBSplinesRTheta.md) <br> |
| struct | [**R**](structR.md) <br>_Define non periodic real contravariant_ [_**R**_](structR.md) _dimension._ |
| struct | [**R\_cov**](structR__cov.md) <br>_Define non periodic real covariant_ [_**R**_](structR.md) _dimension._ |
| struct | [**Theta**](structTheta.md) <br>_Define periodic real contravariant_ [_**Theta**_](structTheta.md) _dimension._ |
| struct | [**Theta\_cov**](structTheta__cov.md) <br>_Define periodic real covariant_ [_**Theta**_](structTheta.md) _dimension._ |
| struct | [**Vr**](structVr.md) <br>_Define non periodic real_ [_**R**_](structR.md) _velocity dimension._ |
| struct | [**Vtheta**](structVtheta.md) <br>_Define periodic real_ [_**Theta**_](structTheta.md) _velocity dimension._ |
| struct | [**Vx**](structVx.md) <br>_Define non periodic real_ [_**X**_](structX.md) _velocity dimension._ |
| struct | [**Vy**](structVy.md) <br>_Define non periodic real_ [_**Y**_](structY.md) _velocity dimension._ |
| struct | [**X**](structX.md) <br>_Define non periodic real_ [_**X**_](structX.md) _dimension._ |
| struct | [**Y**](structY.md) <br>_Define non periodic real_ [_**Y**_](structY.md) _dimension._ |


## Public Types

| Type | Name |
| ---: | :--- |
| typedef ConstField&lt; ElementType, IdxRangeR &gt; | [**ConstFieldR**](#typedef-constfieldr)  <br> |
| typedef ConstField&lt; ElementType, IdxRangeRTheta &gt; | [**ConstFieldRTheta**](#typedef-constfieldrtheta)  <br> |
| typedef ConstField&lt; ElementType, IdxRangeTheta &gt; | [**ConstFieldTheta**](#typedef-constfieldtheta)  <br> |
| typedef DConstField&lt; IdxRangeBSRTheta &gt; | [**ConstSpline2D**](#typedef-constspline2d)  <br> |
| typedef [**VectorConstField**](classVectorField.md)&lt; double, IdxRangeBSRTheta, VectorIndexSet&lt; Dim1, Dim2 &gt; &gt; | [**ConstVectorSplineCoeffs2D**](#typedef-constvectorsplinecoeffs2d)  <br> |
| typedef Coord&lt; [**R**](structR.md) &gt; | [**CoordR**](#typedef-coordr)  <br> |
| typedef Coord&lt; [**R**](structR.md), [**Theta**](structTheta.md) &gt; | [**CoordRTheta**](#typedef-coordrtheta)  <br> |
| typedef Coord&lt; [**Theta**](structTheta.md) &gt; | [**CoordTheta**](#typedef-coordtheta)  <br> |
| typedef Coord&lt; [**Vr**](structVr.md) &gt; | [**CoordVr**](#typedef-coordvr)  <br> |
| typedef Coord&lt; [**Vtheta**](structVtheta.md) &gt; | [**CoordVtheta**](#typedef-coordvtheta)  <br> |
| typedef Coord&lt; [**Vx**](structVx.md) &gt; | [**CoordVx**](#typedef-coordvx)  <br> |
| typedef Coord&lt; [**Vy**](structVy.md) &gt; | [**CoordVy**](#typedef-coordvy)  <br> |
| typedef Coord&lt; [**X**](structX.md) &gt; | [**CoordX**](#typedef-coordx)  <br> |
| typedef Coord&lt; [**X**](structX.md), [**Y**](structY.md) &gt; | [**CoordXY**](#typedef-coordxy)  <br> |
| typedef Coord&lt; [**Y**](structY.md) &gt; | [**CoordY**](#typedef-coordy)  <br> |
| typedef ConstFieldR&lt; double &gt; | [**DConstFieldR**](#typedef-dconstfieldr)  <br> |
| typedef ConstFieldRTheta&lt; double &gt; | [**DConstFieldRTheta**](#typedef-dconstfieldrtheta)  <br> |
| typedef ConstFieldTheta&lt; double &gt; | [**DConstFieldTheta**](#typedef-dconstfieldtheta)  <br> |
| typedef [**VectorConstField**](classVectorField.md)&lt; double, IdxRangeRTheta, VectorIndexSet&lt; Dim1, Dim2 &gt; &gt; | [**DConstVectorFieldRTheta**](#typedef-dconstvectorfieldrtheta)  <br> |
| typedef FieldMemR&lt; double &gt; | [**DFieldMemR**](#typedef-dfieldmemr)  <br> |
| typedef FieldMemRTheta&lt; double &gt; | [**DFieldMemRTheta**](#typedef-dfieldmemrtheta)  <br> |
| typedef FieldMemTheta&lt; double &gt; | [**DFieldMemTheta**](#typedef-dfieldmemtheta)  <br> |
| typedef FieldR&lt; double &gt; | [**DFieldR**](#typedef-dfieldr)  <br> |
| typedef FieldRTheta&lt; double &gt; | [**DFieldRTheta**](#typedef-dfieldrtheta)  <br> |
| typedef FieldTheta&lt; double &gt; | [**DFieldTheta**](#typedef-dfieldtheta)  <br> |
| typedef [**VectorFieldMem**](classVectorFieldMem.md)&lt; double, IdxRangeRTheta, VectorIndexSet&lt; Dim1, Dim2 &gt; &gt; | [**DVectorFieldMemRTheta**](#typedef-dvectorfieldmemrtheta)  <br> |
| typedef [**VectorField**](classVectorField.md)&lt; double, IdxRangeRTheta, VectorIndexSet&lt; Dim1, Dim2 &gt; &gt; | [**DVectorFieldRTheta**](#typedef-dvectorfieldrtheta)  <br> |
| typedef FieldMem&lt; ElementType, IdxRangeR &gt; | [**FieldMemR**](#typedef-fieldmemr)  <br> |
| typedef FieldMem&lt; ElementType, IdxRangeRTheta &gt; | [**FieldMemRTheta**](#typedef-fieldmemrtheta)  <br> |
| typedef FieldMem&lt; ElementType, IdxRangeTheta &gt; | [**FieldMemTheta**](#typedef-fieldmemtheta)  <br> |
| typedef Field&lt; ElementType, IdxRangeR &gt; | [**FieldR**](#typedef-fieldr)  <br> |
| typedef Field&lt; ElementType, IdxRangeRTheta &gt; | [**FieldRTheta**](#typedef-fieldrtheta)  <br> |
| typedef Field&lt; ElementType, IdxRangeTheta &gt; | [**FieldTheta**](#typedef-fieldtheta)  <br> |
| typedef Idx&lt; [**PolarBSplinesRTheta**](structPolarBSplinesRTheta.md) &gt; | [**IdxPolarBspl**](#typedef-idxpolarbspl)  <br>_Type of the index of an element of polar B-splines._  |
| typedef Idx&lt; [**GridR**](structGridR.md) &gt; | [**IdxR**](#typedef-idxr)  <br> |
| typedef Idx&lt; [**GridR**](structGridR.md), [**GridTheta**](structGridTheta.md) &gt; | [**IdxRTheta**](#typedef-idxrtheta)  <br> |
| typedef IdxRange&lt; [**PolarBSplinesRTheta**](structPolarBSplinesRTheta.md) &gt; | [**IdxRangeBSPolar**](#typedef-idxrangebspolar)  <br> |
| typedef IdxRange&lt; [**BSplinesR**](structBSplinesR.md) &gt; | [**IdxRangeBSR**](#typedef-idxrangebsr)  <br> |
| typedef IdxRange&lt; [**BSplinesR**](structBSplinesR.md), [**BSplinesTheta**](structBSplinesTheta.md) &gt; | [**IdxRangeBSRTheta**](#typedef-idxrangebsrtheta)  <br> |
| typedef IdxRange&lt; [**BSplinesTheta**](structBSplinesTheta.md) &gt; | [**IdxRangeBSTheta**](#typedef-idxrangebstheta)  <br> |
| typedef IdxRange&lt; [**GridR**](structGridR.md) &gt; | [**IdxRangeR**](#typedef-idxranger)  <br> |
| typedef IdxRange&lt; [**GridR**](structGridR.md), [**GridTheta**](structGridTheta.md) &gt; | [**IdxRangeRTheta**](#typedef-idxrangertheta)  <br> |
| typedef IdxRange&lt; [**GridTheta**](structGridTheta.md) &gt; | [**IdxRangeTheta**](#typedef-idxrangetheta)  <br> |
| typedef IdxStep&lt; [**GridR**](structGridR.md) &gt; | [**IdxStepR**](#typedef-idxstepr)  <br> |
| typedef IdxStep&lt; [**GridR**](structGridR.md), [**GridTheta**](structGridTheta.md) &gt; | [**IdxStepRTheta**](#typedef-idxsteprtheta)  <br> |
| typedef IdxStep&lt; [**GridTheta**](structGridTheta.md) &gt; | [**IdxStepTheta**](#typedef-idxsteptheta)  <br> |
| typedef Idx&lt; [**GridTheta**](structGridTheta.md) &gt; | [**IdxTheta**](#typedef-idxtheta)  <br> |
| typedef [**PolarSplineMem**](structPolarSplineMem.md)&lt; [**PolarBSplinesRTheta**](structPolarBSplinesRTheta.md) &gt; | [**PolarSplineMemRTheta**](#typedef-polarsplinememrtheta)  <br>_Tag the polar B-splines decomposition of a function._  |
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




### typedef ConstFieldR 

```C++
using ConstFieldR =  ConstField<ElementType, IdxRangeR>;
```




<hr>



### typedef ConstFieldRTheta 

```C++
using ConstFieldRTheta =  ConstField<ElementType, IdxRangeRTheta>;
```




<hr>



### typedef ConstFieldTheta 

```C++
using ConstFieldTheta =  ConstField<ElementType, IdxRangeTheta>;
```




<hr>



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



### typedef CoordR 

```C++
using CoordR =  Coord<R>;
```




<hr>



### typedef CoordRTheta 

```C++
using CoordRTheta =  Coord<R, Theta>;
```




<hr>



### typedef CoordTheta 

```C++
using CoordTheta =  Coord<Theta>;
```




<hr>



### typedef CoordVr 

```C++
using CoordVr =  Coord<Vr>;
```




<hr>



### typedef CoordVtheta 

```C++
using CoordVtheta =  Coord<Vtheta>;
```




<hr>



### typedef CoordVx 

```C++
using CoordVx =  Coord<Vx>;
```




<hr>



### typedef CoordVy 

```C++
using CoordVy =  Coord<Vy>;
```




<hr>



### typedef CoordX 

```C++
using CoordX =  Coord<X>;
```




<hr>



### typedef CoordXY 

```C++
using CoordXY =  Coord<X, Y>;
```




<hr>



### typedef CoordY 

```C++
using CoordY =  Coord<Y>;
```




<hr>



### typedef DConstFieldR 

```C++
using DConstFieldR =  ConstFieldR<double>;
```




<hr>



### typedef DConstFieldRTheta 

```C++
using DConstFieldRTheta =  ConstFieldRTheta<double>;
```




<hr>



### typedef DConstFieldTheta 

```C++
using DConstFieldTheta =  ConstFieldTheta<double>;
```




<hr>



### typedef DConstVectorFieldRTheta 

```C++
using DConstVectorFieldRTheta =  VectorConstField<double, IdxRangeRTheta, VectorIndexSet<Dim1, Dim2> >;
```




<hr>



### typedef DFieldMemR 

```C++
using DFieldMemR =  FieldMemR<double>;
```




<hr>



### typedef DFieldMemRTheta 

```C++
using DFieldMemRTheta =  FieldMemRTheta<double>;
```




<hr>



### typedef DFieldMemTheta 

```C++
using DFieldMemTheta =  FieldMemTheta<double>;
```




<hr>



### typedef DFieldR 

```C++
using DFieldR =  FieldR<double>;
```




<hr>



### typedef DFieldRTheta 

```C++
using DFieldRTheta =  FieldRTheta<double>;
```




<hr>



### typedef DFieldTheta 

```C++
using DFieldTheta =  FieldTheta<double>;
```




<hr>



### typedef DVectorFieldMemRTheta 

```C++
using DVectorFieldMemRTheta =  VectorFieldMem<double, IdxRangeRTheta, VectorIndexSet<Dim1, Dim2> >;
```




<hr>



### typedef DVectorFieldRTheta 

```C++
using DVectorFieldRTheta =  VectorField<double, IdxRangeRTheta, VectorIndexSet<Dim1, Dim2> >;
```




<hr>



### typedef FieldMemR 

```C++
using FieldMemR =  FieldMem<ElementType, IdxRangeR>;
```




<hr>



### typedef FieldMemRTheta 

```C++
using FieldMemRTheta =  FieldMem<ElementType, IdxRangeRTheta>;
```




<hr>



### typedef FieldMemTheta 

```C++
using FieldMemTheta =  FieldMem<ElementType, IdxRangeTheta>;
```




<hr>



### typedef FieldR 

```C++
using FieldR =  Field<ElementType, IdxRangeR>;
```




<hr>



### typedef FieldRTheta 

```C++
using FieldRTheta =  Field<ElementType, IdxRangeRTheta>;
```




<hr>



### typedef FieldTheta 

```C++
using FieldTheta =  Field<ElementType, IdxRangeTheta>;
```




<hr>



### typedef IdxPolarBspl 

_Type of the index of an element of polar B-splines._ 
```C++
using IdxPolarBspl =  Idx<PolarBSplinesRTheta>;
```




<hr>



### typedef IdxR 

```C++
using IdxR =  Idx<GridR>;
```




<hr>



### typedef IdxRTheta 

```C++
using IdxRTheta =  Idx<GridR, GridTheta>;
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



### typedef IdxRangeR 

```C++
using IdxRangeR =  IdxRange<GridR>;
```




<hr>



### typedef IdxRangeRTheta 

```C++
using IdxRangeRTheta =  IdxRange<GridR, GridTheta>;
```




<hr>



### typedef IdxRangeTheta 

```C++
using IdxRangeTheta =  IdxRange<GridTheta>;
```




<hr>



### typedef IdxStepR 

```C++
using IdxStepR =  IdxStep<GridR>;
```




<hr>



### typedef IdxStepRTheta 

```C++
using IdxStepRTheta =  IdxStep<GridR, GridTheta>;
```




<hr>



### typedef IdxStepTheta 

```C++
using IdxStepTheta =  IdxStep<GridTheta>;
```




<hr>



### typedef IdxTheta 

```C++
using IdxTheta =  Idx<GridTheta>;
```




<hr>



### typedef PolarSplineMemRTheta 

_Tag the polar B-splines decomposition of a function._ 
```C++
using PolarSplineMemRTheta =  PolarSplineMem<PolarBSplinesRTheta>;
```



Store the polar B-splines coefficients of the function. 


        

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
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/geometryRTheta/geometry/geometry.hpp`

