

# File geometry.hpp



[**FileList**](files.md) **>** [**geometry**](dir_807bff9d645a62665fbe11aeed095652.md) **>** [**geometry.hpp**](geometryVparMu_2geometry_2geometry_8hpp.md)

[Go to the source code of this file](geometryVparMu_2geometry_2geometry_8hpp_source.md)



* `#include <ddc/ddc.hpp>`
* `#include <ddc/kernels/splines.hpp>`
* `#include "ddc_alias_inline_functions.hpp"`
* `#include "ddc_aliases.hpp"`
* `#include "ddc_helper.hpp"`
* `#include "species_info.hpp"`















## Classes

| Type | Name |
| ---: | :--- |
| struct | [**BSplinesMu**](structBSplinesMu.md) <br> |
| struct | [**BSplinesVpar**](structBSplinesVpar.md) <br> |
| struct | [**GridMu**](structGridMu.md) <br> |
| struct | [**GridVpar**](structGridVpar.md) <br> |
| struct | [**Mu**](structMu.md) <br>_Define non periodic magnetic momentum_  _._ |
| struct | [**Vpar**](structVpar.md) <br>_Define non periodic parallel velocity_  _._ |


## Public Types

| Type | Name |
| ---: | :--- |
| typedef ConstField&lt; ElementType, IdxRangeMu &gt; | [**ConstFieldMu**](#typedef-constfieldmu)  <br> |
| typedef ConstField&lt; ElementType, IdxRangeSpVpar &gt; | [**ConstFieldSpVpar**](#typedef-constfieldspvpar)  <br> |
| typedef ConstField&lt; ElementType, IdxRangeSpVparMu &gt; | [**ConstFieldSpVparMu**](#typedef-constfieldspvparmu)  <br> |
| typedef ConstField&lt; ElementType, IdxRangeVpar &gt; | [**ConstFieldVpar**](#typedef-constfieldvpar)  <br> |
| typedef ConstField&lt; ElementType, IdxRangeVparMu &gt; | [**ConstFieldVparMu**](#typedef-constfieldvparmu)  <br> |
| typedef Coord&lt; [**Mu**](structMu.md) &gt; | [**CoordMu**](#typedef-coordmu)  <br> |
| typedef Coord&lt; [**Vpar**](structVpar.md) &gt; | [**CoordVpar**](#typedef-coordvpar)  <br> |
| typedef ConstFieldMu&lt; double &gt; | [**DConstFieldMu**](#typedef-dconstfieldmu)  <br> |
| typedef ConstFieldSpVpar&lt; double &gt; | [**DConstFieldSpVpar**](#typedef-dconstfieldspvpar)  <br> |
| typedef ConstFieldSpVparMu&lt; double &gt; | [**DConstFieldSpVparMu**](#typedef-dconstfieldspvparmu)  <br> |
| typedef ConstFieldVpar&lt; double &gt; | [**DConstFieldVpar**](#typedef-dconstfieldvpar)  <br> |
| typedef ConstFieldVparMu&lt; double &gt; | [**DConstFieldVparMu**](#typedef-dconstfieldvparmu)  <br> |
| typedef FieldMemMu&lt; double &gt; | [**DFieldMemMu**](#typedef-dfieldmemmu)  <br> |
| typedef FieldMemSpVpar&lt; double &gt; | [**DFieldMemSpVpar**](#typedef-dfieldmemspvpar)  <br> |
| typedef FieldMemSpVparMu&lt; double &gt; | [**DFieldMemSpVparMu**](#typedef-dfieldmemspvparmu)  <br> |
| typedef FieldMemVpar&lt; double &gt; | [**DFieldMemVpar**](#typedef-dfieldmemvpar)  <br> |
| typedef FieldMemVparMu&lt; double &gt; | [**DFieldMemVparMu**](#typedef-dfieldmemvparmu)  <br> |
| typedef FieldMu&lt; double &gt; | [**DFieldMu**](#typedef-dfieldmu)  <br> |
| typedef FieldSpVpar&lt; double &gt; | [**DFieldSpVpar**](#typedef-dfieldspvpar)  <br> |
| typedef FieldSpVparMu&lt; double &gt; | [**DFieldSpVparMu**](#typedef-dfieldspvparmu)  <br> |
| typedef FieldVpar&lt; double &gt; | [**DFieldVpar**](#typedef-dfieldvpar)  <br> |
| typedef FieldVparMu&lt; double &gt; | [**DFieldVparMu**](#typedef-dfieldvparmu)  <br> |
| typedef FieldMem&lt; ElementType, IdxRangeMu &gt; | [**FieldMemMu**](#typedef-fieldmemmu)  <br> |
| typedef FieldMem&lt; ElementType, IdxRangeSpVpar &gt; | [**FieldMemSpVpar**](#typedef-fieldmemspvpar)  <br> |
| typedef FieldMem&lt; ElementType, IdxRangeSpVparMu &gt; | [**FieldMemSpVparMu**](#typedef-fieldmemspvparmu)  <br> |
| typedef FieldMem&lt; ElementType, IdxRangeVpar &gt; | [**FieldMemVpar**](#typedef-fieldmemvpar)  <br> |
| typedef FieldMem&lt; ElementType, IdxRangeVparMu &gt; | [**FieldMemVparMu**](#typedef-fieldmemvparmu)  <br> |
| typedef Field&lt; ElementType, IdxRangeMu &gt; | [**FieldMu**](#typedef-fieldmu)  <br> |
| typedef Field&lt; ElementType, IdxRangeSpVpar &gt; | [**FieldSpVpar**](#typedef-fieldspvpar)  <br> |
| typedef Field&lt; ElementType, IdxRangeSpVparMu &gt; | [**FieldSpVparMu**](#typedef-fieldspvparmu)  <br> |
| typedef Field&lt; ElementType, IdxRangeVpar &gt; | [**FieldVpar**](#typedef-fieldvpar)  <br> |
| typedef Field&lt; ElementType, IdxRangeVparMu &gt; | [**FieldVparMu**](#typedef-fieldvparmu)  <br> |
| typedef Idx&lt; [**GridMu**](structGridMu.md) &gt; | [**IdxMu**](#typedef-idxmu)  <br> |
| typedef IdxRange&lt; [**GridMu**](structGridMu.md) &gt; | [**IdxRangeMu**](#typedef-idxrangemu)  <br> |
| typedef IdxRange&lt; [**Species**](structSpecies.md), [**GridVpar**](structGridVpar.md) &gt; | [**IdxRangeSpVpar**](#typedef-idxrangespvpar)  <br> |
| typedef IdxRange&lt; [**Species**](structSpecies.md), [**GridVpar**](structGridVpar.md), [**GridMu**](structGridMu.md) &gt; | [**IdxRangeSpVparMu**](#typedef-idxrangespvparmu)  <br> |
| typedef IdxRange&lt; [**GridVpar**](structGridVpar.md) &gt; | [**IdxRangeVpar**](#typedef-idxrangevpar)  <br> |
| typedef IdxRange&lt; [**GridVpar**](structGridVpar.md), [**GridMu**](structGridMu.md) &gt; | [**IdxRangeVparMu**](#typedef-idxrangevparmu)  <br> |
| typedef Idx&lt; [**Species**](structSpecies.md), [**GridVpar**](structGridVpar.md) &gt; | [**IdxSpVpar**](#typedef-idxspvpar)  <br> |
| typedef Idx&lt; [**Species**](structSpecies.md), [**GridVpar**](structGridVpar.md), [**GridMu**](structGridMu.md) &gt; | [**IdxSpVparMu**](#typedef-idxspvparmu)  <br> |
| typedef IdxStep&lt; [**GridMu**](structGridMu.md) &gt; | [**IdxStepMu**](#typedef-idxstepmu)  <br> |
| typedef IdxStep&lt; [**Species**](structSpecies.md), [**GridVpar**](structGridVpar.md), [**GridMu**](structGridMu.md) &gt; | [**IdxStepSpVparMu**](#typedef-idxstepspvparmu)  <br> |
| typedef IdxStep&lt; [**GridVpar**](structGridVpar.md) &gt; | [**IdxStepVpar**](#typedef-idxstepvpar)  <br> |
| typedef IdxStep&lt; [**GridVpar**](structGridVpar.md), [**GridMu**](structGridMu.md) &gt; | [**IdxStepVparMu**](#typedef-idxstepvparmu)  <br> |
| typedef Idx&lt; [**GridVpar**](structGridVpar.md) &gt; | [**IdxVpar**](#typedef-idxvpar)  <br> |
| typedef Idx&lt; [**GridVpar**](structGridVpar.md), [**GridMu**](structGridMu.md) &gt; | [**IdxVparMu**](#typedef-idxvparmu)  <br> |
| typedef ddc::GrevilleInterpolationPoints&lt; [**BSplinesMu**](structBSplinesMu.md), SplineMuBoundary, SplineMuBoundary &gt; | [**SplineInterpPointsMu**](#typedef-splineinterppointsmu)  <br> |
| typedef ddc::GrevilleInterpolationPoints&lt; [**BSplinesVpar**](structBSplinesVpar.md), SplineVparBoundary, SplineVparBoundary &gt; | [**SplineInterpPointsVpar**](#typedef-splineinterppointsvpar)  <br> |
| typedef ddc::SplineBuilder&lt; Kokkos::DefaultExecutionSpace, Kokkos::DefaultExecutionSpace::memory\_space, [**BSplinesMu**](structBSplinesMu.md), [**GridMu**](structGridMu.md), SplineMuBoundary, SplineMuBoundary, ddc::SplineSolver::LAPACK &gt; | [**SplineMuBuilder**](#typedef-splinemubuilder)  <br> |
| typedef ddc::SplineEvaluator&lt; Kokkos::DefaultExecutionSpace, Kokkos::DefaultExecutionSpace::memory\_space, [**BSplinesMu**](structBSplinesMu.md), [**GridMu**](structGridMu.md), ddc::ConstantExtrapolationRule&lt; [**Mu**](structMu.md) &gt;, ddc::ConstantExtrapolationRule&lt; [**Mu**](structMu.md) &gt; &gt; | [**SplineMuEvaluator**](#typedef-splinemuevaluator)  <br> |
| typedef ddc::SplineBuilder&lt; Kokkos::DefaultExecutionSpace, Kokkos::DefaultExecutionSpace::memory\_space, [**BSplinesVpar**](structBSplinesVpar.md), [**GridVpar**](structGridVpar.md), SplineVparBoundary, SplineVparBoundary, ddc::SplineSolver::LAPACK &gt; | [**SplineVparBuilder**](#typedef-splinevparbuilder)  <br> |
| typedef ddc::SplineEvaluator&lt; Kokkos::DefaultExecutionSpace, Kokkos::DefaultExecutionSpace::memory\_space, [**BSplinesVpar**](structBSplinesVpar.md), [**GridVpar**](structGridVpar.md), ddc::ConstantExtrapolationRule&lt; [**Vpar**](structVpar.md) &gt;, ddc::ConstantExtrapolationRule&lt; [**Vpar**](structVpar.md) &gt; &gt; | [**SplineVparEvaluator**](#typedef-splinevparevaluator)  <br> |




## Public Attributes

| Type | Name |
| ---: | :--- |
|  int constexpr | [**BSDegreeMu**](#variable-bsdegreemu)   = `3`<br> |
|  int constexpr | [**BSDegreeVpar**](#variable-bsdegreevpar)   = `3`<br> |
|  bool constexpr | [**BsplineOnUniformCellsMu**](#variable-bsplineonuniformcellsmu)   = `true`<br> |
|  bool constexpr | [**BsplineOnUniformCellsVpar**](#variable-bsplineonuniformcellsvpar)   = `true`<br> |
|  ddc::BoundCond constexpr | [**SplineMuBoundary**](#variable-splinemuboundary)   = `ddc::BoundCond::HERMITE`<br> |
|  ddc::BoundCond constexpr | [**SplineVparBoundary**](#variable-splinevparboundary)   = `ddc::BoundCond::HERMITE`<br> |












































## Public Types Documentation




### typedef ConstFieldMu 

```C++
using ConstFieldMu =  ConstField<ElementType, IdxRangeMu>;
```




<hr>



### typedef ConstFieldSpVpar 

```C++
using ConstFieldSpVpar =  ConstField<ElementType, IdxRangeSpVpar>;
```




<hr>



### typedef ConstFieldSpVparMu 

```C++
using ConstFieldSpVparMu =  ConstField<ElementType, IdxRangeSpVparMu>;
```




<hr>



### typedef ConstFieldVpar 

```C++
using ConstFieldVpar =  ConstField<ElementType, IdxRangeVpar>;
```




<hr>



### typedef ConstFieldVparMu 

```C++
using ConstFieldVparMu =  ConstField<ElementType, IdxRangeVparMu>;
```




<hr>



### typedef CoordMu 

```C++
using CoordMu =  Coord<Mu>;
```




<hr>



### typedef CoordVpar 

```C++
using CoordVpar =  Coord<Vpar>;
```




<hr>



### typedef DConstFieldMu 

```C++
using DConstFieldMu =  ConstFieldMu<double>;
```




<hr>



### typedef DConstFieldSpVpar 

```C++
using DConstFieldSpVpar =  ConstFieldSpVpar<double>;
```




<hr>



### typedef DConstFieldSpVparMu 

```C++
using DConstFieldSpVparMu =  ConstFieldSpVparMu<double>;
```




<hr>



### typedef DConstFieldVpar 

```C++
using DConstFieldVpar =  ConstFieldVpar<double>;
```




<hr>



### typedef DConstFieldVparMu 

```C++
using DConstFieldVparMu =  ConstFieldVparMu<double>;
```




<hr>



### typedef DFieldMemMu 

```C++
using DFieldMemMu =  FieldMemMu<double>;
```




<hr>



### typedef DFieldMemSpVpar 

```C++
using DFieldMemSpVpar =  FieldMemSpVpar<double>;
```




<hr>



### typedef DFieldMemSpVparMu 

```C++
using DFieldMemSpVparMu =  FieldMemSpVparMu<double>;
```




<hr>



### typedef DFieldMemVpar 

```C++
using DFieldMemVpar =  FieldMemVpar<double>;
```




<hr>



### typedef DFieldMemVparMu 

```C++
using DFieldMemVparMu =  FieldMemVparMu<double>;
```




<hr>



### typedef DFieldMu 

```C++
using DFieldMu =  FieldMu<double>;
```




<hr>



### typedef DFieldSpVpar 

```C++
using DFieldSpVpar =  FieldSpVpar<double>;
```




<hr>



### typedef DFieldSpVparMu 

```C++
using DFieldSpVparMu =  FieldSpVparMu<double>;
```




<hr>



### typedef DFieldVpar 

```C++
using DFieldVpar =  FieldVpar<double>;
```




<hr>



### typedef DFieldVparMu 

```C++
using DFieldVparMu =  FieldVparMu<double>;
```




<hr>



### typedef FieldMemMu 

```C++
using FieldMemMu =  FieldMem<ElementType, IdxRangeMu>;
```




<hr>



### typedef FieldMemSpVpar 

```C++
using FieldMemSpVpar =  FieldMem<ElementType, IdxRangeSpVpar>;
```




<hr>



### typedef FieldMemSpVparMu 

```C++
using FieldMemSpVparMu =  FieldMem<ElementType, IdxRangeSpVparMu>;
```




<hr>



### typedef FieldMemVpar 

```C++
using FieldMemVpar =  FieldMem<ElementType, IdxRangeVpar>;
```




<hr>



### typedef FieldMemVparMu 

```C++
using FieldMemVparMu =  FieldMem<ElementType, IdxRangeVparMu>;
```




<hr>



### typedef FieldMu 

```C++
using FieldMu =  Field<ElementType, IdxRangeMu>;
```




<hr>



### typedef FieldSpVpar 

```C++
using FieldSpVpar =  Field<ElementType, IdxRangeSpVpar>;
```




<hr>



### typedef FieldSpVparMu 

```C++
using FieldSpVparMu =  Field<ElementType, IdxRangeSpVparMu>;
```




<hr>



### typedef FieldVpar 

```C++
using FieldVpar =  Field<ElementType, IdxRangeVpar>;
```




<hr>



### typedef FieldVparMu 

```C++
using FieldVparMu =  Field<ElementType, IdxRangeVparMu>;
```




<hr>



### typedef IdxMu 

```C++
using IdxMu =  Idx<GridMu>;
```




<hr>



### typedef IdxRangeMu 

```C++
using IdxRangeMu =  IdxRange<GridMu>;
```




<hr>



### typedef IdxRangeSpVpar 

```C++
using IdxRangeSpVpar =  IdxRange<Species, GridVpar>;
```




<hr>



### typedef IdxRangeSpVparMu 

```C++
using IdxRangeSpVparMu =  IdxRange<Species, GridVpar, GridMu>;
```




<hr>



### typedef IdxRangeVpar 

```C++
using IdxRangeVpar =  IdxRange<GridVpar>;
```




<hr>



### typedef IdxRangeVparMu 

```C++
using IdxRangeVparMu =  IdxRange<GridVpar, GridMu>;
```




<hr>



### typedef IdxSpVpar 

```C++
using IdxSpVpar =  Idx<Species, GridVpar>;
```




<hr>



### typedef IdxSpVparMu 

```C++
using IdxSpVparMu =  Idx<Species, GridVpar, GridMu>;
```




<hr>



### typedef IdxStepMu 

```C++
using IdxStepMu =  IdxStep<GridMu>;
```




<hr>



### typedef IdxStepSpVparMu 

```C++
using IdxStepSpVparMu =  IdxStep<Species, GridVpar, GridMu>;
```




<hr>



### typedef IdxStepVpar 

```C++
using IdxStepVpar =  IdxStep<GridVpar>;
```




<hr>



### typedef IdxStepVparMu 

```C++
using IdxStepVparMu =  IdxStep<GridVpar, GridMu>;
```




<hr>



### typedef IdxVpar 

```C++
using IdxVpar =  Idx<GridVpar>;
```




<hr>



### typedef IdxVparMu 

```C++
using IdxVparMu =  Idx<GridVpar, GridMu>;
```




<hr>



### typedef SplineInterpPointsMu 

```C++
using SplineInterpPointsMu =  ddc::GrevilleInterpolationPoints<BSplinesMu, SplineMuBoundary, SplineMuBoundary>;
```




<hr>



### typedef SplineInterpPointsVpar 

```C++
using SplineInterpPointsVpar =  ddc::GrevilleInterpolationPoints<BSplinesVpar, SplineVparBoundary, SplineVparBoundary>;
```




<hr>



### typedef SplineMuBuilder 

```C++
using SplineMuBuilder =  ddc::SplineBuilder< Kokkos::DefaultExecutionSpace, Kokkos::DefaultExecutionSpace::memory_space, BSplinesMu, GridMu, SplineMuBoundary, SplineMuBoundary, ddc::SplineSolver::LAPACK>;
```




<hr>



### typedef SplineMuEvaluator 

```C++
using SplineMuEvaluator =  ddc::SplineEvaluator< Kokkos::DefaultExecutionSpace, Kokkos::DefaultExecutionSpace::memory_space, BSplinesMu, GridMu, ddc::ConstantExtrapolationRule<Mu>, ddc::ConstantExtrapolationRule<Mu> >;
```




<hr>



### typedef SplineVparBuilder 

```C++
using SplineVparBuilder =  ddc::SplineBuilder< Kokkos::DefaultExecutionSpace, Kokkos::DefaultExecutionSpace::memory_space, BSplinesVpar, GridVpar, SplineVparBoundary, SplineVparBoundary, ddc::SplineSolver::LAPACK>;
```




<hr>



### typedef SplineVparEvaluator 

```C++
using SplineVparEvaluator =  ddc::SplineEvaluator< Kokkos::DefaultExecutionSpace, Kokkos::DefaultExecutionSpace::memory_space, BSplinesVpar, GridVpar, ddc::ConstantExtrapolationRule<Vpar>, ddc::ConstantExtrapolationRule<Vpar> >;
```




<hr>
## Public Attributes Documentation




### variable BSDegreeMu 

```C++
int constexpr BSDegreeMu;
```




<hr>



### variable BSDegreeVpar 

```C++
int constexpr BSDegreeVpar;
```




<hr>



### variable BsplineOnUniformCellsMu 

```C++
bool constexpr BsplineOnUniformCellsMu;
```




<hr>



### variable BsplineOnUniformCellsVpar 

```C++
bool constexpr BsplineOnUniformCellsVpar;
```




<hr>



### variable SplineMuBoundary 

```C++
ddc::BoundCond constexpr SplineMuBoundary;
```




<hr>



### variable SplineVparBoundary 

```C++
ddc::BoundCond constexpr SplineVparBoundary;
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/geometryVparMu/geometry/geometry.hpp`

