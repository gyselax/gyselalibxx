

# File geometry.hpp



[**FileList**](files.md) **>** [**geometry**](dir_3d8ac113f1c21fd7bc019cd952574dfc.md) **>** [**geometry.hpp**](geometryXVx_2geometry_2geometry_8hpp.md)

[Go to the source code of this file](geometryXVx_2geometry_2geometry_8hpp_source.md)



* `#include <ddc/ddc.hpp>`
* `#include <ddc/kernels/splines.hpp>`
* `#include "ddc_alias_inline_functions.hpp"`
* `#include "ddc_aliases.hpp"`
* `#include "ddc_helper.hpp"`
* `#include "moments.hpp"`
* `#include "non_uniform_interpolation_points.hpp"`
* `#include "species_info.hpp"`















## Classes

| Type | Name |
| ---: | :--- |
| struct | [**BSplinesVx**](structBSplinesVx.md) <br> |
| struct | [**BSplinesX**](structBSplinesX.md) <br> |
| class | [**GeometryXVx**](classGeometryXVx.md) <br>_A class providing aliases for useful subindex ranges of the geometry. It is used as template parameter for generic dimensionality-agnostic operators such as advections._  |
| struct | [**GridMom**](structGridMom.md) <br> |
| struct | [**GridVx**](structGridVx.md) <br> |
| struct | [**GridX**](structGridX.md) <br> |
| struct | [**T**](structT.md) <br>_A class which describes the real space in the temporal direction._  |
| struct | [**Vx**](structVx.md) <br>_Define non periodic real_ [_**X**_](structX.md) _velocity dimension._ |
| struct | [**X**](structX.md) <br>_Define non periodic real_ [_**X**_](structX.md) _dimension._ |


## Public Types

| Type | Name |
| ---: | :--- |
| typedef ConstField&lt; ElementType, IdxRangeBSX &gt; | [**BSConstFieldX**](#typedef-bsconstfieldx)  <br> |
| typedef FieldMem&lt; ElementType, IdxRangeBSX &gt; | [**BSFieldMemX**](#typedef-bsfieldmemx)  <br> |
| typedef Field&lt; ElementType, IdxRangeBSX &gt; | [**BSFieldX**](#typedef-bsfieldx)  <br> |
| typedef ConstField&lt; ElementType, IdxRangeSpMom &gt; | [**ConstFieldSpMom**](#typedef-constfieldspmom)  <br> |
| typedef ConstField&lt; ElementType, IdxRangeSpMomX &gt; | [**ConstFieldSpMomX**](#typedef-constfieldspmomx)  <br> |
| typedef ConstField&lt; ElementType, IdxRangeSpVx &gt; | [**ConstFieldSpVx**](#typedef-constfieldspvx)  <br> |
| typedef ConstField&lt; ElementType, IdxRangeSpX &gt; | [**ConstFieldSpX**](#typedef-constfieldspx)  <br> |
| typedef ConstField&lt; ElementType, IdxRangeSpXVx &gt; | [**ConstFieldSpXVx**](#typedef-constfieldspxvx)  <br> |
| typedef Field&lt; ElementType const, IdxRangeVx &gt; | [**ConstFieldVx**](#typedef-constfieldvx)  <br> |
| typedef Field&lt; ElementType const, IdxRangeX &gt; | [**ConstFieldX**](#typedef-constfieldx)  <br> |
| typedef Coord&lt; [**T**](structT.md) &gt; | [**CoordT**](#typedef-coordt)  <br> |
| typedef Coord&lt; [**Vx**](structVx.md) &gt; | [**CoordVx**](#typedef-coordvx)  <br> |
| typedef Coord&lt; [**X**](structX.md) &gt; | [**CoordX**](#typedef-coordx)  <br> |
| typedef Coord&lt; [**X**](structX.md), [**Vx**](structVx.md) &gt; | [**CoordXVx**](#typedef-coordxvx)  <br> |
| typedef BSConstFieldX&lt; double &gt; | [**DBSConstFieldX**](#typedef-dbsconstfieldx)  <br> |
| typedef BSFieldMemX&lt; double &gt; | [**DBSFieldMemX**](#typedef-dbsfieldmemx)  <br> |
| typedef BSFieldX&lt; double &gt; | [**DBSFieldX**](#typedef-dbsfieldx)  <br> |
| typedef ConstFieldSpMom&lt; double &gt; | [**DConstFieldSpMom**](#typedef-dconstfieldspmom)  <br> |
| typedef ConstFieldSpMomX&lt; double &gt; | [**DConstFieldSpMomX**](#typedef-dconstfieldspmomx)  <br> |
| typedef ConstFieldSpVx&lt; double &gt; | [**DConstFieldSpVx**](#typedef-dconstfieldspvx)  <br> |
| typedef ConstFieldSpX&lt; double &gt; | [**DConstFieldSpX**](#typedef-dconstfieldspx)  <br> |
| typedef ConstFieldSpXVx&lt; double &gt; | [**DConstFieldSpXVx**](#typedef-dconstfieldspxvx)  <br> |
| typedef ConstFieldVx&lt; double &gt; | [**DConstFieldVx**](#typedef-dconstfieldvx)  <br> |
| typedef ConstFieldX&lt; double &gt; | [**DConstFieldX**](#typedef-dconstfieldx)  <br> |
| typedef FieldMemSpMom&lt; double &gt; | [**DFieldMemSpMom**](#typedef-dfieldmemspmom)  <br> |
| typedef FieldMemSpMomX&lt; double &gt; | [**DFieldMemSpMomX**](#typedef-dfieldmemspmomx)  <br> |
| typedef FieldMemSpVx&lt; double &gt; | [**DFieldMemSpVx**](#typedef-dfieldmemspvx)  <br> |
| typedef FieldMemSpX&lt; double &gt; | [**DFieldMemSpX**](#typedef-dfieldmemspx)  <br> |
| typedef FieldMemSpXVx&lt; double &gt; | [**DFieldMemSpXVx**](#typedef-dfieldmemspxvx)  <br> |
| typedef FieldMemVx&lt; double &gt; | [**DFieldMemVx**](#typedef-dfieldmemvx)  <br> |
| typedef FieldMemX&lt; double &gt; | [**DFieldMemX**](#typedef-dfieldmemx)  <br> |
| typedef FieldSpMom&lt; double &gt; | [**DFieldSpMom**](#typedef-dfieldspmom)  <br> |
| typedef FieldSpMomX&lt; double &gt; | [**DFieldSpMomX**](#typedef-dfieldspmomx)  <br> |
| typedef FieldSpVx&lt; double &gt; | [**DFieldSpVx**](#typedef-dfieldspvx)  <br> |
| typedef FieldSpX&lt; double &gt; | [**DFieldSpX**](#typedef-dfieldspx)  <br> |
| typedef FieldSpXVx&lt; double &gt; | [**DFieldSpXVx**](#typedef-dfieldspxvx)  <br> |
| typedef FieldVx&lt; double &gt; | [**DFieldVx**](#typedef-dfieldvx)  <br> |
| typedef FieldX&lt; double &gt; | [**DFieldX**](#typedef-dfieldx)  <br> |
| typedef FieldMem&lt; ElementType, IdxRangeSpMom &gt; | [**FieldMemSpMom**](#typedef-fieldmemspmom)  <br> |
| typedef FieldMem&lt; ElementType, IdxRangeSpMomX &gt; | [**FieldMemSpMomX**](#typedef-fieldmemspmomx)  <br> |
| typedef FieldMem&lt; ElementType, IdxRangeSpVx &gt; | [**FieldMemSpVx**](#typedef-fieldmemspvx)  <br> |
| typedef FieldMem&lt; ElementType, IdxRangeSpX &gt; | [**FieldMemSpX**](#typedef-fieldmemspx)  <br> |
| typedef FieldMem&lt; ElementType, IdxRangeSpXVx &gt; | [**FieldMemSpXVx**](#typedef-fieldmemspxvx)  <br> |
| typedef FieldMem&lt; ElementType, IdxRangeVx &gt; | [**FieldMemVx**](#typedef-fieldmemvx)  <br> |
| typedef FieldMem&lt; ElementType, IdxRangeX &gt; | [**FieldMemX**](#typedef-fieldmemx)  <br> |
| typedef Field&lt; ElementType, IdxRangeSpMom &gt; | [**FieldSpMom**](#typedef-fieldspmom)  <br> |
| typedef Field&lt; ElementType, IdxRangeSpMomX &gt; | [**FieldSpMomX**](#typedef-fieldspmomx)  <br> |
| typedef Field&lt; ElementType, IdxRangeSpVx &gt; | [**FieldSpVx**](#typedef-fieldspvx)  <br> |
| typedef Field&lt; ElementType, IdxRangeSpX &gt; | [**FieldSpX**](#typedef-fieldspx)  <br> |
| typedef Field&lt; ElementType, IdxRangeSpXVx &gt; | [**FieldSpXVx**](#typedef-fieldspxvx)  <br> |
| typedef Field&lt; ElementType, IdxRangeVx &gt; | [**FieldVx**](#typedef-fieldvx)  <br> |
| typedef Field&lt; ElementType, IdxRangeX &gt; | [**FieldX**](#typedef-fieldx)  <br> |
| typedef Idx&lt; [**GridMom**](structGridMom.md) &gt; | [**IdxMom**](#typedef-idxmom)  <br> |
| typedef IdxRange&lt; [**BSplinesVx**](structBSplinesVx.md) &gt; | [**IdxRangeBSVx**](#typedef-idxrangebsvx)  <br> |
| typedef IdxRange&lt; [**BSplinesX**](structBSplinesX.md) &gt; | [**IdxRangeBSX**](#typedef-idxrangebsx)  <br> |
| typedef IdxRange&lt; [**GridMom**](structGridMom.md) &gt; | [**IdxRangeMom**](#typedef-idxrangemom)  <br> |
| typedef IdxRange&lt; [**Species**](structSpecies.md), [**GridMom**](structGridMom.md) &gt; | [**IdxRangeSpMom**](#typedef-idxrangespmom)  <br> |
| typedef IdxRange&lt; [**Species**](structSpecies.md), [**GridMom**](structGridMom.md), [**GridX**](structGridX.md) &gt; | [**IdxRangeSpMomX**](#typedef-idxrangespmomx)  <br> |
| typedef IdxRange&lt; [**Species**](structSpecies.md), [**GridVx**](structGridVx.md) &gt; | [**IdxRangeSpVx**](#typedef-idxrangespvx)  <br> |
| typedef IdxRange&lt; [**Species**](structSpecies.md), [**GridX**](structGridX.md) &gt; | [**IdxRangeSpX**](#typedef-idxrangespx)  <br> |
| typedef IdxRange&lt; [**Species**](structSpecies.md), [**GridX**](structGridX.md), [**GridVx**](structGridVx.md) &gt; | [**IdxRangeSpXVx**](#typedef-idxrangespxvx)  <br> |
| typedef IdxRange&lt; [**GridVx**](structGridVx.md) &gt; | [**IdxRangeVx**](#typedef-idxrangevx)  <br> |
| typedef IdxRange&lt; [**GridX**](structGridX.md) &gt; | [**IdxRangeX**](#typedef-idxrangex)  <br> |
| typedef IdxRange&lt; [**GridX**](structGridX.md), [**GridVx**](structGridVx.md) &gt; | [**IdxRangeXVx**](#typedef-idxrangexvx)  <br> |
| typedef Idx&lt; [**Species**](structSpecies.md), [**GridMom**](structGridMom.md) &gt; | [**IdxSpMom**](#typedef-idxspmom)  <br> |
| typedef Idx&lt; [**Species**](structSpecies.md), [**GridMom**](structGridMom.md), [**GridX**](structGridX.md) &gt; | [**IdxSpMomX**](#typedef-idxspmomx)  <br> |
| typedef Idx&lt; [**Species**](structSpecies.md), [**GridVx**](structGridVx.md) &gt; | [**IdxSpVx**](#typedef-idxspvx)  <br> |
| typedef Idx&lt; [**Species**](structSpecies.md), [**GridX**](structGridX.md) &gt; | [**IdxSpX**](#typedef-idxspx)  <br> |
| typedef Idx&lt; [**Species**](structSpecies.md), [**GridX**](structGridX.md), [**GridVx**](structGridVx.md) &gt; | [**IdxSpXVx**](#typedef-idxspxvx)  <br> |
| typedef IdxStep&lt; [**GridMom**](structGridMom.md) &gt; | [**IdxStepMom**](#typedef-idxstepmom)  <br> |
| typedef IdxStep&lt; [**Species**](structSpecies.md), [**GridMom**](structGridMom.md) &gt; | [**IdxStepSpMom**](#typedef-idxstepspmom)  <br> |
| typedef IdxStep&lt; [**Species**](structSpecies.md), [**GridMom**](structGridMom.md), [**GridX**](structGridX.md) &gt; | [**IdxStepSpMomX**](#typedef-idxstepspmomx)  <br> |
| typedef IdxStep&lt; [**Species**](structSpecies.md), [**GridVx**](structGridVx.md) &gt; | [**IdxStepSpVx**](#typedef-idxstepspvx)  <br> |
| typedef IdxStep&lt; [**Species**](structSpecies.md), [**GridX**](structGridX.md) &gt; | [**IdxStepSpX**](#typedef-idxstepspx)  <br> |
| typedef IdxStep&lt; [**Species**](structSpecies.md), [**GridX**](structGridX.md), [**GridVx**](structGridVx.md) &gt; | [**IdxStepSpXVx**](#typedef-idxstepspxvx)  <br> |
| typedef IdxStep&lt; [**GridVx**](structGridVx.md) &gt; | [**IdxStepVx**](#typedef-idxstepvx)  <br> |
| typedef IdxStep&lt; [**GridX**](structGridX.md) &gt; | [**IdxStepX**](#typedef-idxstepx)  <br> |
| typedef IdxStep&lt; [**GridX**](structGridX.md), [**GridVx**](structGridVx.md) &gt; | [**IdxStepXVx**](#typedef-idxstepxvx)  <br> |
| typedef Idx&lt; [**GridVx**](structGridVx.md) &gt; | [**IdxVx**](#typedef-idxvx)  <br> |
| typedef Idx&lt; [**GridX**](structGridX.md) &gt; | [**IdxX**](#typedef-idxx)  <br> |
| typedef Idx&lt; [**GridX**](structGridX.md), [**GridVx**](structGridVx.md) &gt; | [**IdxXVx**](#typedef-idxxvx)  <br> |
| typedef ddc::GrevilleInterpolationPoints&lt; [**BSplinesVx**](structBSplinesVx.md), SplineVxBoundary, SplineVxBoundary &gt; | [**SplineInterpPointsVx**](#typedef-splineinterppointsvx)  <br> |
| typedef ddc::GrevilleInterpolationPoints&lt; [**BSplinesX**](structBSplinesX.md), SplineXBoundary, SplineXBoundary &gt; | [**SplineInterpPointsX**](#typedef-splineinterppointsx)  <br> |
| typedef ddc::SplineBuilder&lt; Kokkos::DefaultExecutionSpace, Kokkos::DefaultExecutionSpace::memory\_space, [**BSplinesVx**](structBSplinesVx.md), [**GridVx**](structGridVx.md), SplineVxBoundary, SplineVxBoundary, ddc::SplineSolver::LAPACK, [**GridX**](structGridX.md), [**GridVx**](structGridVx.md) &gt; | [**SplineVxBuilder**](#typedef-splinevxbuilder)  <br> |
| typedef ddc::SplineBuilder&lt; Kokkos::DefaultHostExecutionSpace, Kokkos::HostSpace, [**BSplinesVx**](structBSplinesVx.md), [**GridVx**](structGridVx.md), SplineVxBoundary, SplineVxBoundary, ddc::SplineSolver::LAPACK, [**GridVx**](structGridVx.md) &gt; | [**SplineVxBuilder\_1d**](#typedef-splinevxbuilder_1d)  <br> |
| typedef ddc::SplineEvaluator&lt; Kokkos::DefaultExecutionSpace, Kokkos::DefaultExecutionSpace::memory\_space, [**BSplinesVx**](structBSplinesVx.md), [**GridVx**](structGridVx.md), ddc::ConstantExtrapolationRule&lt; [**Vx**](structVx.md) &gt;, ddc::ConstantExtrapolationRule&lt; [**Vx**](structVx.md) &gt;, [**GridX**](structGridX.md), [**GridVx**](structGridVx.md) &gt; | [**SplineVxEvaluator**](#typedef-splinevxevaluator)  <br> |
| typedef ddc::SplineEvaluator&lt; Kokkos::DefaultHostExecutionSpace, Kokkos::HostSpace, [**BSplinesVx**](structBSplinesVx.md), [**GridVx**](structGridVx.md), ddc::ConstantExtrapolationRule&lt; [**Vx**](structVx.md) &gt;, ddc::ConstantExtrapolationRule&lt; [**Vx**](structVx.md) &gt;, [**GridVx**](structGridVx.md) &gt; | [**SplineVxEvaluator\_1d**](#typedef-splinevxevaluator_1d)  <br> |
| typedef ddc::SplineBuilder&lt; Kokkos::DefaultExecutionSpace, Kokkos::DefaultExecutionSpace::memory\_space, [**BSplinesX**](structBSplinesX.md), [**GridX**](structGridX.md), SplineXBoundary, SplineXBoundary, ddc::SplineSolver::LAPACK, [**GridX**](structGridX.md), [**GridVx**](structGridVx.md) &gt; | [**SplineXBuilder**](#typedef-splinexbuilder)  <br> |
| typedef ddc::SplineBuilder&lt; Kokkos::DefaultExecutionSpace, Kokkos::DefaultExecutionSpace::memory\_space, [**BSplinesX**](structBSplinesX.md), [**GridX**](structGridX.md), SplineXBoundary, SplineXBoundary, ddc::SplineSolver::LAPACK, [**GridX**](structGridX.md) &gt; | [**SplineXBuilder\_1d**](#typedef-splinexbuilder_1d)  <br> |
| typedef ddc::SplineEvaluator&lt; Kokkos::DefaultExecutionSpace, Kokkos::DefaultExecutionSpace::memory\_space, [**BSplinesX**](structBSplinesX.md), [**GridX**](structGridX.md), ddc::ConstantExtrapolationRule&lt; [**X**](structX.md) &gt;, ddc::ConstantExtrapolationRule&lt; [**X**](structX.md) &gt;, [**GridX**](structGridX.md), [**GridVx**](structGridVx.md) &gt; | [**SplineXEvaluator**](#typedef-splinexevaluator)  <br> |
| typedef ddc::SplineEvaluator&lt; Kokkos::DefaultExecutionSpace, Kokkos::DefaultExecutionSpace::memory\_space, [**BSplinesX**](structBSplinesX.md), [**GridX**](structGridX.md), ddc::ConstantExtrapolationRule&lt; [**X**](structX.md) &gt;, ddc::ConstantExtrapolationRule&lt; [**X**](structX.md) &gt;, [**GridX**](structGridX.md) &gt; | [**SplineXEvaluator\_1d**](#typedef-splinexevaluator_1d)  <br> |




## Public Attributes

| Type | Name |
| ---: | :--- |
|  int constexpr | [**BSDegreeVx**](#variable-bsdegreevx)   = `3`<br> |
|  int constexpr | [**BSDegreeX**](#variable-bsdegreex)   = `3`<br> |
|  bool constexpr | [**BsplineOnUniformCellsVx**](#variable-bsplineonuniformcellsvx)   = `true`<br> |
|  bool constexpr | [**BsplineOnUniformCellsX**](#variable-bsplineonuniformcellsx)   = `true`<br> |
|  auto constexpr | [**SplineVxBoundary**](#variable-splinevxboundary)   = `ddc::BoundCond::HERMITE`<br> |
|  auto constexpr | [**SplineXBoundary**](#variable-splinexboundary)   = `[**X::PERIODIC**](structX.md#variable-periodic) ? ddc::BoundCond::PERIODIC : ddc::BoundCond::GREVILLE`<br> |












































## Public Types Documentation




### typedef BSConstFieldX 

```C++
using BSConstFieldX =  ConstField<ElementType, IdxRangeBSX>;
```




<hr>



### typedef BSFieldMemX 

```C++
using BSFieldMemX =  FieldMem<ElementType, IdxRangeBSX>;
```




<hr>



### typedef BSFieldX 

```C++
using BSFieldX =  Field<ElementType, IdxRangeBSX>;
```




<hr>



### typedef ConstFieldSpMom 

```C++
using ConstFieldSpMom =  ConstField<ElementType, IdxRangeSpMom>;
```




<hr>



### typedef ConstFieldSpMomX 

```C++
using ConstFieldSpMomX =  ConstField<ElementType, IdxRangeSpMomX>;
```




<hr>



### typedef ConstFieldSpVx 

```C++
using ConstFieldSpVx =  ConstField<ElementType, IdxRangeSpVx>;
```




<hr>



### typedef ConstFieldSpX 

```C++
using ConstFieldSpX =  ConstField<ElementType, IdxRangeSpX>;
```




<hr>



### typedef ConstFieldSpXVx 

```C++
using ConstFieldSpXVx =  ConstField<ElementType, IdxRangeSpXVx>;
```




<hr>



### typedef ConstFieldVx 

```C++
using ConstFieldVx =  Field<ElementType const, IdxRangeVx>;
```




<hr>



### typedef ConstFieldX 

```C++
using ConstFieldX =  Field<ElementType const, IdxRangeX>;
```




<hr>



### typedef CoordT 

```C++
using CoordT =  Coord<T>;
```




<hr>



### typedef CoordVx 

```C++
using CoordVx =  Coord<Vx>;
```




<hr>



### typedef CoordX 

```C++
using CoordX =  Coord<X>;
```




<hr>



### typedef CoordXVx 

```C++
using CoordXVx =  Coord<X, Vx>;
```




<hr>



### typedef DBSConstFieldX 

```C++
using DBSConstFieldX =  BSConstFieldX<double>;
```




<hr>



### typedef DBSFieldMemX 

```C++
using DBSFieldMemX =  BSFieldMemX<double>;
```




<hr>



### typedef DBSFieldX 

```C++
using DBSFieldX =  BSFieldX<double>;
```




<hr>



### typedef DConstFieldSpMom 

```C++
typedef ConstFieldSpMom< double > DConstFieldSpMom;
```




<hr>



### typedef DConstFieldSpMomX 

```C++
using DConstFieldSpMomX =  ConstFieldSpMomX<double>;
```




<hr>



### typedef DConstFieldSpVx 

```C++
using DConstFieldSpVx =  ConstFieldSpVx<double>;
```




<hr>



### typedef DConstFieldSpX 

```C++
using DConstFieldSpX =  ConstFieldSpX<double>;
```




<hr>



### typedef DConstFieldSpXVx 

```C++
using DConstFieldSpXVx =  ConstFieldSpXVx<double>;
```




<hr>



### typedef DConstFieldVx 

```C++
using DConstFieldVx =  ConstFieldVx<double>;
```




<hr>



### typedef DConstFieldX 

```C++
using DConstFieldX =  ConstFieldX<double>;
```




<hr>



### typedef DFieldMemSpMom 

```C++
using DFieldMemSpMom =  FieldMemSpMom<double>;
```




<hr>



### typedef DFieldMemSpMomX 

```C++
using DFieldMemSpMomX =  FieldMemSpMomX<double>;
```




<hr>



### typedef DFieldMemSpVx 

```C++
using DFieldMemSpVx =  FieldMemSpVx<double>;
```




<hr>



### typedef DFieldMemSpX 

```C++
using DFieldMemSpX =  FieldMemSpX<double>;
```




<hr>



### typedef DFieldMemSpXVx 

```C++
using DFieldMemSpXVx =  FieldMemSpXVx<double>;
```




<hr>



### typedef DFieldMemVx 

```C++
using DFieldMemVx =  FieldMemVx<double>;
```




<hr>



### typedef DFieldMemX 

```C++
using DFieldMemX =  FieldMemX<double>;
```




<hr>



### typedef DFieldSpMom 

```C++
using DFieldSpMom =  FieldSpMom<double>;
```




<hr>



### typedef DFieldSpMomX 

```C++
using DFieldSpMomX =  FieldSpMomX<double>;
```




<hr>



### typedef DFieldSpVx 

```C++
using DFieldSpVx =  FieldSpVx<double>;
```




<hr>



### typedef DFieldSpX 

```C++
using DFieldSpX =  FieldSpX<double>;
```




<hr>



### typedef DFieldSpXVx 

```C++
using DFieldSpXVx =  FieldSpXVx<double>;
```




<hr>



### typedef DFieldVx 

```C++
using DFieldVx =  FieldVx<double>;
```




<hr>



### typedef DFieldX 

```C++
using DFieldX =  FieldX<double>;
```




<hr>



### typedef FieldMemSpMom 

```C++
using FieldMemSpMom =  FieldMem<ElementType, IdxRangeSpMom>;
```




<hr>



### typedef FieldMemSpMomX 

```C++
using FieldMemSpMomX =  FieldMem<ElementType, IdxRangeSpMomX>;
```




<hr>



### typedef FieldMemSpVx 

```C++
using FieldMemSpVx =  FieldMem<ElementType, IdxRangeSpVx>;
```




<hr>



### typedef FieldMemSpX 

```C++
using FieldMemSpX =  FieldMem<ElementType, IdxRangeSpX>;
```




<hr>



### typedef FieldMemSpXVx 

```C++
using FieldMemSpXVx =  FieldMem<ElementType, IdxRangeSpXVx>;
```




<hr>



### typedef FieldMemVx 

```C++
using FieldMemVx =  FieldMem<ElementType, IdxRangeVx>;
```




<hr>



### typedef FieldMemX 

```C++
using FieldMemX =  FieldMem<ElementType, IdxRangeX>;
```




<hr>



### typedef FieldSpMom 

```C++
using FieldSpMom =  Field<ElementType, IdxRangeSpMom>;
```




<hr>



### typedef FieldSpMomX 

```C++
using FieldSpMomX =  Field<ElementType, IdxRangeSpMomX>;
```




<hr>



### typedef FieldSpVx 

```C++
using FieldSpVx =  Field<ElementType, IdxRangeSpVx>;
```




<hr>



### typedef FieldSpX 

```C++
using FieldSpX =  Field<ElementType, IdxRangeSpX>;
```




<hr>



### typedef FieldSpXVx 

```C++
using FieldSpXVx =  Field<ElementType, IdxRangeSpXVx>;
```




<hr>



### typedef FieldVx 

```C++
using FieldVx =  Field<ElementType, IdxRangeVx>;
```




<hr>



### typedef FieldX 

```C++
using FieldX =  Field<ElementType, IdxRangeX>;
```




<hr>



### typedef IdxMom 

```C++
using IdxMom =  Idx<GridMom>;
```




<hr>



### typedef IdxRangeBSVx 

```C++
using IdxRangeBSVx =  IdxRange<BSplinesVx>;
```




<hr>



### typedef IdxRangeBSX 

```C++
using IdxRangeBSX =  IdxRange<BSplinesX>;
```




<hr>



### typedef IdxRangeMom 

```C++
using IdxRangeMom =  IdxRange<GridMom>;
```




<hr>



### typedef IdxRangeSpMom 

```C++
using IdxRangeSpMom =  IdxRange<Species, GridMom>;
```




<hr>



### typedef IdxRangeSpMomX 

```C++
using IdxRangeSpMomX =  IdxRange<Species, GridMom, GridX>;
```




<hr>



### typedef IdxRangeSpVx 

```C++
using IdxRangeSpVx =  IdxRange<Species, GridVx>;
```




<hr>



### typedef IdxRangeSpX 

```C++
using IdxRangeSpX =  IdxRange<Species, GridX>;
```




<hr>



### typedef IdxRangeSpXVx 

```C++
using IdxRangeSpXVx =  IdxRange<Species, GridX, GridVx>;
```




<hr>



### typedef IdxRangeVx 

```C++
using IdxRangeVx =  IdxRange<GridVx>;
```




<hr>



### typedef IdxRangeX 

```C++
using IdxRangeX =  IdxRange<GridX>;
```




<hr>



### typedef IdxRangeXVx 

```C++
using IdxRangeXVx =  IdxRange<GridX, GridVx>;
```




<hr>



### typedef IdxSpMom 

```C++
using IdxSpMom =  Idx<Species, GridMom>;
```




<hr>



### typedef IdxSpMomX 

```C++
using IdxSpMomX =  Idx<Species, GridMom, GridX>;
```




<hr>



### typedef IdxSpVx 

```C++
using IdxSpVx =  Idx<Species, GridVx>;
```




<hr>



### typedef IdxSpX 

```C++
using IdxSpX =  Idx<Species, GridX>;
```




<hr>



### typedef IdxSpXVx 

```C++
using IdxSpXVx =  Idx<Species, GridX, GridVx>;
```




<hr>



### typedef IdxStepMom 

```C++
using IdxStepMom =  IdxStep<GridMom>;
```




<hr>



### typedef IdxStepSpMom 

```C++
using IdxStepSpMom =  IdxStep<Species, GridMom>;
```




<hr>



### typedef IdxStepSpMomX 

```C++
using IdxStepSpMomX =  IdxStep<Species, GridMom, GridX>;
```




<hr>



### typedef IdxStepSpVx 

```C++
using IdxStepSpVx =  IdxStep<Species, GridVx>;
```




<hr>



### typedef IdxStepSpX 

```C++
using IdxStepSpX =  IdxStep<Species, GridX>;
```




<hr>



### typedef IdxStepSpXVx 

```C++
using IdxStepSpXVx =  IdxStep<Species, GridX, GridVx>;
```




<hr>



### typedef IdxStepVx 

```C++
using IdxStepVx =  IdxStep<GridVx>;
```




<hr>



### typedef IdxStepX 

```C++
using IdxStepX =  IdxStep<GridX>;
```




<hr>



### typedef IdxStepXVx 

```C++
using IdxStepXVx =  IdxStep<GridX, GridVx>;
```




<hr>



### typedef IdxVx 

```C++
using IdxVx =  Idx<GridVx>;
```




<hr>



### typedef IdxX 

```C++
using IdxX =  Idx<GridX>;
```




<hr>



### typedef IdxXVx 

```C++
using IdxXVx =  Idx<GridX, GridVx>;
```




<hr>



### typedef SplineInterpPointsVx 

```C++
using SplineInterpPointsVx =  ddc::GrevilleInterpolationPoints<BSplinesVx, SplineVxBoundary, SplineVxBoundary>;
```




<hr>



### typedef SplineInterpPointsX 

```C++
using SplineInterpPointsX =  ddc::GrevilleInterpolationPoints<BSplinesX, SplineXBoundary, SplineXBoundary>;
```




<hr>



### typedef SplineVxBuilder 

```C++
using SplineVxBuilder =  ddc::SplineBuilder< Kokkos::DefaultExecutionSpace, Kokkos::DefaultExecutionSpace::memory_space, BSplinesVx, GridVx, SplineVxBoundary, SplineVxBoundary, ddc::SplineSolver::LAPACK, GridX, GridVx>;
```




<hr>



### typedef SplineVxBuilder\_1d 

```C++
using SplineVxBuilder_1d =  ddc::SplineBuilder< Kokkos::DefaultHostExecutionSpace, Kokkos::HostSpace, BSplinesVx, GridVx, SplineVxBoundary, SplineVxBoundary, ddc::SplineSolver::LAPACK, GridVx>;
```




<hr>



### typedef SplineVxEvaluator 

```C++
using SplineVxEvaluator =  ddc::SplineEvaluator< Kokkos::DefaultExecutionSpace, Kokkos::DefaultExecutionSpace::memory_space, BSplinesVx, GridVx, ddc::ConstantExtrapolationRule<Vx>, ddc::ConstantExtrapolationRule<Vx>, GridX, GridVx>;
```




<hr>



### typedef SplineVxEvaluator\_1d 

```C++
using SplineVxEvaluator_1d =  ddc::SplineEvaluator< Kokkos::DefaultHostExecutionSpace, Kokkos::HostSpace, BSplinesVx, GridVx, ddc::ConstantExtrapolationRule<Vx>, ddc::ConstantExtrapolationRule<Vx>, GridVx>;
```




<hr>



### typedef SplineXBuilder 

```C++
using SplineXBuilder =  ddc::SplineBuilder< Kokkos::DefaultExecutionSpace, Kokkos::DefaultExecutionSpace::memory_space, BSplinesX, GridX, SplineXBoundary, SplineXBoundary, ddc::SplineSolver::LAPACK, GridX, GridVx>;
```




<hr>



### typedef SplineXBuilder\_1d 

```C++
using SplineXBuilder_1d =  ddc::SplineBuilder< Kokkos::DefaultExecutionSpace, Kokkos::DefaultExecutionSpace::memory_space, BSplinesX, GridX, SplineXBoundary, SplineXBoundary, ddc::SplineSolver::LAPACK, GridX>;
```




<hr>



### typedef SplineXEvaluator 

```C++
using SplineXEvaluator =  ddc::SplineEvaluator< Kokkos::DefaultExecutionSpace, Kokkos::DefaultExecutionSpace::memory_space, BSplinesX, GridX, ddc::ConstantExtrapolationRule<X>, ddc::ConstantExtrapolationRule<X>, GridX, GridVx>;
```




<hr>



### typedef SplineXEvaluator\_1d 

```C++
using SplineXEvaluator_1d =  ddc::SplineEvaluator< Kokkos::DefaultExecutionSpace, Kokkos::DefaultExecutionSpace::memory_space, BSplinesX, GridX, ddc::ConstantExtrapolationRule<X>, ddc::ConstantExtrapolationRule<X>, GridX>;
```




<hr>
## Public Attributes Documentation




### variable BSDegreeVx 

```C++
int constexpr BSDegreeVx;
```




<hr>



### variable BSDegreeX 

```C++
int constexpr BSDegreeX;
```




<hr>



### variable BsplineOnUniformCellsVx 

```C++
bool constexpr BsplineOnUniformCellsVx;
```




<hr>



### variable BsplineOnUniformCellsX 

```C++
bool constexpr BsplineOnUniformCellsX;
```




<hr>



### variable SplineVxBoundary 

```C++
auto constexpr SplineVxBoundary;
```




<hr>



### variable SplineXBoundary 

```C++
auto constexpr SplineXBoundary;
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/geometryXVx/geometry/geometry.hpp`

