

# File geometry.hpp



[**FileList**](files.md) **>** [**geometry**](dir_7ddd2963f3e4609fce61e92aa9c5ff14.md) **>** [**geometry.hpp**](geometryXYVxVy_2geometry_2geometry_8hpp.md)

[Go to the source code of this file](geometryXYVxVy_2geometry_2geometry_8hpp_source.md)



* `#include <ddc/ddc.hpp>`
* `#include <ddc/kernels/splines.hpp>`
* `#include "ddc_alias_inline_functions.hpp"`
* `#include "ddc_aliases.hpp"`
* `#include "ddc_helper.hpp"`
* `#include "mpilayout.hpp"`
* `#include "species_info.hpp"`















## Classes

| Type | Name |
| ---: | :--- |
| struct | [**BSplinesVx**](structBSplinesVx.md) <br> |
| struct | [**BSplinesVy**](structBSplinesVy.md) <br> |
| struct | [**BSplinesX**](structBSplinesX.md) <br> |
| struct | [**BSplinesY**](structBSplinesY.md) <br> |
| class | [**GeometryVxVyXY**](classGeometryVxVyXY.md) <br>_A class providing aliases for useful subindex ranges of the geometry when the data is saved with the velocity dimensions distributed across MPI ranks. It is used as template parameter for generic dimensionality-agnostic operators such as advections._  |
| class | [**GeometryXYVxVy**](classGeometryXYVxVy.md) <br>_A class providing aliases for useful subindex ranges of the geometry when the data is saved with the spatial dimensions distributed across MPI ranks. It is used as template parameter for generic dimensionality-agnostic operators such as advections._  |
| struct | [**GridVx**](structGridVx.md) <br> |
| struct | [**GridVy**](structGridVy.md) <br> |
| struct | [**GridX**](structGridX.md) <br> |
| struct | [**GridY**](structGridY.md) <br> |
| struct | [**Vx**](structVx.md) <br>_Define non periodic real_ [_**X**_](structX.md) _velocity dimension._ |
| struct | [**Vy**](structVy.md) <br>_Define non periodic real_ [_**Y**_](structY.md) _velocity dimension._ |
| struct | [**X**](structX.md) <br>_Define non periodic real_ [_**X**_](structX.md) _dimension._ |
| struct | [**Y**](structY.md) <br>_Define non periodic real_ [_**Y**_](structY.md) _dimension._ |


## Public Types

| Type | Name |
| ---: | :--- |
| typedef Field&lt; ElementType const, IdxRangeBSXY &gt; | [**BSConstFieldXY**](#typedef-bsconstfieldxy)  <br> |
| typedef Field&lt; ElementType const, IdxRangeSpVxVy &gt; | [**ConstFieldSpVxVy**](#typedef-constfieldspvxvy)  <br> |
| typedef Field&lt; ElementType const, IdxRangeSpVxVyXY &gt; | [**ConstFieldSpVxVyXY**](#typedef-constfieldspvxvyxy)  <br> |
| typedef Field&lt; ElementType const, IdxRangeSpXYVxVy &gt; | [**ConstFieldSpXYVxVy**](#typedef-constfieldspxyvxvy)  <br> |
| typedef Field&lt; ElementType const, IdxRangeVx &gt; | [**ConstFieldVx**](#typedef-constfieldvx)  <br> |
| typedef Field&lt; ElementType const, IdxRangeVxVy &gt; | [**ConstFieldVxVy**](#typedef-constfieldvxvy)  <br> |
| typedef Field&lt; ElementType const, IdxRangeVy &gt; | [**ConstFieldVy**](#typedef-constfieldvy)  <br> |
| typedef Field&lt; ElementType const, IdxRangeX &gt; | [**ConstFieldX**](#typedef-constfieldx)  <br> |
| typedef Field&lt; ElementType const, IdxRangeXY &gt; | [**ConstFieldXY**](#typedef-constfieldxy)  <br> |
| typedef Field&lt; ElementType const, IdxRangeY &gt; | [**ConstFieldY**](#typedef-constfieldy)  <br> |
| typedef Coord&lt; [**Vx**](structVx.md) &gt; | [**CoordVx**](#typedef-coordvx)  <br> |
| typedef Coord&lt; [**Vy**](structVy.md) &gt; | [**CoordVy**](#typedef-coordvy)  <br> |
| typedef Coord&lt; [**X**](structX.md) &gt; | [**CoordX**](#typedef-coordx)  <br> |
| typedef Coord&lt; [**X**](structX.md), [**Y**](structY.md) &gt; | [**CoordXY**](#typedef-coordxy)  <br> |
| typedef Coord&lt; [**Y**](structY.md) &gt; | [**CoordY**](#typedef-coordy)  <br> |
| typedef BSConstFieldXY&lt; double &gt; | [**DBSConstFieldXY**](#typedef-dbsconstfieldxy)  <br> |
| typedef ConstFieldSpVxVy&lt; double &gt; | [**DConstFieldSpVxVy**](#typedef-dconstfieldspvxvy)  <br> |
| typedef ConstFieldSpVxVyXY&lt; double &gt; | [**DConstFieldSpVxVyXY**](#typedef-dconstfieldspvxvyxy)  <br> |
| typedef ConstFieldSpXYVxVy&lt; double &gt; | [**DConstFieldSpXYVxVy**](#typedef-dconstfieldspxyvxvy)  <br> |
| typedef ConstFieldVxVy&lt; double &gt; | [**DConstFieldVxVy**](#typedef-dconstfieldvxvy)  <br> |
| typedef ConstFieldXY&lt; double &gt; | [**DConstFieldXY**](#typedef-dconstfieldxy)  <br> |
| typedef FieldMemSpVxVy&lt; double &gt; | [**DFieldMemSpVxVy**](#typedef-dfieldmemspvxvy)  <br> |
| typedef FieldMemSpVxVyXY&lt; double &gt; | [**DFieldMemSpVxVyXY**](#typedef-dfieldmemspvxvyxy)  <br> |
| typedef FieldMemSpXYVxVy&lt; double &gt; | [**DFieldMemSpXYVxVy**](#typedef-dfieldmemspxyvxvy)  <br> |
| typedef FieldMemVxVy&lt; double &gt; | [**DFieldMemVxVy**](#typedef-dfieldmemvxvy)  <br> |
| typedef FieldMemX&lt; double &gt; | [**DFieldMemX**](#typedef-dfieldmemx)  <br> |
| typedef FieldMemXY&lt; double &gt; | [**DFieldMemXY**](#typedef-dfieldmemxy)  <br> |
| typedef FieldMemXYVxVy&lt; double &gt; | [**DFieldMemXYVxVy**](#typedef-dfieldmemxyvxvy)  <br> |
| typedef FieldMemY&lt; double &gt; | [**DFieldMemY**](#typedef-dfieldmemy)  <br> |
| typedef FieldSpVxVy&lt; double &gt; | [**DFieldSpVxVy**](#typedef-dfieldspvxvy)  <br> |
| typedef FieldSpVxVyXY&lt; double &gt; | [**DFieldSpVxVyXY**](#typedef-dfieldspvxvyxy)  <br> |
| typedef FieldSpXYVxVy&lt; double &gt; | [**DFieldSpXYVxVy**](#typedef-dfieldspxyvxvy)  <br> |
| typedef FieldVx&lt; double &gt; | [**DFieldVx**](#typedef-dfieldvx)  <br> |
| typedef FieldVxVy&lt; double &gt; | [**DFieldVxVy**](#typedef-dfieldvxvy)  <br> |
| typedef FieldVy&lt; double &gt; | [**DFieldVy**](#typedef-dfieldvy)  <br> |
| typedef FieldX&lt; double &gt; | [**DFieldX**](#typedef-dfieldx)  <br> |
| typedef FieldXY&lt; double &gt; | [**DFieldXY**](#typedef-dfieldxy)  <br> |
| typedef FieldY&lt; double &gt; | [**DFieldY**](#typedef-dfieldy)  <br> |
| typedef FieldMem&lt; ElementType, IdxRangeSpVxVy &gt; | [**FieldMemSpVxVy**](#typedef-fieldmemspvxvy)  <br> |
| typedef FieldMem&lt; ElementType, IdxRangeSpVxVyXY &gt; | [**FieldMemSpVxVyXY**](#typedef-fieldmemspvxvyxy)  <br> |
| typedef FieldMem&lt; ElementType, IdxRangeSpXYVxVy &gt; | [**FieldMemSpXYVxVy**](#typedef-fieldmemspxyvxvy)  <br> |
| typedef FieldMem&lt; ElementType, IdxRangeVx &gt; | [**FieldMemVx**](#typedef-fieldmemvx)  <br> |
| typedef FieldMem&lt; ElementType, IdxRangeVxVy &gt; | [**FieldMemVxVy**](#typedef-fieldmemvxvy)  <br> |
| typedef FieldMem&lt; ElementType, IdxRangeVy &gt; | [**FieldMemVy**](#typedef-fieldmemvy)  <br> |
| typedef FieldMem&lt; ElementType, IdxRangeX &gt; | [**FieldMemX**](#typedef-fieldmemx)  <br> |
| typedef FieldMem&lt; ElementType, IdxRangeXY &gt; | [**FieldMemXY**](#typedef-fieldmemxy)  <br> |
| typedef FieldMem&lt; ElementType, IdxRangeXYVxVy &gt; | [**FieldMemXYVxVy**](#typedef-fieldmemxyvxvy)  <br> |
| typedef FieldMem&lt; ElementType, IdxRangeY &gt; | [**FieldMemY**](#typedef-fieldmemy)  <br> |
| typedef Field&lt; ElementType, IdxRangeSpVxVy &gt; | [**FieldSpVxVy**](#typedef-fieldspvxvy)  <br> |
| typedef Field&lt; ElementType, IdxRangeSpVxVyXY &gt; | [**FieldSpVxVyXY**](#typedef-fieldspvxvyxy)  <br> |
| typedef Field&lt; ElementType, IdxRangeSpXYVxVy &gt; | [**FieldSpXYVxVy**](#typedef-fieldspxyvxvy)  <br> |
| typedef Field&lt; ElementType, IdxRangeVx &gt; | [**FieldVx**](#typedef-fieldvx)  <br> |
| typedef Field&lt; ElementType, IdxRangeVxVy &gt; | [**FieldVxVy**](#typedef-fieldvxvy)  <br> |
| typedef Field&lt; ElementType, IdxRangeVy &gt; | [**FieldVy**](#typedef-fieldvy)  <br> |
| typedef Field&lt; ElementType, IdxRangeX &gt; | [**FieldX**](#typedef-fieldx)  <br> |
| typedef Field&lt; ElementType, IdxRangeXY &gt; | [**FieldXY**](#typedef-fieldxy)  <br> |
| typedef Field&lt; ElementType, IdxRangeY &gt; | [**FieldY**](#typedef-fieldy)  <br> |
| typedef IdxRange&lt; [**BSplinesVx**](structBSplinesVx.md) &gt; | [**IdxRangeBSVx**](#typedef-idxrangebsvx)  <br> |
| typedef IdxRange&lt; [**BSplinesVx**](structBSplinesVx.md), [**BSplinesVy**](structBSplinesVy.md) &gt; | [**IdxRangeBSVxVy**](#typedef-idxrangebsvxvy)  <br> |
| typedef IdxRange&lt; [**BSplinesVy**](structBSplinesVy.md) &gt; | [**IdxRangeBSVy**](#typedef-idxrangebsvy)  <br> |
| typedef IdxRange&lt; [**BSplinesX**](structBSplinesX.md) &gt; | [**IdxRangeBSX**](#typedef-idxrangebsx)  <br> |
| typedef IdxRange&lt; [**BSplinesX**](structBSplinesX.md), [**BSplinesY**](structBSplinesY.md) &gt; | [**IdxRangeBSXY**](#typedef-idxrangebsxy)  <br> |
| typedef IdxRange&lt; [**BSplinesY**](structBSplinesY.md) &gt; | [**IdxRangeBSY**](#typedef-idxrangebsy)  <br> |
| typedef IdxRange&lt; [**Species**](structSpecies.md), [**GridVx**](structGridVx.md), [**GridVy**](structGridVy.md) &gt; | [**IdxRangeSpVxVy**](#typedef-idxrangespvxvy)  <br> |
| typedef IdxRange&lt; [**Species**](structSpecies.md), [**GridVx**](structGridVx.md), [**GridVy**](structGridVy.md), [**GridX**](structGridX.md), [**GridY**](structGridY.md) &gt; | [**IdxRangeSpVxVyXY**](#typedef-idxrangespvxvyxy)  <br> |
| typedef IdxRange&lt; [**Species**](structSpecies.md), [**GridX**](structGridX.md), [**GridY**](structGridY.md), [**GridVx**](structGridVx.md), [**GridVy**](structGridVy.md) &gt; | [**IdxRangeSpXYVxVy**](#typedef-idxrangespxyvxvy)  <br> |
| typedef IdxRange&lt; [**GridVx**](structGridVx.md) &gt; | [**IdxRangeVx**](#typedef-idxrangevx)  <br> |
| typedef IdxRange&lt; [**GridVx**](structGridVx.md), [**GridVy**](structGridVy.md) &gt; | [**IdxRangeVxVy**](#typedef-idxrangevxvy)  <br> |
| typedef IdxRange&lt; [**GridVx**](structGridVx.md), [**GridVy**](structGridVy.md), [**GridX**](structGridX.md), [**GridY**](structGridY.md) &gt; | [**IdxRangeVxVyXY**](#typedef-idxrangevxvyxy)  <br> |
| typedef IdxRange&lt; [**GridVy**](structGridVy.md) &gt; | [**IdxRangeVy**](#typedef-idxrangevy)  <br> |
| typedef IdxRange&lt; [**GridX**](structGridX.md) &gt; | [**IdxRangeX**](#typedef-idxrangex)  <br> |
| typedef IdxRange&lt; [**GridX**](structGridX.md), [**GridY**](structGridY.md) &gt; | [**IdxRangeXY**](#typedef-idxrangexy)  <br> |
| typedef IdxRange&lt; [**GridX**](structGridX.md), [**GridY**](structGridY.md), [**GridVx**](structGridVx.md), [**GridVy**](structGridVy.md) &gt; | [**IdxRangeXYVxVy**](#typedef-idxrangexyvxvy)  <br> |
| typedef IdxRange&lt; [**GridY**](structGridY.md) &gt; | [**IdxRangeY**](#typedef-idxrangey)  <br> |
| typedef Idx&lt; [**Species**](structSpecies.md), [**GridX**](structGridX.md), [**GridY**](structGridY.md), [**GridVx**](structGridVx.md), [**GridVy**](structGridVy.md) &gt; | [**IdxSpXYVxVy**](#typedef-idxspxyvxvy)  <br> |
| typedef IdxStep&lt; [**GridVx**](structGridVx.md) &gt; | [**IdxStepVx**](#typedef-idxstepvx)  <br> |
| typedef IdxStep&lt; [**GridVy**](structGridVy.md) &gt; | [**IdxStepVy**](#typedef-idxstepvy)  <br> |
| typedef IdxStep&lt; [**GridX**](structGridX.md) &gt; | [**IdxStepX**](#typedef-idxstepx)  <br> |
| typedef IdxStep&lt; [**GridY**](structGridY.md) &gt; | [**IdxStepY**](#typedef-idxstepy)  <br> |
| typedef Idx&lt; [**GridVx**](structGridVx.md) &gt; | [**IdxVx**](#typedef-idxvx)  <br> |
| typedef Idx&lt; [**GridVx**](structGridVx.md), [**GridVy**](structGridVy.md) &gt; | [**IdxVxVy**](#typedef-idxvxvy)  <br> |
| typedef Idx&lt; [**GridVy**](structGridVy.md) &gt; | [**IdxVy**](#typedef-idxvy)  <br> |
| typedef Idx&lt; [**GridX**](structGridX.md) &gt; | [**IdxX**](#typedef-idxx)  <br> |
| typedef Idx&lt; [**GridX**](structGridX.md), [**GridY**](structGridY.md) &gt; | [**IdxXY**](#typedef-idxxy)  <br> |
| typedef Idx&lt; [**GridX**](structGridX.md), [**GridY**](structGridY.md), [**GridVx**](structGridVx.md), [**GridVy**](structGridVy.md) &gt; | [**IdxXYVxVy**](#typedef-idxxyvxvy)  <br> |
| typedef Idx&lt; [**GridY**](structGridY.md) &gt; | [**IdxY**](#typedef-idxy)  <br> |
| typedef ddc::GrevilleInterpolationPoints&lt; [**BSplinesVx**](structBSplinesVx.md), SplineVxBoundary, SplineVxBoundary &gt; | [**SplineInterpPointsVx**](#typedef-splineinterppointsvx)  <br> |
| typedef ddc::GrevilleInterpolationPoints&lt; [**BSplinesVy**](structBSplinesVy.md), SplineVyBoundary, SplineVyBoundary &gt; | [**SplineInterpPointsVy**](#typedef-splineinterppointsvy)  <br> |
| typedef ddc::GrevilleInterpolationPoints&lt; [**BSplinesX**](structBSplinesX.md), SplineXBoundary, SplineXBoundary &gt; | [**SplineInterpPointsX**](#typedef-splineinterppointsx)  <br> |
| typedef ddc::GrevilleInterpolationPoints&lt; [**BSplinesY**](structBSplinesY.md), SplineYBoundary, SplineYBoundary &gt; | [**SplineInterpPointsY**](#typedef-splineinterppointsy)  <br> |
| typedef ddc::SplineBuilder&lt; Kokkos::DefaultExecutionSpace, Kokkos::DefaultExecutionSpace::memory\_space, [**BSplinesVx**](structBSplinesVx.md), [**GridVx**](structGridVx.md), SplineVxBoundary, SplineVxBoundary, ddc::SplineSolver::LAPACK &gt; | [**SplineVxBuilder**](#typedef-splinevxbuilder)  <br> |
| typedef ddc::SplineEvaluator&lt; Kokkos::DefaultExecutionSpace, Kokkos::DefaultExecutionSpace::memory\_space, [**BSplinesVx**](structBSplinesVx.md), [**GridVx**](structGridVx.md), ddc::ConstantExtrapolationRule&lt; [**Vx**](structVx.md) &gt;, ddc::ConstantExtrapolationRule&lt; [**Vx**](structVx.md) &gt; &gt; | [**SplineVxEvaluator**](#typedef-splinevxevaluator)  <br> |
| typedef ddc::SplineBuilder&lt; Kokkos::DefaultExecutionSpace, Kokkos::DefaultExecutionSpace::memory\_space, [**BSplinesVy**](structBSplinesVy.md), [**GridVy**](structGridVy.md), SplineVyBoundary, SplineVyBoundary, ddc::SplineSolver::LAPACK &gt; | [**SplineVyBuilder**](#typedef-splinevybuilder)  <br> |
| typedef ddc::SplineEvaluator&lt; Kokkos::DefaultExecutionSpace, Kokkos::DefaultExecutionSpace::memory\_space, [**BSplinesVy**](structBSplinesVy.md), [**GridVy**](structGridVy.md), ddc::ConstantExtrapolationRule&lt; [**Vy**](structVy.md) &gt;, ddc::ConstantExtrapolationRule&lt; [**Vy**](structVy.md) &gt; &gt; | [**SplineVyEvaluator**](#typedef-splinevyevaluator)  <br> |
| typedef ddc::SplineBuilder&lt; Kokkos::DefaultExecutionSpace, Kokkos::DefaultExecutionSpace::memory\_space, [**BSplinesX**](structBSplinesX.md), [**GridX**](structGridX.md), SplineXBoundary, SplineXBoundary, ddc::SplineSolver::LAPACK &gt; | [**SplineXBuilder**](#typedef-splinexbuilder)  <br> |
| typedef ddc::SplineEvaluator&lt; Kokkos::DefaultExecutionSpace, Kokkos::DefaultExecutionSpace::memory\_space, [**BSplinesX**](structBSplinesX.md), [**GridX**](structGridX.md), ddc::PeriodicExtrapolationRule&lt; [**X**](structX.md) &gt;, ddc::PeriodicExtrapolationRule&lt; [**X**](structX.md) &gt; &gt; | [**SplineXEvaluator**](#typedef-splinexevaluator)  <br> |
| typedef ddc::SplineBuilder&lt; Kokkos::DefaultExecutionSpace, Kokkos::DefaultExecutionSpace::memory\_space, [**BSplinesY**](structBSplinesY.md), [**GridY**](structGridY.md), SplineYBoundary, SplineYBoundary, ddc::SplineSolver::LAPACK &gt; | [**SplineYBuilder**](#typedef-splineybuilder)  <br> |
| typedef ddc::SplineEvaluator&lt; Kokkos::DefaultExecutionSpace, Kokkos::DefaultExecutionSpace::memory\_space, [**BSplinesY**](structBSplinesY.md), [**GridY**](structGridY.md), ddc::PeriodicExtrapolationRule&lt; [**Y**](structY.md) &gt;, ddc::PeriodicExtrapolationRule&lt; [**Y**](structY.md) &gt; &gt; | [**SplineYEvaluator**](#typedef-splineyevaluator)  <br> |
| typedef [**MPILayout**](classMPILayout.md)&lt; IdxRangeSpVxVyXY, [**GridVx**](structGridVx.md), [**GridVy**](structGridVy.md) &gt; | [**V2DSplit**](#typedef-v2dsplit)  <br> |
| typedef [**MPILayout**](classMPILayout.md)&lt; IdxRangeSpXYVxVy, [**GridX**](structGridX.md), [**GridY**](structGridY.md) &gt; | [**X2DSplit**](#typedef-x2dsplit)  <br> |




## Public Attributes

| Type | Name |
| ---: | :--- |
|  int constexpr | [**BSDegreeVx**](#variable-bsdegreevx)   = `3`<br> |
|  int constexpr | [**BSDegreeVy**](#variable-bsdegreevy)   = `3`<br> |
|  int constexpr | [**BSDegreeX**](#variable-bsdegreex)   = `3`<br> |
|  int constexpr | [**BSDegreeY**](#variable-bsdegreey)   = `3`<br> |
|  bool constexpr | [**BsplineOnUniformCellsVx**](#variable-bsplineonuniformcellsvx)   = `true`<br> |
|  bool constexpr | [**BsplineOnUniformCellsVy**](#variable-bsplineonuniformcellsvy)   = `true`<br> |
|  bool constexpr | [**BsplineOnUniformCellsX**](#variable-bsplineonuniformcellsx)   = `true`<br> |
|  bool constexpr | [**BsplineOnUniformCellsY**](#variable-bsplineonuniformcellsy)   = `true`<br> |
|  ddc::BoundCond constexpr | [**SplineVxBoundary**](#variable-splinevxboundary)   = `ddc::BoundCond::HERMITE`<br> |
|  ddc::BoundCond constexpr | [**SplineVyBoundary**](#variable-splinevyboundary)   = `ddc::BoundCond::HERMITE`<br> |
|  ddc::BoundCond constexpr | [**SplineXBoundary**](#variable-splinexboundary)   = `ddc::BoundCond::PERIODIC`<br> |
|  ddc::BoundCond constexpr | [**SplineYBoundary**](#variable-splineyboundary)   = `ddc::BoundCond::PERIODIC`<br> |












































## Public Types Documentation




### typedef BSConstFieldXY 

```C++
using BSConstFieldXY =  Field<ElementType const, IdxRangeBSXY>;
```




<hr>



### typedef ConstFieldSpVxVy 

```C++
using ConstFieldSpVxVy =  Field<ElementType const, IdxRangeSpVxVy>;
```




<hr>



### typedef ConstFieldSpVxVyXY 

```C++
using ConstFieldSpVxVyXY =  Field<ElementType const, IdxRangeSpVxVyXY>;
```




<hr>



### typedef ConstFieldSpXYVxVy 

```C++
using ConstFieldSpXYVxVy =  Field<ElementType const, IdxRangeSpXYVxVy>;
```




<hr>



### typedef ConstFieldVx 

```C++
using ConstFieldVx =  Field<ElementType const, IdxRangeVx>;
```




<hr>



### typedef ConstFieldVxVy 

```C++
using ConstFieldVxVy =  Field<ElementType const, IdxRangeVxVy>;
```




<hr>



### typedef ConstFieldVy 

```C++
using ConstFieldVy =  Field<ElementType const, IdxRangeVy>;
```




<hr>



### typedef ConstFieldX 

```C++
using ConstFieldX =  Field<ElementType const, IdxRangeX>;
```




<hr>



### typedef ConstFieldXY 

```C++
using ConstFieldXY =  Field<ElementType const, IdxRangeXY>;
```




<hr>



### typedef ConstFieldY 

```C++
using ConstFieldY =  Field<ElementType const, IdxRangeY>;
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



### typedef DBSConstFieldXY 

```C++
using DBSConstFieldXY =  BSConstFieldXY<double>;
```




<hr>



### typedef DConstFieldSpVxVy 

```C++
using DConstFieldSpVxVy =  ConstFieldSpVxVy<double>;
```




<hr>



### typedef DConstFieldSpVxVyXY 

```C++
using DConstFieldSpVxVyXY =  ConstFieldSpVxVyXY<double>;
```




<hr>



### typedef DConstFieldSpXYVxVy 

```C++
using DConstFieldSpXYVxVy =  ConstFieldSpXYVxVy<double>;
```




<hr>



### typedef DConstFieldVxVy 

```C++
using DConstFieldVxVy =  ConstFieldVxVy<double>;
```




<hr>



### typedef DConstFieldXY 

```C++
using DConstFieldXY =  ConstFieldXY<double>;
```




<hr>



### typedef DFieldMemSpVxVy 

```C++
using DFieldMemSpVxVy =  FieldMemSpVxVy<double>;
```




<hr>



### typedef DFieldMemSpVxVyXY 

```C++
using DFieldMemSpVxVyXY =  FieldMemSpVxVyXY<double>;
```




<hr>



### typedef DFieldMemSpXYVxVy 

```C++
using DFieldMemSpXYVxVy =  FieldMemSpXYVxVy<double>;
```




<hr>



### typedef DFieldMemVxVy 

```C++
using DFieldMemVxVy =  FieldMemVxVy<double>;
```




<hr>



### typedef DFieldMemX 

```C++
using DFieldMemX =  FieldMemX<double>;
```




<hr>



### typedef DFieldMemXY 

```C++
using DFieldMemXY =  FieldMemXY<double>;
```




<hr>



### typedef DFieldMemXYVxVy 

```C++
using DFieldMemXYVxVy =  FieldMemXYVxVy<double>;
```




<hr>



### typedef DFieldMemY 

```C++
using DFieldMemY =  FieldMemY<double>;
```




<hr>



### typedef DFieldSpVxVy 

```C++
using DFieldSpVxVy =  FieldSpVxVy<double>;
```




<hr>



### typedef DFieldSpVxVyXY 

```C++
using DFieldSpVxVyXY =  FieldSpVxVyXY<double>;
```




<hr>



### typedef DFieldSpXYVxVy 

```C++
using DFieldSpXYVxVy =  FieldSpXYVxVy<double>;
```




<hr>



### typedef DFieldVx 

```C++
using DFieldVx =  FieldVx<double>;
```




<hr>



### typedef DFieldVxVy 

```C++
using DFieldVxVy =  FieldVxVy<double>;
```




<hr>



### typedef DFieldVy 

```C++
using DFieldVy =  FieldVy<double>;
```




<hr>



### typedef DFieldX 

```C++
using DFieldX =  FieldX<double>;
```




<hr>



### typedef DFieldXY 

```C++
using DFieldXY =  FieldXY<double>;
```




<hr>



### typedef DFieldY 

```C++
using DFieldY =  FieldY<double>;
```




<hr>



### typedef FieldMemSpVxVy 

```C++
using FieldMemSpVxVy =  FieldMem<ElementType, IdxRangeSpVxVy>;
```




<hr>



### typedef FieldMemSpVxVyXY 

```C++
using FieldMemSpVxVyXY =  FieldMem<ElementType, IdxRangeSpVxVyXY>;
```




<hr>



### typedef FieldMemSpXYVxVy 

```C++
using FieldMemSpXYVxVy =  FieldMem<ElementType, IdxRangeSpXYVxVy>;
```




<hr>



### typedef FieldMemVx 

```C++
using FieldMemVx =  FieldMem<ElementType, IdxRangeVx>;
```




<hr>



### typedef FieldMemVxVy 

```C++
using FieldMemVxVy =  FieldMem<ElementType, IdxRangeVxVy>;
```




<hr>



### typedef FieldMemVy 

```C++
using FieldMemVy =  FieldMem<ElementType, IdxRangeVy>;
```




<hr>



### typedef FieldMemX 

```C++
using FieldMemX =  FieldMem<ElementType, IdxRangeX>;
```




<hr>



### typedef FieldMemXY 

```C++
using FieldMemXY =  FieldMem<ElementType, IdxRangeXY>;
```




<hr>



### typedef FieldMemXYVxVy 

```C++
using FieldMemXYVxVy =  FieldMem<ElementType, IdxRangeXYVxVy>;
```




<hr>



### typedef FieldMemY 

```C++
using FieldMemY =  FieldMem<ElementType, IdxRangeY>;
```




<hr>



### typedef FieldSpVxVy 

```C++
using FieldSpVxVy =  Field<ElementType, IdxRangeSpVxVy>;
```




<hr>



### typedef FieldSpVxVyXY 

```C++
using FieldSpVxVyXY =  Field<ElementType, IdxRangeSpVxVyXY>;
```




<hr>



### typedef FieldSpXYVxVy 

```C++
using FieldSpXYVxVy =  Field<ElementType, IdxRangeSpXYVxVy>;
```




<hr>



### typedef FieldVx 

```C++
using FieldVx =  Field<ElementType, IdxRangeVx>;
```




<hr>



### typedef FieldVxVy 

```C++
using FieldVxVy =  Field<ElementType, IdxRangeVxVy>;
```




<hr>



### typedef FieldVy 

```C++
using FieldVy =  Field<ElementType, IdxRangeVy>;
```




<hr>



### typedef FieldX 

```C++
using FieldX =  Field<ElementType, IdxRangeX>;
```




<hr>



### typedef FieldXY 

```C++
using FieldXY =  Field<ElementType, IdxRangeXY>;
```




<hr>



### typedef FieldY 

```C++
using FieldY =  Field<ElementType, IdxRangeY>;
```




<hr>



### typedef IdxRangeBSVx 

```C++
using IdxRangeBSVx =  IdxRange<BSplinesVx>;
```




<hr>



### typedef IdxRangeBSVxVy 

```C++
using IdxRangeBSVxVy =  IdxRange<BSplinesVx, BSplinesVy>;
```




<hr>



### typedef IdxRangeBSVy 

```C++
using IdxRangeBSVy =  IdxRange<BSplinesVy>;
```




<hr>



### typedef IdxRangeBSX 

```C++
using IdxRangeBSX =  IdxRange<BSplinesX>;
```




<hr>



### typedef IdxRangeBSXY 

```C++
using IdxRangeBSXY =  IdxRange<BSplinesX, BSplinesY>;
```




<hr>



### typedef IdxRangeBSY 

```C++
using IdxRangeBSY =  IdxRange<BSplinesY>;
```




<hr>



### typedef IdxRangeSpVxVy 

```C++
using IdxRangeSpVxVy =  IdxRange<Species, GridVx, GridVy>;
```




<hr>



### typedef IdxRangeSpVxVyXY 

```C++
using IdxRangeSpVxVyXY =  IdxRange<Species, GridVx, GridVy, GridX, GridY>;
```




<hr>



### typedef IdxRangeSpXYVxVy 

```C++
using IdxRangeSpXYVxVy =  IdxRange<Species, GridX, GridY, GridVx, GridVy>;
```




<hr>



### typedef IdxRangeVx 

```C++
using IdxRangeVx =  IdxRange<GridVx>;
```




<hr>



### typedef IdxRangeVxVy 

```C++
using IdxRangeVxVy =  IdxRange<GridVx, GridVy>;
```




<hr>



### typedef IdxRangeVxVyXY 

```C++
using IdxRangeVxVyXY =  IdxRange<GridVx, GridVy, GridX, GridY>;
```




<hr>



### typedef IdxRangeVy 

```C++
using IdxRangeVy =  IdxRange<GridVy>;
```




<hr>



### typedef IdxRangeX 

```C++
using IdxRangeX =  IdxRange<GridX>;
```




<hr>



### typedef IdxRangeXY 

```C++
using IdxRangeXY =  IdxRange<GridX, GridY>;
```




<hr>



### typedef IdxRangeXYVxVy 

```C++
using IdxRangeXYVxVy =  IdxRange<GridX, GridY, GridVx, GridVy>;
```




<hr>



### typedef IdxRangeY 

```C++
using IdxRangeY =  IdxRange<GridY>;
```




<hr>



### typedef IdxSpXYVxVy 

```C++
using IdxSpXYVxVy =  Idx<Species, GridX, GridY, GridVx, GridVy>;
```




<hr>



### typedef IdxStepVx 

```C++
using IdxStepVx =  IdxStep<GridVx>;
```




<hr>



### typedef IdxStepVy 

```C++
using IdxStepVy =  IdxStep<GridVy>;
```




<hr>



### typedef IdxStepX 

```C++
using IdxStepX =  IdxStep<GridX>;
```




<hr>



### typedef IdxStepY 

```C++
using IdxStepY =  IdxStep<GridY>;
```




<hr>



### typedef IdxVx 

```C++
using IdxVx =  Idx<GridVx>;
```




<hr>



### typedef IdxVxVy 

```C++
using IdxVxVy =  Idx<GridVx, GridVy>;
```




<hr>



### typedef IdxVy 

```C++
using IdxVy =  Idx<GridVy>;
```




<hr>



### typedef IdxX 

```C++
using IdxX =  Idx<GridX>;
```




<hr>



### typedef IdxXY 

```C++
using IdxXY =  Idx<GridX, GridY>;
```




<hr>



### typedef IdxXYVxVy 

```C++
using IdxXYVxVy =  Idx<GridX, GridY, GridVx, GridVy>;
```




<hr>



### typedef IdxY 

```C++
using IdxY =  Idx<GridY>;
```




<hr>



### typedef SplineInterpPointsVx 

```C++
using SplineInterpPointsVx =  ddc::GrevilleInterpolationPoints<BSplinesVx, SplineVxBoundary, SplineVxBoundary>;
```




<hr>



### typedef SplineInterpPointsVy 

```C++
using SplineInterpPointsVy =  ddc::GrevilleInterpolationPoints<BSplinesVy, SplineVyBoundary, SplineVyBoundary>;
```




<hr>



### typedef SplineInterpPointsX 

```C++
using SplineInterpPointsX =  ddc::GrevilleInterpolationPoints<BSplinesX, SplineXBoundary, SplineXBoundary>;
```




<hr>



### typedef SplineInterpPointsY 

```C++
using SplineInterpPointsY =  ddc::GrevilleInterpolationPoints<BSplinesY, SplineYBoundary, SplineYBoundary>;
```




<hr>



### typedef SplineVxBuilder 

```C++
using SplineVxBuilder =  ddc::SplineBuilder< Kokkos::DefaultExecutionSpace, Kokkos::DefaultExecutionSpace::memory_space, BSplinesVx, GridVx, SplineVxBoundary, SplineVxBoundary, ddc::SplineSolver::LAPACK>;
```




<hr>



### typedef SplineVxEvaluator 

```C++
using SplineVxEvaluator =  ddc::SplineEvaluator< Kokkos::DefaultExecutionSpace, Kokkos::DefaultExecutionSpace::memory_space, BSplinesVx, GridVx, ddc::ConstantExtrapolationRule<Vx>, ddc::ConstantExtrapolationRule<Vx> >;
```




<hr>



### typedef SplineVyBuilder 

```C++
using SplineVyBuilder =  ddc::SplineBuilder< Kokkos::DefaultExecutionSpace, Kokkos::DefaultExecutionSpace::memory_space, BSplinesVy, GridVy, SplineVyBoundary, SplineVyBoundary, ddc::SplineSolver::LAPACK>;
```




<hr>



### typedef SplineVyEvaluator 

```C++
using SplineVyEvaluator =  ddc::SplineEvaluator< Kokkos::DefaultExecutionSpace, Kokkos::DefaultExecutionSpace::memory_space, BSplinesVy, GridVy, ddc::ConstantExtrapolationRule<Vy>, ddc::ConstantExtrapolationRule<Vy> >;
```




<hr>



### typedef SplineXBuilder 

```C++
using SplineXBuilder =  ddc::SplineBuilder< Kokkos::DefaultExecutionSpace, Kokkos::DefaultExecutionSpace::memory_space, BSplinesX, GridX, SplineXBoundary, SplineXBoundary, ddc::SplineSolver::LAPACK>;
```




<hr>



### typedef SplineXEvaluator 

```C++
using SplineXEvaluator =  ddc::SplineEvaluator< Kokkos::DefaultExecutionSpace, Kokkos::DefaultExecutionSpace::memory_space, BSplinesX, GridX, ddc::PeriodicExtrapolationRule<X>, ddc::PeriodicExtrapolationRule<X> >;
```




<hr>



### typedef SplineYBuilder 

```C++
using SplineYBuilder =  ddc::SplineBuilder< Kokkos::DefaultExecutionSpace, Kokkos::DefaultExecutionSpace::memory_space, BSplinesY, GridY, SplineYBoundary, SplineYBoundary, ddc::SplineSolver::LAPACK>;
```




<hr>



### typedef SplineYEvaluator 

```C++
using SplineYEvaluator =  ddc::SplineEvaluator< Kokkos::DefaultExecutionSpace, Kokkos::DefaultExecutionSpace::memory_space, BSplinesY, GridY, ddc::PeriodicExtrapolationRule<Y>, ddc::PeriodicExtrapolationRule<Y> >;
```




<hr>



### typedef V2DSplit 

```C++
using V2DSplit =  MPILayout<IdxRangeSpVxVyXY, GridVx, GridVy>;
```




<hr>



### typedef X2DSplit 

```C++
using X2DSplit =  MPILayout<IdxRangeSpXYVxVy, GridX, GridY>;
```




<hr>
## Public Attributes Documentation




### variable BSDegreeVx 

```C++
int constexpr BSDegreeVx;
```




<hr>



### variable BSDegreeVy 

```C++
int constexpr BSDegreeVy;
```




<hr>



### variable BSDegreeX 

```C++
int constexpr BSDegreeX;
```




<hr>



### variable BSDegreeY 

```C++
int constexpr BSDegreeY;
```




<hr>



### variable BsplineOnUniformCellsVx 

```C++
bool constexpr BsplineOnUniformCellsVx;
```




<hr>



### variable BsplineOnUniformCellsVy 

```C++
bool constexpr BsplineOnUniformCellsVy;
```




<hr>



### variable BsplineOnUniformCellsX 

```C++
bool constexpr BsplineOnUniformCellsX;
```




<hr>



### variable BsplineOnUniformCellsY 

```C++
bool constexpr BsplineOnUniformCellsY;
```




<hr>



### variable SplineVxBoundary 

```C++
ddc::BoundCond constexpr SplineVxBoundary;
```




<hr>



### variable SplineVyBoundary 

```C++
ddc::BoundCond constexpr SplineVyBoundary;
```




<hr>



### variable SplineXBoundary 

```C++
ddc::BoundCond constexpr SplineXBoundary;
```




<hr>



### variable SplineYBoundary 

```C++
ddc::BoundCond constexpr SplineYBoundary;
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/geometryXYVxVy/geometry/geometry.hpp`

