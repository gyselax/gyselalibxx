

# File geometry\_r\_theta.hpp



[**FileList**](files.md) **>** [**geometry**](dir_718520565cc7a7cfd9ba0e7c9c4c6d52.md) **>** [**geometry\_r\_theta.hpp**](geometry__r__theta_8hpp.md)

[Go to the source code of this file](geometry__r__theta_8hpp_source.md)



* `#include <ddc/ddc.hpp>`
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
| struct | [**GridR**](structGridR.md) <br> |
| struct | [**GridTheta**](structGridTheta.md) <br> |
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
| typedef FieldMemR&lt; double &gt; | [**DFieldMemR**](#typedef-dfieldmemr)  <br> |
| typedef FieldMemRTheta&lt; double &gt; | [**DFieldMemRTheta**](#typedef-dfieldmemrtheta)  <br> |
| typedef FieldMemTheta&lt; double &gt; | [**DFieldMemTheta**](#typedef-dfieldmemtheta)  <br> |
| typedef FieldR&lt; double &gt; | [**DFieldR**](#typedef-dfieldr)  <br> |
| typedef FieldRTheta&lt; double &gt; | [**DFieldRTheta**](#typedef-dfieldrtheta)  <br> |
| typedef FieldTheta&lt; double &gt; | [**DFieldTheta**](#typedef-dfieldtheta)  <br> |
| typedef [**VectorConstField**](classVectorField.md)&lt; double, IdxRangeRTheta, VectorIndexSet&lt; Dim1, Dim2 &gt; &gt; | [**DVectorConstFieldRTheta**](#typedef-dvectorconstfieldrtheta)  <br> |
| typedef [**VectorFieldMem**](classVectorFieldMem.md)&lt; double, IdxRangeRTheta, VectorIndexSet&lt; Dim1, Dim2 &gt; &gt; | [**DVectorFieldMemRTheta**](#typedef-dvectorfieldmemrtheta)  <br> |
| typedef [**VectorField**](classVectorField.md)&lt; double, IdxRangeRTheta, VectorIndexSet&lt; Dim1, Dim2 &gt; &gt; | [**DVectorFieldRTheta**](#typedef-dvectorfieldrtheta)  <br> |
| typedef FieldMem&lt; ElementType, IdxRangeR &gt; | [**FieldMemR**](#typedef-fieldmemr)  <br> |
| typedef FieldMem&lt; ElementType, IdxRangeRTheta &gt; | [**FieldMemRTheta**](#typedef-fieldmemrtheta)  <br> |
| typedef FieldMem&lt; ElementType, IdxRangeTheta &gt; | [**FieldMemTheta**](#typedef-fieldmemtheta)  <br> |
| typedef Field&lt; ElementType, IdxRangeR &gt; | [**FieldR**](#typedef-fieldr)  <br> |
| typedef Field&lt; ElementType, IdxRangeRTheta &gt; | [**FieldRTheta**](#typedef-fieldrtheta)  <br> |
| typedef Field&lt; ElementType, IdxRangeTheta &gt; | [**FieldTheta**](#typedef-fieldtheta)  <br> |
| typedef Idx&lt; [**GridR**](structGridR.md) &gt; | [**IdxR**](#typedef-idxr)  <br> |
| typedef Idx&lt; [**GridR**](structGridR.md), [**GridTheta**](structGridTheta.md) &gt; | [**IdxRTheta**](#typedef-idxrtheta)  <br> |
| typedef IdxRange&lt; [**GridR**](structGridR.md) &gt; | [**IdxRangeR**](#typedef-idxranger)  <br> |
| typedef IdxRange&lt; [**GridR**](structGridR.md), [**GridTheta**](structGridTheta.md) &gt; | [**IdxRangeRTheta**](#typedef-idxrangertheta)  <br> |
| typedef IdxRange&lt; [**GridTheta**](structGridTheta.md) &gt; | [**IdxRangeTheta**](#typedef-idxrangetheta)  <br> |
| typedef IdxStep&lt; [**GridR**](structGridR.md) &gt; | [**IdxStepR**](#typedef-idxstepr)  <br> |
| typedef IdxStep&lt; [**GridR**](structGridR.md), [**GridTheta**](structGridTheta.md) &gt; | [**IdxStepRTheta**](#typedef-idxsteprtheta)  <br> |
| typedef IdxStep&lt; [**GridTheta**](structGridTheta.md) &gt; | [**IdxStepTheta**](#typedef-idxsteptheta)  <br> |
| typedef Idx&lt; [**GridTheta**](structGridTheta.md) &gt; | [**IdxTheta**](#typedef-idxtheta)  <br> |
















































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



### typedef DVectorConstFieldRTheta 

```C++
using DVectorConstFieldRTheta =  VectorConstField<double, IdxRangeRTheta, VectorIndexSet<Dim1, Dim2> >;
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

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/geometryRTheta/geometry/geometry_r_theta.hpp`

