

# File geometry.hpp



[**FileList**](files.md) **>** [**geometry**](dir_6ef3b5c953c12640e6eb10271de0236d.md) **>** [**geometry.hpp**](geometryXY_2geometry_2geometry_8hpp.md)

[Go to the source code of this file](geometryXY_2geometry_2geometry_8hpp_source.md)



* `#include <ddc/ddc.hpp>`
* `#include "ddc_alias_inline_functions.hpp"`
* `#include "ddc_aliases.hpp"`
* `#include "ddc_helper.hpp"`
* `#include "vector_field.hpp"`
* `#include "vector_field_mem.hpp"`
* `#include "vector_index_tools.hpp"`















## Classes

| Type | Name |
| ---: | :--- |
| struct | [**GridX**](structGridX.md) <br> |
| struct | [**GridY**](structGridY.md) <br> |
| struct | [**X**](structX.md) <br>_Define non periodic real_ [_**X**_](structX.md) _dimension._ |
| struct | [**Y**](structY.md) <br>_Define non periodic real_ [_**Y**_](structY.md) _dimension._ |


## Public Types

| Type | Name |
| ---: | :--- |
| typedef Field&lt; ElementType const, IdxRangeX &gt; | [**ConstFieldX**](#typedef-constfieldx)  <br> |
| typedef Field&lt; ElementType const, IdxRangeXY &gt; | [**ConstFieldXY**](#typedef-constfieldxy)  <br> |
| typedef Field&lt; ElementType const, IdxRangeY &gt; | [**ConstFieldY**](#typedef-constfieldy)  <br> |
| typedef Coord&lt; [**X**](structX.md) &gt; | [**CoordX**](#typedef-coordx)  <br> |
| typedef Coord&lt; [**X**](structX.md), [**Y**](structY.md) &gt; | [**CoordXY**](#typedef-coordxy)  <br> |
| typedef Coord&lt; [**Y**](structY.md) &gt; | [**CoordY**](#typedef-coordy)  <br> |
| typedef ConstFieldXY&lt; double &gt; | [**DConstFieldXY**](#typedef-dconstfieldxy)  <br> |
| typedef FieldMemX&lt; double &gt; | [**DFieldMemX**](#typedef-dfieldmemx)  <br> |
| typedef FieldMemXY&lt; double &gt; | [**DFieldMemXY**](#typedef-dfieldmemxy)  <br> |
| typedef FieldMemY&lt; double &gt; | [**DFieldMemY**](#typedef-dfieldmemy)  <br> |
| typedef FieldX&lt; double &gt; | [**DFieldX**](#typedef-dfieldx)  <br> |
| typedef FieldXY&lt; double &gt; | [**DFieldXY**](#typedef-dfieldxy)  <br> |
| typedef FieldY&lt; double &gt; | [**DFieldY**](#typedef-dfieldy)  <br> |
| typedef FieldMem&lt; ElementType, IdxRangeX &gt; | [**FieldMemX**](#typedef-fieldmemx)  <br> |
| typedef FieldMem&lt; ElementType, IdxRangeXY &gt; | [**FieldMemXY**](#typedef-fieldmemxy)  <br> |
| typedef FieldMem&lt; ElementType, IdxRangeY &gt; | [**FieldMemY**](#typedef-fieldmemy)  <br> |
| typedef Field&lt; ElementType, IdxRangeX &gt; | [**FieldX**](#typedef-fieldx)  <br> |
| typedef Field&lt; ElementType, IdxRangeXY &gt; | [**FieldXY**](#typedef-fieldxy)  <br> |
| typedef Field&lt; ElementType, IdxRangeY &gt; | [**FieldY**](#typedef-fieldy)  <br> |
| typedef IdxRange&lt; [**GridX**](structGridX.md) &gt; | [**IdxRangeX**](#typedef-idxrangex)  <br> |
| typedef IdxRange&lt; [**GridX**](structGridX.md), [**GridY**](structGridY.md) &gt; | [**IdxRangeXY**](#typedef-idxrangexy)  <br> |
| typedef IdxRange&lt; [**GridY**](structGridY.md) &gt; | [**IdxRangeY**](#typedef-idxrangey)  <br> |
| typedef IdxStep&lt; [**GridX**](structGridX.md) &gt; | [**IdxStepX**](#typedef-idxstepx)  <br> |
| typedef IdxStep&lt; [**GridY**](structGridY.md) &gt; | [**IdxStepY**](#typedef-idxstepy)  <br> |
| typedef Idx&lt; [**GridX**](structGridX.md) &gt; | [**IdxX**](#typedef-idxx)  <br> |
| typedef Idx&lt; [**GridX**](structGridX.md), [**GridY**](structGridY.md) &gt; | [**IdxXY**](#typedef-idxxy)  <br> |
| typedef Idx&lt; [**GridY**](structGridY.md) &gt; | [**IdxY**](#typedef-idxy)  <br> |
| typedef typename [**VectorFieldMemXY\_XY::view\_type**](classVectorFieldMem.md#typedef-view_type) | [**VectorConstFieldXY\_XY**](#typedef-vectorconstfieldxy_xy)  <br> |
| typedef [**VectorFieldMem**](classVectorFieldMem.md)&lt; double, IdxRangeXY, VectorIndexSet&lt; [**X**](structX.md), [**Y**](structY.md) &gt;, Kokkos::DefaultExecutionSpace::memory\_space &gt; | [**VectorFieldMemXY\_XY**](#typedef-vectorfieldmemxy_xy)  <br> |
| typedef typename [**VectorFieldMemXY\_XY::span\_type**](classVectorFieldMem.md#typedef-span_type) | [**VectorFieldXY\_XY**](#typedef-vectorfieldxy_xy)  <br> |
















































## Public Types Documentation




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



### typedef DConstFieldXY 

```C++
using DConstFieldXY =  ConstFieldXY<double>;
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



### typedef DFieldMemY 

```C++
using DFieldMemY =  FieldMemY<double>;
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



### typedef FieldMemY 

```C++
using FieldMemY =  FieldMem<ElementType, IdxRangeY>;
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



### typedef IdxRangeY 

```C++
using IdxRangeY =  IdxRange<GridY>;
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



### typedef IdxY 

```C++
using IdxY =  Idx<GridY>;
```




<hr>



### typedef VectorConstFieldXY\_XY 

```C++
using VectorConstFieldXY_XY =  typename VectorFieldMemXY_XY::view_type;
```




<hr>



### typedef VectorFieldMemXY\_XY 

```C++
using VectorFieldMemXY_XY =  VectorFieldMem< double, IdxRangeXY, VectorIndexSet<X, Y>, Kokkos::DefaultExecutionSpace::memory_space>;
```




<hr>



### typedef VectorFieldXY\_XY 

```C++
using VectorFieldXY_XY =  typename VectorFieldMemXY_XY::span_type;
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/geometryXY/geometry/geometry.hpp`

