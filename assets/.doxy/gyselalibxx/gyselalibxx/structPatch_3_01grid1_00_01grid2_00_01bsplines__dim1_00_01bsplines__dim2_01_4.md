

# Struct Patch&lt; grid1, grid2, bsplines\_dim1, bsplines\_dim2 &gt;

**template &lt;class grid1, class grid2, class bsplines\_dim1, class bsplines\_dim2&gt;**



[**ClassList**](annotated.md) **>** [**Patch&lt; grid1, grid2, bsplines\_dim1, bsplines\_dim2 &gt;**](structPatch_3_01grid1_00_01grid2_00_01bsplines__dim1_00_01bsplines__dim2_01_4.md)



_Tag for a patch._ [More...](#detailed-description)

* `#include <patch.hpp>`

















## Public Types

| Type | Name |
| ---: | :--- |
| typedef bsplines\_dim1 | [**BSplines1**](#typedef-bsplines1)  <br>_B-splines defined along the first dimension._  |
| typedef bsplines\_dim2 | [**BSplines2**](#typedef-bsplines2)  <br>_B-splines defined along the second dimension._  |
| typedef Coord&lt; [**Dim1**](structPatch_3_01grid1_00_01grid2_00_01bsplines__dim1_00_01bsplines__dim2_01_4.md#typedef-dim1) &gt; | [**Coord1**](#typedef-coord1)  <br>_Coordinate type on the first continuous dimension._  |
| typedef Coord&lt; [**Dim1**](structPatch_3_01grid1_00_01grid2_00_01bsplines__dim1_00_01bsplines__dim2_01_4.md#typedef-dim1), [**Dim2**](structPatch_3_01grid1_00_01grid2_00_01bsplines__dim1_00_01bsplines__dim2_01_4.md#typedef-dim2) &gt; | [**Coord12**](#typedef-coord12)  <br>_Coordinate type on the 2D continuous dimensions._  |
| typedef Coord&lt; [**Dim2**](structPatch_3_01grid1_00_01grid2_00_01bsplines__dim1_00_01bsplines__dim2_01_4.md#typedef-dim2) &gt; | [**Coord2**](#typedef-coord2)  <br>_Coordinate type on the second continuous dimension._  |
| typedef typename grid1::continuous\_dimension\_type | [**Dim1**](#typedef-dim1)  <br>_First continuous dimension._  |
| typedef typename grid2::continuous\_dimension\_type | [**Dim2**](#typedef-dim2)  <br>_Second continuous dimension._  |
| typedef grid1 | [**Grid1**](#typedef-grid1)  <br>_Grid on the first logical dimension._  |
| typedef grid2 | [**Grid2**](#typedef-grid2)  <br>_Grid on the second logical dimension._  |
| typedef Idx&lt; [**Grid1**](structPatch_3_01grid1_00_01grid2_00_01bsplines__dim1_00_01bsplines__dim2_01_4.md#typedef-grid1) &gt; | [**Idx1**](#typedef-idx1)  <br>_1D index of a grid point along the first dimension._  |
| typedef Idx&lt; [**Grid1**](structPatch_3_01grid1_00_01grid2_00_01bsplines__dim1_00_01bsplines__dim2_01_4.md#typedef-grid1), [**Grid2**](structPatch_3_01grid1_00_01grid2_00_01bsplines__dim1_00_01bsplines__dim2_01_4.md#typedef-grid2) &gt; | [**Idx12**](#typedef-idx12)  <br>_2D index of a grid point along the first and second dimensions._  |
| typedef Idx&lt; [**Grid2**](structPatch_3_01grid1_00_01grid2_00_01bsplines__dim1_00_01bsplines__dim2_01_4.md#typedef-grid2) &gt; | [**Idx2**](#typedef-idx2)  <br>_1D index of a grid point along the second dimension._  |
| typedef IdxRange&lt; [**Grid1**](structPatch_3_01grid1_00_01grid2_00_01bsplines__dim1_00_01bsplines__dim2_01_4.md#typedef-grid1) &gt; | [**IdxRange1**](#typedef-idxrange1)  <br>_Index range of a grids over the first dimension._  |
| typedef IdxRange&lt; [**Grid1**](structPatch_3_01grid1_00_01grid2_00_01bsplines__dim1_00_01bsplines__dim2_01_4.md#typedef-grid1), [**Grid2**](structPatch_3_01grid1_00_01grid2_00_01bsplines__dim1_00_01bsplines__dim2_01_4.md#typedef-grid2) &gt; | [**IdxRange12**](#typedef-idxrange12)  <br>_Index range of a grids over the first and second dimension._  |
| typedef IdxRange&lt; [**Grid2**](structPatch_3_01grid1_00_01grid2_00_01bsplines__dim1_00_01bsplines__dim2_01_4.md#typedef-grid2) &gt; | [**IdxRange2**](#typedef-idxrange2)  <br>_Index range of a grids over the second dimension._  |
| typedef IdxRange&lt; [**BSplines1**](structPatch_3_01grid1_00_01grid2_00_01bsplines__dim1_00_01bsplines__dim2_01_4.md#typedef-bsplines1) &gt; | [**IdxRangeBS1**](#typedef-idxrangebs1)  <br>_Index range of a grids over the first spline dimension._  |
| typedef IdxRange&lt; [**BSplines1**](structPatch_3_01grid1_00_01grid2_00_01bsplines__dim1_00_01bsplines__dim2_01_4.md#typedef-bsplines1), [**BSplines2**](structPatch_3_01grid1_00_01grid2_00_01bsplines__dim1_00_01bsplines__dim2_01_4.md#typedef-bsplines2) &gt; | [**IdxRangeBS12**](#typedef-idxrangebs12)  <br>_Index range of a grids over the first and second spline dimension._  |
| typedef IdxRange&lt; [**BSplines2**](structPatch_3_01grid1_00_01grid2_00_01bsplines__dim1_00_01bsplines__dim2_01_4.md#typedef-bsplines2) &gt; | [**IdxRangeBS2**](#typedef-idxrangebs2)  <br>_Index range of a grids over the second spline dimension._  |
| typedef IdxStep&lt; [**Grid1**](structPatch_3_01grid1_00_01grid2_00_01bsplines__dim1_00_01bsplines__dim2_01_4.md#typedef-grid1) &gt; | [**IdxStep1**](#typedef-idxstep1)  <br>_1D index step between grid points along the first dimension._  |
| typedef IdxStep&lt; [**Grid1**](structPatch_3_01grid1_00_01grid2_00_01bsplines__dim1_00_01bsplines__dim2_01_4.md#typedef-grid1), [**Grid2**](structPatch_3_01grid1_00_01grid2_00_01bsplines__dim1_00_01bsplines__dim2_01_4.md#typedef-grid2) &gt; | [**IdxStep12**](#typedef-idxstep12)  <br>_2D index step between grid points along the first and second dimensions._  |
| typedef IdxStep&lt; [**Grid2**](structPatch_3_01grid1_00_01grid2_00_01bsplines__dim1_00_01bsplines__dim2_01_4.md#typedef-grid2) &gt; | [**IdxStep2**](#typedef-idxstep2)  <br>_1D index step between grid points along the second dimension._  |






## Public Static Attributes

| Type | Name |
| ---: | :--- |
|  int constexpr | [**n\_dims**](#variable-n_dims)   = `2`<br>_The number of dimensions of the patch._  |










































## Detailed Description




**Template parameters:**


* `grid1` Grid on the first logical dimension associated to the patch. 
* `grid2` Grid on the second logical dimension associated to the patch. 
* `bsplines_dim1` Bspline dimension defined on the first logical continuous dimension associated to the patch. 
* `bsplines_dim2` Bspline dimension defined on the second logical continuous dimension associated to the patch. 




    
## Public Types Documentation




### typedef BSplines1 

_B-splines defined along the first dimension._ 
```C++
using Patch< grid1, grid2, bsplines_dim1, bsplines_dim2 >::BSplines1 =  bsplines_dim1;
```




<hr>



### typedef BSplines2 

_B-splines defined along the second dimension._ 
```C++
using Patch< grid1, grid2, bsplines_dim1, bsplines_dim2 >::BSplines2 =  bsplines_dim2;
```




<hr>



### typedef Coord1 

_Coordinate type on the first continuous dimension._ 
```C++
using Patch< grid1, grid2, bsplines_dim1, bsplines_dim2 >::Coord1 =  Coord<Dim1>;
```




<hr>



### typedef Coord12 

_Coordinate type on the 2D continuous dimensions._ 
```C++
using Patch< grid1, grid2, bsplines_dim1, bsplines_dim2 >::Coord12 =  Coord<Dim1, Dim2>;
```




<hr>



### typedef Coord2 

_Coordinate type on the second continuous dimension._ 
```C++
using Patch< grid1, grid2, bsplines_dim1, bsplines_dim2 >::Coord2 =  Coord<Dim2>;
```




<hr>



### typedef Dim1 

_First continuous dimension._ 
```C++
using Patch< grid1, grid2, bsplines_dim1, bsplines_dim2 >::Dim1 =  typename grid1::continuous_dimension_type;
```




<hr>



### typedef Dim2 

_Second continuous dimension._ 
```C++
using Patch< grid1, grid2, bsplines_dim1, bsplines_dim2 >::Dim2 =  typename grid2::continuous_dimension_type;
```




<hr>



### typedef Grid1 

_Grid on the first logical dimension._ 
```C++
using Patch< grid1, grid2, bsplines_dim1, bsplines_dim2 >::Grid1 =  grid1;
```




<hr>



### typedef Grid2 

_Grid on the second logical dimension._ 
```C++
using Patch< grid1, grid2, bsplines_dim1, bsplines_dim2 >::Grid2 =  grid2;
```




<hr>



### typedef Idx1 

_1D index of a grid point along the first dimension._ 
```C++
using Patch< grid1, grid2, bsplines_dim1, bsplines_dim2 >::Idx1 =  Idx<Grid1>;
```




<hr>



### typedef Idx12 

_2D index of a grid point along the first and second dimensions._ 
```C++
using Patch< grid1, grid2, bsplines_dim1, bsplines_dim2 >::Idx12 =  Idx<Grid1, Grid2>;
```




<hr>



### typedef Idx2 

_1D index of a grid point along the second dimension._ 
```C++
using Patch< grid1, grid2, bsplines_dim1, bsplines_dim2 >::Idx2 =  Idx<Grid2>;
```




<hr>



### typedef IdxRange1 

_Index range of a grids over the first dimension._ 
```C++
using Patch< grid1, grid2, bsplines_dim1, bsplines_dim2 >::IdxRange1 =  IdxRange<Grid1>;
```




<hr>



### typedef IdxRange12 

_Index range of a grids over the first and second dimension._ 
```C++
using Patch< grid1, grid2, bsplines_dim1, bsplines_dim2 >::IdxRange12 =  IdxRange<Grid1, Grid2>;
```




<hr>



### typedef IdxRange2 

_Index range of a grids over the second dimension._ 
```C++
using Patch< grid1, grid2, bsplines_dim1, bsplines_dim2 >::IdxRange2 =  IdxRange<Grid2>;
```




<hr>



### typedef IdxRangeBS1 

_Index range of a grids over the first spline dimension._ 
```C++
using Patch< grid1, grid2, bsplines_dim1, bsplines_dim2 >::IdxRangeBS1 =  IdxRange<BSplines1>;
```




<hr>



### typedef IdxRangeBS12 

_Index range of a grids over the first and second spline dimension._ 
```C++
using Patch< grid1, grid2, bsplines_dim1, bsplines_dim2 >::IdxRangeBS12 =  IdxRange<BSplines1, BSplines2>;
```




<hr>



### typedef IdxRangeBS2 

_Index range of a grids over the second spline dimension._ 
```C++
using Patch< grid1, grid2, bsplines_dim1, bsplines_dim2 >::IdxRangeBS2 =  IdxRange<BSplines2>;
```




<hr>



### typedef IdxStep1 

_1D index step between grid points along the first dimension._ 
```C++
using Patch< grid1, grid2, bsplines_dim1, bsplines_dim2 >::IdxStep1 =  IdxStep<Grid1>;
```




<hr>



### typedef IdxStep12 

_2D index step between grid points along the first and second dimensions._ 
```C++
using Patch< grid1, grid2, bsplines_dim1, bsplines_dim2 >::IdxStep12 =  IdxStep<Grid1, Grid2>;
```




<hr>



### typedef IdxStep2 

_1D index step between grid points along the second dimension._ 
```C++
using Patch< grid1, grid2, bsplines_dim1, bsplines_dim2 >::IdxStep2 =  IdxStep<Grid2>;
```




<hr>
## Public Static Attributes Documentation




### variable n\_dims 

_The number of dimensions of the patch._ 
```C++
int constexpr Patch< grid1, grid2, bsplines_dim1, bsplines_dim2 >::n_dims;
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/multipatch/connectivity/patch.hpp`

