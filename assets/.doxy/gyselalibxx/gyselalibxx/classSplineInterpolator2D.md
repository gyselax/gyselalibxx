

# Class SplineInterpolator2D

**template &lt;class Spline2DBuilder, class Spline2DEvaluator&gt;**



[**ClassList**](annotated.md) **>** [**SplineInterpolator2D**](classSplineInterpolator2D.md)



_A class for interpolating a function using splines in polar coordinates._ [More...](#detailed-description)

* `#include <spline_interpolator_2d.hpp>`



Inherits the following classes: [IInterpolator2D](classIInterpolator2D.md)














## Public Types

| Type | Name |
| ---: | :--- |
| typedef ConstField&lt; [**CoordType**](classSplineInterpolator2D.md#typedef-coordtype), IdxRangeBatched &gt; | [**CConstFieldType**](#typedef-cconstfieldtype)  <br>_The type of a field containing the coordinates where a field should be evaluated._  |
| typedef Coord&lt; Dim1, Dim2 &gt; | [**CoordType**](#typedef-coordtype)  <br>_The type of a coordinate on the 2D plane._  |
| typedef DField&lt; IdxRangeBatched &gt; | [**DFieldType**](#typedef-dfieldtype)  <br>_The type of a field which can be interpolated._  |


## Public Types inherited from IInterpolator2D

See [IInterpolator2D](classIInterpolator2D.md)

| Type | Name |
| ---: | :--- |
| typedef ConstField&lt; [**CoordType**](classIInterpolator2D.md#typedef-coordtype), IdxRangeBatched &gt; | [**CConstFieldType**](classIInterpolator2D.md#typedef-cconstfieldtype)  <br>_The type of a field containing the coordinates where a field should be evaluated._  |
| typedef Coord&lt; Dim1, Dim2 &gt; | [**CoordType**](classIInterpolator2D.md#typedef-coordtype)  <br>_The type of a coordinate on the 2D plane._  |
| typedef DField&lt; IdxRangeBatched &gt; | [**DFieldType**](classIInterpolator2D.md#typedef-dfieldtype)  <br>_The type of a field which can be interpolated._  |






































## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**SplineInterpolator2D**](#function-splineinterpolator2d) (Spline2DBuilder const & builder, Spline2DEvaluator const & evaluator) <br>_Create a spline interpolator object._  |
|  [**DFieldType**](classSplineInterpolator2D.md#typedef-dfieldtype) | [**operator()**](#function-operator) ([**DFieldType**](classSplineInterpolator2D.md#typedef-dfieldtype) const inout\_data, [**CConstFieldType**](classSplineInterpolator2D.md#typedef-cconstfieldtype) const coordinates) override const<br>_Approximate the value of a function at a set of polar coordinates using the current values at a known set of interpolation points._  |


## Public Functions inherited from IInterpolator2D

See [IInterpolator2D](classIInterpolator2D.md)

| Type | Name |
| ---: | :--- |
| virtual DField&lt; IdxRangeBatched &gt; | [**operator()**](classIInterpolator2D.md#function-operator) (DField&lt; IdxRangeBatched &gt; inout\_data, ConstField&lt; [**CoordType**](classIInterpolator2D.md#typedef-coordtype), IdxRangeBatched &gt; coordinates) const = 0<br>_Approximate the value of a function at a set of coordinates using the current values at a known set of interpolation points._  |
| virtual  | [**~IInterpolator2D**](classIInterpolator2D.md#function-iinterpolator2d) () = default<br> |






















































## Detailed Description




**Template parameters:**


* `RadialExtrapolationRule` The extrapolation rule applied at the outer radial bound. 




    
## Public Types Documentation




### typedef CConstFieldType 

_The type of a field containing the coordinates where a field should be evaluated._ 
```C++
using IInterpolator2D< IdxRange2D, IdxRangeBatched >::CConstFieldType =  ConstField<CoordType, IdxRangeBatched>;
```




<hr>



### typedef CoordType 

_The type of a coordinate on the 2D plane._ 
```C++
using IInterpolator2D< IdxRange2D, IdxRangeBatched >::CoordType =  Coord<Dim1, Dim2>;
```




<hr>



### typedef DFieldType 

_The type of a field which can be interpolated._ 
```C++
using IInterpolator2D< IdxRange2D, IdxRangeBatched >::DFieldType =  DField<IdxRangeBatched>;
```




<hr>
## Public Functions Documentation




### function SplineInterpolator2D 

_Create a spline interpolator object._ 
```C++
inline SplineInterpolator2D::SplineInterpolator2D (
    Spline2DBuilder const & builder,
    Spline2DEvaluator const & evaluator
) 
```





**Parameters:**


* `builder` An operator which builds spline coefficients from the values of a function at known interpolation points. 
* `evaluator` An operator which evaluates the value of a spline at requested coordinates. 




        

<hr>



### function operator() 

_Approximate the value of a function at a set of polar coordinates using the current values at a known set of interpolation points._ 
```C++
inline DFieldType SplineInterpolator2D::operator() (
    DFieldType const inout_data,
    CConstFieldType const coordinates
) override const
```





**Parameters:**


* `inout_data` On input: an array containing the value of the function at the interpolation points. On output: an array containing the value of the function at the coordinates. 
* `coordinates` The polar coordinates where the function should be evaluated.



**Returns:**

A reference to the inout\_data array containing the value of the function at the coordinates. 





        

<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/interpolation/spline_interpolator_2d.hpp`

