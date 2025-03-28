

# Class IInterpolator2D

**template &lt;class IdxRange2D, class IdxRangeBatched&gt;**



[**ClassList**](annotated.md) **>** [**IInterpolator2D**](classIInterpolator2D.md)



_A class which provides an interpolating function._ [More...](#detailed-description)

* `#include <i_interpolator_2d.hpp>`





Inherited by the following classes: [IPreallocatableInterpolator2D](classIPreallocatableInterpolator2D.md),  [IPreallocatableInterpolator2D](classIPreallocatableInterpolator2D.md)












## Public Types

| Type | Name |
| ---: | :--- |
| typedef ConstField&lt; [**CoordType**](classIInterpolator2D.md#typedef-coordtype), IdxRangeBatched &gt; | [**CConstFieldType**](#typedef-cconstfieldtype)  <br>_The type of a field containing the coordinates where a field should be evaluated._  |
| typedef Coord&lt; Dim1, Dim2 &gt; | [**CoordType**](#typedef-coordtype)  <br>_The type of a coordinate on the 2D plane._  |
| typedef DField&lt; IdxRangeBatched &gt; | [**DFieldType**](#typedef-dfieldtype)  <br>_The type of a field which can be interpolated._  |




















## Public Functions

| Type | Name |
| ---: | :--- |
| virtual DField&lt; IdxRangeBatched &gt; | [**operator()**](#function-operator) (DField&lt; IdxRangeBatched &gt; inout\_data, ConstField&lt; [**CoordType**](classIInterpolator2D.md#typedef-coordtype), IdxRangeBatched &gt; coordinates) const = 0<br>_Approximate the value of a function at a set of coordinates using the current values at a known set of interpolation points._  |
| virtual  | [**~IInterpolator2D**](#function-iinterpolator2d) () = default<br> |




























## Detailed Description


An abstract class which implements a function allowing the value of a function to be approximated at a set of coordinates from a set of known values of the function. 


    
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




### function operator() 

_Approximate the value of a function at a set of coordinates using the current values at a known set of interpolation points._ 
```C++
virtual DField< IdxRangeBatched > IInterpolator2D::operator() (
    DField< IdxRangeBatched > inout_data,
    ConstField< CoordType , IdxRangeBatched > coordinates
) const = 0
```





**Parameters:**


* `inout_data` On input: an array containing the value of the function at the interpolation points. On output: an array containing the value of the function at the coordinates. 
* `coordinates` The coordinates where the function should be evaluated.



**Returns:**

A reference to the inout\_data array containing the value of the function at the coordinates. 





        

<hr>



### function ~IInterpolator2D 

```C++
virtual IInterpolator2D::~IInterpolator2D () = default
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/interpolation/i_interpolator_2d.hpp`

