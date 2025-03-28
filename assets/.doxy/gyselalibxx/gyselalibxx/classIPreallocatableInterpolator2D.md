

# Class IPreallocatableInterpolator2D

**template &lt;class IdxRange2D, class IdxRangeBatched&gt;**



[**ClassList**](annotated.md) **>** [**IPreallocatableInterpolator2D**](classIPreallocatableInterpolator2D.md)



_A class which provides access to an interpolating function which can be preallocated where useful._ [More...](#detailed-description)

* `#include <i_interpolator_2d.hpp>`



Inherits the following classes: [IInterpolator2D](classIInterpolator2D.md)
















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
| virtual DField&lt; IdxRangeBatched &gt; | [**operator()**](#function-operator) (DField&lt; IdxRangeBatched &gt; inout\_data, ConstField&lt; [**CoordType**](classIInterpolator2D.md#typedef-coordtype), IdxRangeBatched &gt; coordinates) override const<br>_Approximate the value of a function at a set of coordinates using the current values at a known set of interpolation points._  |
| virtual std::unique\_ptr&lt; [**IInterpolator2D**](classIInterpolator2D.md)&lt; IdxRange2D, IdxRangeBatched &gt; &gt; | [**preallocate**](#function-preallocate) () const = 0<br>_Allocate an instance of a pointer to an InterpolatorRTheta._  |


## Public Functions inherited from IInterpolator2D

See [IInterpolator2D](classIInterpolator2D.md)

| Type | Name |
| ---: | :--- |
| virtual DField&lt; IdxRangeBatched &gt; | [**operator()**](classIInterpolator2D.md#function-operator) (DField&lt; IdxRangeBatched &gt; inout\_data, ConstField&lt; [**CoordType**](classIInterpolator2D.md#typedef-coordtype), IdxRangeBatched &gt; coordinates) const = 0<br>_Approximate the value of a function at a set of coordinates using the current values at a known set of interpolation points._  |
| virtual  | [**~IInterpolator2D**](classIInterpolator2D.md#function-iinterpolator2d) () = default<br> |






















































## Detailed Description


An abstract class which implements a preallocate function returning a pointer to an InterpolatorRTheta. A pointer to an InterpolatorRTheta is used so that the returned object can be any sub-class of [**IInterpolator2D**](classIInterpolator2D.md). The type (and thus the implementation of the operator) will be determined when the pointer is dereferenced.


The preallocate function should be used to allocate an instance of the [**IInterpolator2D**](classIInterpolator2D.md) before using it repeatedly. Once the preallocated object goes out of scope it will be deallocated. This means that objects of this class take up little or no space in memory. 


    
## Public Functions Documentation




### function operator() 

_Approximate the value of a function at a set of coordinates using the current values at a known set of interpolation points._ 
```C++
inline virtual DField< IdxRangeBatched > IPreallocatableInterpolator2D::operator() (
    DField< IdxRangeBatched > inout_data,
    ConstField< CoordType , IdxRangeBatched > coordinates
) override const
```





**Parameters:**


* `inout_data` On input: an array containing the value of the function at the interpolation points. On output: an array containing the value of the function at the coordinates. 
* `coordinates` The coordinates where the function should be evaluated.



**Returns:**

A reference to the inout\_data array containing the value of the function at the coordinates. 





        
Implements [*IInterpolator2D::operator()*](classIInterpolator2D.md#function-operator)


<hr>



### function preallocate 

_Allocate an instance of a pointer to an InterpolatorRTheta._ 
```C++
virtual std::unique_ptr< IInterpolator2D < IdxRange2D, IdxRangeBatched > > IPreallocatableInterpolator2D::preallocate () const = 0
```



Allocate and return an instance of a pointer to a sub-class of [**IInterpolator2D**](classIInterpolator2D.md).




**Returns:**

A pointer to an InterpolatorRTheta.




**See also:** [**IInterpolator2D**](classIInterpolator2D.md) 



        

<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/interpolation/i_interpolator_2d.hpp`

