

# Class PreallocatableSplineInterpolator2D

**template &lt;class Spline2DBuilder, class Spline2DEvaluator, class IdxRangeBatched&gt;**



[**ClassList**](annotated.md) **>** [**PreallocatableSplineInterpolator2D**](classPreallocatableSplineInterpolator2D.md)



_A class which stores information necessary to create a pointer to an instance of the_ [_**SplineInterpolator2D**_](classSplineInterpolator2D.md) _class._[More...](#detailed-description)

* `#include <spline_interpolator_2d.hpp>`



Inherits the following classes: [IPreallocatableInterpolator2D](classIPreallocatableInterpolator2D.md)


















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
|   | [**PreallocatableSplineInterpolator2D**](#function-preallocatablesplineinterpolator2d) (Spline2DBuilder const & builder, Spline2DEvaluator const & evaluator, IdxRangeBatched idx\_range\_batched) <br>_Create an object capable of creating_ [_**SplineInterpolator2D**_](classSplineInterpolator2D.md) _objects._ |
| virtual std::unique\_ptr&lt; [**IInterpolator2D**](classIInterpolator2D.md)&lt; typename Spline2DBuilder::interpolation\_domain\_type, IdxRangeBatched &gt; &gt; | [**preallocate**](#function-preallocate) () override const<br> |


## Public Functions inherited from IPreallocatableInterpolator2D

See [IPreallocatableInterpolator2D](classIPreallocatableInterpolator2D.md)

| Type | Name |
| ---: | :--- |
| virtual DField&lt; IdxRangeBatched &gt; | [**operator()**](classIPreallocatableInterpolator2D.md#function-operator) (DField&lt; IdxRangeBatched &gt; inout\_data, ConstField&lt; [**CoordType**](classIInterpolator2D.md#typedef-coordtype), IdxRangeBatched &gt; coordinates) override const<br>_Approximate the value of a function at a set of coordinates using the current values at a known set of interpolation points._  |
| virtual std::unique\_ptr&lt; [**IInterpolator2D**](classIInterpolator2D.md)&lt; IdxRange2D, IdxRangeBatched &gt; &gt; | [**preallocate**](classIPreallocatableInterpolator2D.md#function-preallocate) () const = 0<br>_Allocate an instance of a pointer to an InterpolatorRTheta._  |


## Public Functions inherited from IInterpolator2D

See [IInterpolator2D](classIInterpolator2D.md)

| Type | Name |
| ---: | :--- |
| virtual DField&lt; IdxRangeBatched &gt; | [**operator()**](classIInterpolator2D.md#function-operator) (DField&lt; IdxRangeBatched &gt; inout\_data, ConstField&lt; [**CoordType**](classIInterpolator2D.md#typedef-coordtype), IdxRangeBatched &gt; coordinates) const = 0<br>_Approximate the value of a function at a set of coordinates using the current values at a known set of interpolation points._  |
| virtual  | [**~IInterpolator2D**](classIInterpolator2D.md#function-iinterpolator2d) () = default<br> |
















































































## Detailed Description


This class allows an instance of the [**SplineInterpolator2D**](classSplineInterpolator2D.md) class to be instantiated where necessary. This allows the memory allocated in the private members of the [**SplineInterpolator2D**](classSplineInterpolator2D.md) to be freed when the object is not in use. These objects are: m\_coefs.


The class is parametrised by multiple template parameters. Please note that CTAD will deduce all these template parameters from the Builder and Evaluator passed as constructor arguments.




**Template parameters:**


* `Spline2DBuilder` The type of the 2D spline builder. 
* `Spline2DEvaluator` The type of the 2D spline evaluator. 
* `IdxRangeBatched` The type af the index range over which this operator will operate. This is necessary to define the internal spline representation. 




    
## Public Functions Documentation




### function PreallocatableSplineInterpolator2D 

_Create an object capable of creating_ [_**SplineInterpolator2D**_](classSplineInterpolator2D.md) _objects._
```C++
inline PreallocatableSplineInterpolator2D::PreallocatableSplineInterpolator2D (
    Spline2DBuilder const & builder,
    Spline2DEvaluator const & evaluator,
    IdxRangeBatched idx_range_batched
) 
```





**Parameters:**


* `builder` An operator which builds spline coefficients from the values of a function at known interpolation points. 
* `evaluator` An operator which evaluates the value of a spline at requested coordinates. 
* `idx_range_batched` The index range on which this operator operates. 




        

<hr>



### function preallocate 

```C++
inline virtual std::unique_ptr< IInterpolator2D < typename Spline2DBuilder::interpolation_domain_type, IdxRangeBatched > > PreallocatableSplineInterpolator2D::preallocate () override const
```



Create a pointer to an instance of the [**SplineInterpolator2D**](classSplineInterpolator2D.md) class.




**Returns:**

A pointer to an instance of the [**SplineInterpolator2D**](classSplineInterpolator2D.md) class. 





        
Implements [*IPreallocatableInterpolator2D::preallocate*](classIPreallocatableInterpolator2D.md#function-preallocate)


<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/interpolation/spline_interpolator_2d.hpp`

