

# Class PreallocatableLagrangeInterpolator

**template &lt;class GridInterp, BCond BcMin, BCond BcMax, class... Grid1D&gt;**



[**ClassList**](annotated.md) **>** [**PreallocatableLagrangeInterpolator**](classPreallocatableLagrangeInterpolator.md)



_A class which stores information necessary to create an instance of the_ [_**LagrangeInterpolator**_](classLagrangeInterpolator.md) _class._[More...](#detailed-description)

* `#include <Lagrange_interpolator.hpp>`



Inherits the following classes: [IPreallocatableInterpolator](classIPreallocatableInterpolator.md)
















## Public Types inherited from IPreallocatableInterpolator

See [IPreallocatableInterpolator](classIPreallocatableInterpolator.md)

| Type | Name |
| ---: | :--- |
| typedef typename [**IInterpolator**](classIInterpolator.md)&lt; GridInterp, Grid1D... &gt;::batched\_derivs\_idx\_range\_type | [**batched\_derivs\_idx\_range\_type**](classIPreallocatableInterpolator.md#typedef-batched_derivs_idx_range_type)  <br>_The type of the whole index range on which derivatives are defined._  |
| typedef typename [**IInterpolator**](classIInterpolator.md)&lt; GridInterp, Grid1D... &gt;::deriv\_type | [**deriv\_type**](classIPreallocatableInterpolator.md#typedef-deriv_type)  <br>_The type of the dimension representing derivatives._  |


## Public Types inherited from IInterpolator

See [IInterpolator](classIInterpolator.md)

| Type | Name |
| ---: | :--- |
| typedef ddc::replace\_dim\_of\_t&lt; IdxRange&lt; Grid1D... &gt;, GridInterp, [**deriv\_type**](classIInterpolator.md#typedef-deriv_type) &gt; | [**batched\_derivs\_idx\_range\_type**](classIInterpolator.md#typedef-batched_derivs_idx_range_type)  <br>_The type of the whole index range on which derivatives are defined._  |
| typedef ddc::Deriv&lt; typename GridInterp::continuous\_dimension\_type &gt; | [**deriv\_type**](classIInterpolator.md#typedef-deriv_type)  <br>_The type of the dimension representing derivatives._  |
























































## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**PreallocatableLagrangeInterpolator**](#function-preallocatablelagrangeinterpolator) ([**LagrangeInterpolator**](classLagrangeInterpolator.md)&lt; GridInterp, BcMin, BcMax, Grid1D... &gt; const & evaluator) <br>_Create an object capable of creating_ [_**LagrangeInterpolator**_](classLagrangeInterpolator.md) _objects._ |
| virtual std::unique\_ptr&lt; [**IInterpolator**](classIInterpolator.md)&lt; GridInterp, Grid1D... &gt; &gt; | [**preallocate**](#function-preallocate) () override const<br> |
|   | [**~PreallocatableLagrangeInterpolator**](#function-preallocatablelagrangeinterpolator) () override<br> |


## Public Functions inherited from IPreallocatableInterpolator

See [IPreallocatableInterpolator](classIPreallocatableInterpolator.md)

| Type | Name |
| ---: | :--- |
| virtual [**batched\_derivs\_idx\_range\_type**](classIPreallocatableInterpolator.md#typedef-batched_derivs_idx_range_type) | [**batched\_derivs\_idx\_range\_xmax**](classIPreallocatableInterpolator.md#function-batched_derivs_idx_range_xmax) (IdxRange&lt; Grid1D... &gt; idx\_range) override const<br>_Get the batched derivs index range on upper boundaries._  |
| virtual [**batched\_derivs\_idx\_range\_type**](classIPreallocatableInterpolator.md#typedef-batched_derivs_idx_range_type) | [**batched\_derivs\_idx\_range\_xmin**](classIPreallocatableInterpolator.md#function-batched_derivs_idx_range_xmin) (IdxRange&lt; Grid1D... &gt; idx\_range) override const<br>_Get the batched derivs index range on lower boundaries._  |
|  Field&lt; double, IdxRange&lt; Grid1D... &gt; &gt; | [**operator()**](classIPreallocatableInterpolator.md#function-operator) (Field&lt; double, IdxRange&lt; Grid1D... &gt; &gt; const inout\_data, Field&lt; const Coord&lt; typename GridInterp::continuous\_dimension\_type &gt;, IdxRange&lt; Grid1D... &gt; &gt; const coordinates, std::optional&lt; Field&lt; double const, [**batched\_derivs\_idx\_range\_type**](classIPreallocatableInterpolator.md#typedef-batched_derivs_idx_range_type) &gt; &gt; derivs\_xmin=std::nullopt, std::optional&lt; Field&lt; double const, [**batched\_derivs\_idx\_range\_type**](classIPreallocatableInterpolator.md#typedef-batched_derivs_idx_range_type) &gt; &gt; derivs\_xmax=std::nullopt) override const<br>_Approximate the value of a function at a set of coordinates using the current values at a known set of interpolation points by temporarily preallocating an_ [_**IInterpolator**_](classIInterpolator.md) _._ |
| virtual std::unique\_ptr&lt; [**IInterpolator**](classIInterpolator.md)&lt; GridInterp, Grid1D... &gt; &gt; | [**preallocate**](classIPreallocatableInterpolator.md#function-preallocate) () const = 0<br>_Allocate an instance of an InterpolatorProxy to use as an_ [_**IInterpolator**_](classIInterpolator.md) _._ |
|   | [**~IPreallocatableInterpolator**](classIPreallocatableInterpolator.md#function-ipreallocatableinterpolator) () override<br> |


## Public Functions inherited from IInterpolator

See [IInterpolator](classIInterpolator.md)

| Type | Name |
| ---: | :--- |
| virtual [**batched\_derivs\_idx\_range\_type**](classIInterpolator.md#typedef-batched_derivs_idx_range_type) | [**batched\_derivs\_idx\_range\_xmax**](classIInterpolator.md#function-batched_derivs_idx_range_xmax) (IdxRange&lt; Grid1D... &gt; idx\_range) const = 0<br>_Get the batched derivs index range on upper boundaries._  |
| virtual [**batched\_derivs\_idx\_range\_type**](classIInterpolator.md#typedef-batched_derivs_idx_range_type) | [**batched\_derivs\_idx\_range\_xmin**](classIInterpolator.md#function-batched_derivs_idx_range_xmin) (IdxRange&lt; Grid1D... &gt; idx\_range) const = 0<br>_Get the batched derivs index range on lower boundaries._  |
| virtual Field&lt; double, IdxRange&lt; Grid1D... &gt; &gt; | [**operator()**](classIInterpolator.md#function-operator) (Field&lt; double, IdxRange&lt; Grid1D... &gt; &gt; inout\_data, Field&lt; const Coord&lt; typename GridInterp::continuous\_dimension\_type &gt;, IdxRange&lt; Grid1D... &gt; &gt; coordinates, std::optional&lt; Field&lt; double const, [**batched\_derivs\_idx\_range\_type**](classIInterpolator.md#typedef-batched_derivs_idx_range_type) &gt; &gt; derivs\_xmin=std::nullopt, std::optional&lt; Field&lt; double const, [**batched\_derivs\_idx\_range\_type**](classIInterpolator.md#typedef-batched_derivs_idx_range_type) &gt; &gt; derivs\_xmax=std::nullopt) const = 0<br>_Approximate the value of a function at a set of coordinates using the current values at a known set of interpolation points._  |
| virtual  | [**~IInterpolator**](classIInterpolator.md#function-iinterpolator) () = default<br> |
















































































## Detailed Description


This class allows an instance of the [**LagrangeInterpolator**](classLagrangeInterpolator.md) class where necessary. This allows the memory allocated in the private members of the Interpolator to be freed when the object is not in use. 


    
## Public Functions Documentation




### function PreallocatableLagrangeInterpolator 

_Create an object capable of creating_ [_**LagrangeInterpolator**_](classLagrangeInterpolator.md) _objects._
```C++
inline explicit PreallocatableLagrangeInterpolator::PreallocatableLagrangeInterpolator (
    LagrangeInterpolator < GridInterp, BcMin, BcMax, Grid1D... > const & evaluator
) 
```





**Parameters:**


* `evaluator` An operator which evaluates the value of the interpolation polynomial at requested coordinates. 




        

<hr>



### function preallocate 

```C++
inline virtual std::unique_ptr< IInterpolator < GridInterp, Grid1D... > > PreallocatableLagrangeInterpolator::preallocate () override const
```



Create an instance of the [**LagrangeInterpolator**](classLagrangeInterpolator.md) class.




**Returns:**

A unique pointer to an instance of the [**LagrangeInterpolator**](classLagrangeInterpolator.md) class. 





        
Implements [*IPreallocatableInterpolator::preallocate*](classIPreallocatableInterpolator.md#function-preallocate)


<hr>



### function ~PreallocatableLagrangeInterpolator 

```C++
PreallocatableLagrangeInterpolator::~PreallocatableLagrangeInterpolator () override
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/interpolation/Lagrange_interpolator.hpp`

