

# Class IPreallocatableInterpolator

**template &lt;class GridInterp, class... Grid1D&gt;**



[**ClassList**](annotated.md) **>** [**IPreallocatableInterpolator**](classIPreallocatableInterpolator.md)



_A class which provides access to an interpolating function which can be preallocated where useful._ [More...](#detailed-description)

* `#include <iinterpolator.hpp>`



Inherits the following classes: [IInterpolator](classIInterpolator.md)














## Public Types

| Type | Name |
| ---: | :--- |
| typedef typename [**IInterpolator**](classIInterpolator.md)&lt; GridInterp, Grid1D... &gt;::batched\_derivs\_idx\_range\_type | [**batched\_derivs\_idx\_range\_type**](#typedef-batched_derivs_idx_range_type)  <br>_The type of the whole index range on which derivatives are defined._  |
| typedef typename [**IInterpolator**](classIInterpolator.md)&lt; GridInterp, Grid1D... &gt;::deriv\_type | [**deriv\_type**](#typedef-deriv_type)  <br>_The type of the dimension representing derivatives._  |


## Public Types inherited from IInterpolator

See [IInterpolator](classIInterpolator.md)

| Type | Name |
| ---: | :--- |
| typedef ddc::replace\_dim\_of\_t&lt; IdxRange&lt; Grid1D... &gt;, GridInterp, [**deriv\_type**](classIInterpolator.md#typedef-deriv_type) &gt; | [**batched\_derivs\_idx\_range\_type**](classIInterpolator.md#typedef-batched_derivs_idx_range_type)  <br>_The type of the whole index range on which derivatives are defined._  |
| typedef ddc::Deriv&lt; typename GridInterp::continuous\_dimension\_type &gt; | [**deriv\_type**](classIInterpolator.md#typedef-deriv_type)  <br>_The type of the dimension representing derivatives._  |






































## Public Functions

| Type | Name |
| ---: | :--- |
| virtual [**batched\_derivs\_idx\_range\_type**](classIPreallocatableInterpolator.md#typedef-batched_derivs_idx_range_type) | [**batched\_derivs\_idx\_range\_xmax**](#function-batched_derivs_idx_range_xmax) (IdxRange&lt; Grid1D... &gt; idx\_range) override const<br>_Get the batched derivs index range on upper boundaries._  |
| virtual [**batched\_derivs\_idx\_range\_type**](classIPreallocatableInterpolator.md#typedef-batched_derivs_idx_range_type) | [**batched\_derivs\_idx\_range\_xmin**](#function-batched_derivs_idx_range_xmin) (IdxRange&lt; Grid1D... &gt; idx\_range) override const<br>_Get the batched derivs index range on lower boundaries._  |
|  Field&lt; double, IdxRange&lt; Grid1D... &gt; &gt; | [**operator()**](#function-operator) (Field&lt; double, IdxRange&lt; Grid1D... &gt; &gt; const inout\_data, Field&lt; const Coord&lt; typename GridInterp::continuous\_dimension\_type &gt;, IdxRange&lt; Grid1D... &gt; &gt; const coordinates, std::optional&lt; Field&lt; double const, [**batched\_derivs\_idx\_range\_type**](classIPreallocatableInterpolator.md#typedef-batched_derivs_idx_range_type) &gt; &gt; derivs\_xmin=std::nullopt, std::optional&lt; Field&lt; double const, [**batched\_derivs\_idx\_range\_type**](classIPreallocatableInterpolator.md#typedef-batched_derivs_idx_range_type) &gt; &gt; derivs\_xmax=std::nullopt) override const<br>_Approximate the value of a function at a set of coordinates using the current values at a known set of interpolation points by temporarily preallocating an_ [_**IInterpolator**_](classIInterpolator.md) _._ |
| virtual std::unique\_ptr&lt; [**IInterpolator**](classIInterpolator.md)&lt; GridInterp, Grid1D... &gt; &gt; | [**preallocate**](#function-preallocate) () const = 0<br>_Allocate an instance of an InterpolatorProxy to use as an_ [_**IInterpolator**_](classIInterpolator.md) _._ |
|   | [**~IPreallocatableInterpolator**](#function-ipreallocatableinterpolator) () override<br> |


## Public Functions inherited from IInterpolator

See [IInterpolator](classIInterpolator.md)

| Type | Name |
| ---: | :--- |
| virtual [**batched\_derivs\_idx\_range\_type**](classIInterpolator.md#typedef-batched_derivs_idx_range_type) | [**batched\_derivs\_idx\_range\_xmax**](classIInterpolator.md#function-batched_derivs_idx_range_xmax) (IdxRange&lt; Grid1D... &gt; idx\_range) const = 0<br>_Get the batched derivs index range on upper boundaries._  |
| virtual [**batched\_derivs\_idx\_range\_type**](classIInterpolator.md#typedef-batched_derivs_idx_range_type) | [**batched\_derivs\_idx\_range\_xmin**](classIInterpolator.md#function-batched_derivs_idx_range_xmin) (IdxRange&lt; Grid1D... &gt; idx\_range) const = 0<br>_Get the batched derivs index range on lower boundaries._  |
| virtual Field&lt; double, IdxRange&lt; Grid1D... &gt; &gt; | [**operator()**](classIInterpolator.md#function-operator) (Field&lt; double, IdxRange&lt; Grid1D... &gt; &gt; inout\_data, Field&lt; const Coord&lt; typename GridInterp::continuous\_dimension\_type &gt;, IdxRange&lt; Grid1D... &gt; &gt; coordinates, std::optional&lt; Field&lt; double const, [**batched\_derivs\_idx\_range\_type**](classIInterpolator.md#typedef-batched_derivs_idx_range_type) &gt; &gt; derivs\_xmin=std::nullopt, std::optional&lt; Field&lt; double const, [**batched\_derivs\_idx\_range\_type**](classIInterpolator.md#typedef-batched_derivs_idx_range_type) &gt; &gt; derivs\_xmax=std::nullopt) const = 0<br>_Approximate the value of a function at a set of coordinates using the current values at a known set of interpolation points._  |
| virtual  | [**~IInterpolator**](classIInterpolator.md#function-iinterpolator) () = default<br> |






















































## Detailed Description


An abstract class which implements a preallocate function returning an unique pointer to an [**IInterpolator**](classIInterpolator.md). A pointer is used so that the returned object can be any sub-class of [**IInterpolator**](classIInterpolator.md). The type (and thus the implementation of the operator) will be determined when the pointer is dereferenced.


The preallocate function should be used to allocate an instance of the [**IInterpolator**](classIInterpolator.md) before using it repeatedly. Once the preallocated object goes out of scope it will be deallocated. This means that objects of this class take up little or no space in memory.


An example of this is seen in [**BslAdvectionVelocity**](classBslAdvectionVelocity.md). The [**IPreallocatableInterpolator**](classIPreallocatableInterpolator.md) stored in the [**BslAdvectionVelocity**](classBslAdvectionVelocity.md) takes up no memory between advections, however during the execution of the BslAdvectionVelocity::operator() function the [**IPreallocatableInterpolator::preallocate()**](classIPreallocatableInterpolator.md#function-preallocate) function is called. This leads to the creation of an [**IInterpolator**](classIInterpolator.md) instance, ensuring that all buffers necessary for the interpolation during the advection are allocated before the [**IInterpolator**](classIInterpolator.md) is used for interpolation in the advection loop. This ensures that these buffers are only allocated once per advection at the start of the BslAdvectionVelocity::operator() function. At the end of this function the unique pointer goes out of scope and the buffers are deallocated. 


    
## Public Types Documentation




### typedef batched\_derivs\_idx\_range\_type 

_The type of the whole index range on which derivatives are defined._ 
```C++
using IPreallocatableInterpolator< GridInterp, Grid1D >::batched_derivs_idx_range_type =  typename IInterpolator<GridInterp, Grid1D...>::batched_derivs_idx_range_type;
```




<hr>



### typedef deriv\_type 

_The type of the dimension representing derivatives._ 
```C++
using IPreallocatableInterpolator< GridInterp, Grid1D >::deriv_type =  typename IInterpolator<GridInterp, Grid1D...>::deriv_type;
```




<hr>
## Public Functions Documentation




### function batched\_derivs\_idx\_range\_xmax 

_Get the batched derivs index range on upper boundaries._ 
```C++
inline virtual batched_derivs_idx_range_type IPreallocatableInterpolator::batched_derivs_idx_range_xmax (
    IdxRange< Grid1D... > idx_range
) override const
```



Dimension of interest IDimI is replaced with ddc::Deriv&lt;IDimI::continuous\_dimensions\_type&gt;. This is the index range on which derivatives on upper boundaries are defined.




**Parameters:**


* `idx_range` The index range of a single-species distribution function. 



**Returns:**

idx\_range The upper boundaries of this index range. 





        
Implements [*IInterpolator::batched\_derivs\_idx\_range\_xmax*](classIInterpolator.md#function-batched_derivs_idx_range_xmax)


<hr>



### function batched\_derivs\_idx\_range\_xmin 

_Get the batched derivs index range on lower boundaries._ 
```C++
inline virtual batched_derivs_idx_range_type IPreallocatableInterpolator::batched_derivs_idx_range_xmin (
    IdxRange< Grid1D... > idx_range
) override const
```



Dimension of interest IDimI is replaced with ddc::Deriv&lt;IDimI::continuous\_dimensions\_type&gt;. This is the index range on which derivatives on lower boundaries are defined.




**Parameters:**


* `idx_range` The index range of a single-species distribution function. 



**Returns:**

idx\_range The lower boundaries of this index range. 





        
Implements [*IInterpolator::batched\_derivs\_idx\_range\_xmin*](classIInterpolator.md#function-batched_derivs_idx_range_xmin)


<hr>



### function operator() 

_Approximate the value of a function at a set of coordinates using the current values at a known set of interpolation points by temporarily preallocating an_ [_**IInterpolator**_](classIInterpolator.md) _._
```C++
inline Field< double, IdxRange< Grid1D... > > IPreallocatableInterpolator::operator() (
    Field< double, IdxRange< Grid1D... > > const inout_data,
    Field< const Coord< typename GridInterp::continuous_dimension_type >, IdxRange< Grid1D... > > const coordinates,
    std::optional< Field< double const, batched_derivs_idx_range_type > > derivs_xmin=std::nullopt,
    std::optional< Field< double const, batched_derivs_idx_range_type > > derivs_xmax=std::nullopt
) override const
```





**Parameters:**


* `inout_data` On input: an array containing the value of the function at the interpolation points. On output: an array containing the value of the function at the coordinates. 
* `coordinates` The coordinates where the function should be evaluated. 
* `derivs_xmin` The values of the derivatives at the lower boundary (used only with splines and ddc::BoundCond::HERMITE lower boundary condition). 
* `derivs_xmax` The values of the derivatives at the upper boundary (used only with splines and ddc::BoundCond::HERMITE upper boundary condition).



**Returns:**

A reference to the inout\_data array containing the value of the function at the coordinates. 





        

<hr>



### function preallocate 

_Allocate an instance of an InterpolatorProxy to use as an_ [_**IInterpolator**_](classIInterpolator.md) _._
```C++
virtual std::unique_ptr< IInterpolator < GridInterp, Grid1D... > > IPreallocatableInterpolator::preallocate () const = 0
```



Allocate and return an unique pointer to an instance of an [**IInterpolator**](classIInterpolator.md).




**Returns:**

An allocated instance of an InterpolatorProxy.




**See also:** InterpolatorProxy 


**See also:** [**IInterpolator**](classIInterpolator.md) 



        

<hr>



### function ~IPreallocatableInterpolator 

```C++
IPreallocatableInterpolator::~IPreallocatableInterpolator () override
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/interpolation/iinterpolator.hpp`

