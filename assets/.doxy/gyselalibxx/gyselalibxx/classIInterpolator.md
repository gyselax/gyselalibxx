

# Class IInterpolator

**template &lt;class GridInterp, class... Grid1D&gt;**



[**ClassList**](annotated.md) **>** [**IInterpolator**](classIInterpolator.md)



_A class which provides an interpolating function._ [More...](#detailed-description)

* `#include <iinterpolator.hpp>`

















## Public Types

| Type | Name |
| ---: | :--- |
| typedef ddc::replace\_dim\_of\_t&lt; IdxRange&lt; Grid1D... &gt;, GridInterp, [**deriv\_type**](classIInterpolator.md#typedef-deriv_type) &gt; | [**batched\_derivs\_idx\_range\_type**](#typedef-batched_derivs_idx_range_type)  <br>_The type of the whole index range on which derivatives are defined._  |
| typedef ddc::Deriv&lt; typename GridInterp::continuous\_dimension\_type &gt; | [**deriv\_type**](#typedef-deriv_type)  <br>_The type of the dimension representing derivatives._  |




















## Public Functions

| Type | Name |
| ---: | :--- |
| virtual [**batched\_derivs\_idx\_range\_type**](classIInterpolator.md#typedef-batched_derivs_idx_range_type) | [**batched\_derivs\_idx\_range\_xmax**](#function-batched_derivs_idx_range_xmax) (IdxRange&lt; Grid1D... &gt; idx\_range) const = 0<br>_Get the batched derivs index range on upper boundaries._  |
| virtual [**batched\_derivs\_idx\_range\_type**](classIInterpolator.md#typedef-batched_derivs_idx_range_type) | [**batched\_derivs\_idx\_range\_xmin**](#function-batched_derivs_idx_range_xmin) (IdxRange&lt; Grid1D... &gt; idx\_range) const = 0<br>_Get the batched derivs index range on lower boundaries._  |
| virtual Field&lt; double, IdxRange&lt; Grid1D... &gt; &gt; | [**operator()**](#function-operator) (Field&lt; double, IdxRange&lt; Grid1D... &gt; &gt; inout\_data, Field&lt; const Coord&lt; typename GridInterp::continuous\_dimension\_type &gt;, IdxRange&lt; Grid1D... &gt; &gt; coordinates, std::optional&lt; Field&lt; double const, [**batched\_derivs\_idx\_range\_type**](classIInterpolator.md#typedef-batched_derivs_idx_range_type) &gt; &gt; derivs\_xmin=std::nullopt, std::optional&lt; Field&lt; double const, [**batched\_derivs\_idx\_range\_type**](classIInterpolator.md#typedef-batched_derivs_idx_range_type) &gt; &gt; derivs\_xmax=std::nullopt) const = 0<br>_Approximate the value of a function at a set of coordinates using the current values at a known set of interpolation points._  |
| virtual  | [**~IInterpolator**](#function-iinterpolator) () = default<br> |




























## Detailed Description


An abstract class which implements a function allowing the value of a function to be approximated at a set of coordinates from a set of known values of the function. 


    
## Public Types Documentation




### typedef batched\_derivs\_idx\_range\_type 

_The type of the whole index range on which derivatives are defined._ 
```C++
using IInterpolator< GridInterp, Grid1D >::batched_derivs_idx_range_type =  ddc::replace_dim_of_t<IdxRange<Grid1D...>, GridInterp, deriv_type>;
```




<hr>



### typedef deriv\_type 

_The type of the dimension representing derivatives._ 
```C++
using IInterpolator< GridInterp, Grid1D >::deriv_type =  ddc::Deriv<typename GridInterp::continuous_dimension_type>;
```




<hr>
## Public Functions Documentation




### function batched\_derivs\_idx\_range\_xmax 

_Get the batched derivs index range on upper boundaries._ 
```C++
virtual batched_derivs_idx_range_type IInterpolator::batched_derivs_idx_range_xmax (
    IdxRange< Grid1D... > idx_range
) const = 0
```



Dimension of interest IDimI is replaced with ddc::Deriv&lt;IDimI::continuous\_dimensions\_type&gt;. This is the index range on which derivatives on upper boundaries are defined.




**Parameters:**


* `idx_range` The index range of a single-species distribution function. 



**Returns:**

idx\_range The upper boundaries of this index range. 





        

<hr>



### function batched\_derivs\_idx\_range\_xmin 

_Get the batched derivs index range on lower boundaries._ 
```C++
virtual batched_derivs_idx_range_type IInterpolator::batched_derivs_idx_range_xmin (
    IdxRange< Grid1D... > idx_range
) const = 0
```



Dimension of interest IDimI is replaced with ddc::Deriv&lt;IDimI::continuous\_dimensions\_type&gt;. This is the index range on which derivatives on lower boundaries are defined.




**Parameters:**


* `idx_range` The index range of a single-species distribution function. 



**Returns:**

idx\_range The lower boundaries of this index range. 





        

<hr>



### function operator() 

_Approximate the value of a function at a set of coordinates using the current values at a known set of interpolation points._ 
```C++
virtual Field< double, IdxRange< Grid1D... > > IInterpolator::operator() (
    Field< double, IdxRange< Grid1D... > > inout_data,
    Field< const Coord< typename GridInterp::continuous_dimension_type >, IdxRange< Grid1D... > > coordinates,
    std::optional< Field< double const, batched_derivs_idx_range_type > > derivs_xmin=std::nullopt,
    std::optional< Field< double const, batched_derivs_idx_range_type > > derivs_xmax=std::nullopt
) const = 0
```





**Parameters:**


* `inout_data` On input: an array containing the value of the function at the interpolation points. On output: an array containing the value of the function at the coordinates. 
* `coordinates` The coordinates where the function should be evaluated. 
* `derivs_xmin` The values of the derivatives at the lower boundary (used only with splines and ddc::BoundCond::HERMITE lower boundary condition). 
* `derivs_xmax` The values of the derivatives at the upper boundary (used only with splines and ddc::BoundCond::HERMITE upper boundary condition).



**Returns:**

A reference to the inout\_data array containing the value of the function at the coordinates. 





        

<hr>



### function ~IInterpolator 

```C++
virtual IInterpolator::~IInterpolator () = default
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/interpolation/iinterpolator.hpp`

