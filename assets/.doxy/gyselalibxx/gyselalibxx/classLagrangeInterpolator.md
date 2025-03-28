

# Class LagrangeInterpolator

**template &lt;class GridInterp, BCond BcMin, BCond BcMax, class... Grid1D&gt;**



[**ClassList**](annotated.md) **>** [**LagrangeInterpolator**](classLagrangeInterpolator.md)



_A class for interpolating a function using_ [_**Lagrange**_](classLagrange.md) _polynomials. It is designed to work with both uniform and non-uniform mesh, and have the advantage to be local._

* `#include <Lagrange_interpolator.hpp>`



Inherits the following classes: [IInterpolator](classIInterpolator.md)
















## Public Types inherited from IInterpolator

See [IInterpolator](classIInterpolator.md)

| Type | Name |
| ---: | :--- |
| typedef ddc::replace\_dim\_of\_t&lt; IdxRange&lt; Grid1D... &gt;, GridInterp, [**deriv\_type**](classIInterpolator.md#typedef-deriv_type) &gt; | [**batched\_derivs\_idx\_range\_type**](classIInterpolator.md#typedef-batched_derivs_idx_range_type)  <br>_The type of the whole index range on which derivatives are defined._  |
| typedef ddc::Deriv&lt; typename GridInterp::continuous\_dimension\_type &gt; | [**deriv\_type**](classIInterpolator.md#typedef-deriv_type)  <br>_The type of the dimension representing derivatives._  |






































## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**LagrangeInterpolator**](#function-lagrangeinterpolator) (int degree, IdxStep&lt; GridInterp &gt; ghost) <br>_Create a_ [_**Lagrange**_](classLagrange.md) _interpolator object._ |
| virtual batched\_derivs\_idx\_range\_type | [**batched\_derivs\_idx\_range\_xmax**](#function-batched_derivs_idx_range_xmax) (IdxRange&lt; Grid1D... &gt; idx\_range) override const<br>_Get the batched derivs index range on upper boundaries._  |
| virtual batched\_derivs\_idx\_range\_type | [**batched\_derivs\_idx\_range\_xmin**](#function-batched_derivs_idx_range_xmin) (IdxRange&lt; Grid1D... &gt; idx\_range) override const<br>_Get the batched derivs index range on lower boundaries._  |
|  Field&lt; double, IdxRange&lt; Grid1D... &gt; &gt; | [**operator()**](#function-operator) (Field&lt; double, IdxRange&lt; Grid1D... &gt; &gt; inout\_data, Field&lt; const Coord&lt; typename GridInterp::continuous\_dimension\_type &gt;, IdxRange&lt; Grid1D... &gt; &gt; coordinates, std::optional&lt; Field&lt; double const, typename [**IInterpolator**](classIInterpolator.md)&lt; GridInterp, Grid1D... &gt;::batched\_derivs\_idx\_range\_type &gt; &gt; derivs\_xmin=std::nullopt, std::optional&lt; Field&lt; double const, typename [**IInterpolator**](classIInterpolator.md)&lt; GridInterp, Grid1D... &gt;::batched\_derivs\_idx\_range\_type &gt; &gt; derivs\_xmax=std::nullopt) override const<br>_Approximate the value of a function at a set of coordinates using the current values at a known set of interpolation points._  |
|   | [**~LagrangeInterpolator**](#function-lagrangeinterpolator) () override<br> |


## Public Functions inherited from IInterpolator

See [IInterpolator](classIInterpolator.md)

| Type | Name |
| ---: | :--- |
| virtual [**batched\_derivs\_idx\_range\_type**](classIInterpolator.md#typedef-batched_derivs_idx_range_type) | [**batched\_derivs\_idx\_range\_xmax**](classIInterpolator.md#function-batched_derivs_idx_range_xmax) (IdxRange&lt; Grid1D... &gt; idx\_range) const = 0<br>_Get the batched derivs index range on upper boundaries._  |
| virtual [**batched\_derivs\_idx\_range\_type**](classIInterpolator.md#typedef-batched_derivs_idx_range_type) | [**batched\_derivs\_idx\_range\_xmin**](classIInterpolator.md#function-batched_derivs_idx_range_xmin) (IdxRange&lt; Grid1D... &gt; idx\_range) const = 0<br>_Get the batched derivs index range on lower boundaries._  |
| virtual Field&lt; double, IdxRange&lt; Grid1D... &gt; &gt; | [**operator()**](classIInterpolator.md#function-operator) (Field&lt; double, IdxRange&lt; Grid1D... &gt; &gt; inout\_data, Field&lt; const Coord&lt; typename GridInterp::continuous\_dimension\_type &gt;, IdxRange&lt; Grid1D... &gt; &gt; coordinates, std::optional&lt; Field&lt; double const, [**batched\_derivs\_idx\_range\_type**](classIInterpolator.md#typedef-batched_derivs_idx_range_type) &gt; &gt; derivs\_xmin=std::nullopt, std::optional&lt; Field&lt; double const, [**batched\_derivs\_idx\_range\_type**](classIInterpolator.md#typedef-batched_derivs_idx_range_type) &gt; &gt; derivs\_xmax=std::nullopt) const = 0<br>_Approximate the value of a function at a set of coordinates using the current values at a known set of interpolation points._  |
| virtual  | [**~IInterpolator**](classIInterpolator.md#function-iinterpolator) () = default<br> |






















































## Public Functions Documentation




### function LagrangeInterpolator 

_Create a_ [_**Lagrange**_](classLagrange.md) _interpolator object._
```C++
inline LagrangeInterpolator::LagrangeInterpolator (
    int degree,
    IdxStep< GridInterp > ghost
) 
```





**Parameters:**


* `degree` Degree of polynomials 
* `ghost` Discrete vector which gives the number of ghost points. By default choose 2. 




        

<hr>



### function batched\_derivs\_idx\_range\_xmax 

_Get the batched derivs index range on upper boundaries._ 
```C++
inline virtual batched_derivs_idx_range_type LagrangeInterpolator::batched_derivs_idx_range_xmax (
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
inline virtual batched_derivs_idx_range_type LagrangeInterpolator::batched_derivs_idx_range_xmin (
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

_Approximate the value of a function at a set of coordinates using the current values at a known set of interpolation points._ 
```C++
inline Field< double, IdxRange< Grid1D... > > LagrangeInterpolator::operator() (
    Field< double, IdxRange< Grid1D... > > inout_data,
    Field< const Coord< typename GridInterp::continuous_dimension_type >, IdxRange< Grid1D... > > coordinates,
    std::optional< Field< double const, typename IInterpolator < GridInterp, Grid1D... >::batched_derivs_idx_range_type > > derivs_xmin=std::nullopt,
    std::optional< Field< double const, typename IInterpolator < GridInterp, Grid1D... >::batched_derivs_idx_range_type > > derivs_xmax=std::nullopt
) override const
```





**Parameters:**


* `inout_data` On input: an array containing the value of the function at the interpolation points. On output: an array containing the value of the function at the coordinates. 
* `coordinates` The coordinates where the function should be evaluated. 
* `derivs_xmin` The values of the derivatives at the lower boundary (unused in this class). 
* `derivs_xmax` The values of the derivatives at the upper boundary (unused in this class).



**Returns:**

A reference to the inout\_data array containing the value of the function at the coordinates. 





        

<hr>



### function ~LagrangeInterpolator 

```C++
LagrangeInterpolator::~LagrangeInterpolator () override
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/interpolation/Lagrange_interpolator.hpp`

