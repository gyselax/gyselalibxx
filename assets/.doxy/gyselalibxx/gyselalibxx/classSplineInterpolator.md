

# Class SplineInterpolator

**template &lt;class GridInterp, class BSplines, ddc::BoundCond BcMin, ddc::BoundCond BcMax, class LeftExtrapolationRule, class RightExtrapolationRule, ddc::SplineSolver Solver, class... Grid1D&gt;**



[**ClassList**](annotated.md) **>** [**SplineInterpolator**](classSplineInterpolator.md)



_A class for interpolating a function using splines._ [More...](#detailed-description)

* `#include <spline_interpolator.hpp>`



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
|   | [**SplineInterpolator**](#function-splineinterpolator) (BuilderType const & builder, EvaluatorType const & evaluator) <br>_Create a spline interpolator object._  |
| virtual batched\_derivs\_idx\_range\_type | [**batched\_derivs\_idx\_range\_xmax**](#function-batched_derivs_idx_range_xmax) (IdxRange&lt; Grid1D... &gt; idx\_range) override const<br>_Get the batched derivs index range on upper boundaries._  |
| virtual batched\_derivs\_idx\_range\_type | [**batched\_derivs\_idx\_range\_xmin**](#function-batched_derivs_idx_range_xmin) (IdxRange&lt; Grid1D... &gt; idx\_range) override const<br>_Get the batched derivs index range on lower boundaries._  |
|  Field&lt; double, IdxRange&lt; Grid1D... &gt; &gt; | [**operator()**](#function-operator) (Field&lt; double, IdxRange&lt; Grid1D... &gt; &gt; const inout\_data, ConstField&lt; Coord&lt; typename GridInterp::continuous\_dimension\_type &gt;, IdxRange&lt; Grid1D... &gt; &gt; const coordinates, std::optional&lt; batched\_deriv\_field\_type &gt; derivs\_xmin=std::nullopt, std::optional&lt; batched\_deriv\_field\_type &gt; derivs\_xmax=std::nullopt) override const<br>_Approximate the value of a function at a set of coordinates using the current values at a known set of interpolation points._  |
|   | [**~SplineInterpolator**](#function-splineinterpolator) () override<br> |


## Public Functions inherited from IInterpolator

See [IInterpolator](classIInterpolator.md)

| Type | Name |
| ---: | :--- |
| virtual [**batched\_derivs\_idx\_range\_type**](classIInterpolator.md#typedef-batched_derivs_idx_range_type) | [**batched\_derivs\_idx\_range\_xmax**](classIInterpolator.md#function-batched_derivs_idx_range_xmax) (IdxRange&lt; Grid1D... &gt; idx\_range) const = 0<br>_Get the batched derivs index range on upper boundaries._  |
| virtual [**batched\_derivs\_idx\_range\_type**](classIInterpolator.md#typedef-batched_derivs_idx_range_type) | [**batched\_derivs\_idx\_range\_xmin**](classIInterpolator.md#function-batched_derivs_idx_range_xmin) (IdxRange&lt; Grid1D... &gt; idx\_range) const = 0<br>_Get the batched derivs index range on lower boundaries._  |
| virtual Field&lt; double, IdxRange&lt; Grid1D... &gt; &gt; | [**operator()**](classIInterpolator.md#function-operator) (Field&lt; double, IdxRange&lt; Grid1D... &gt; &gt; inout\_data, Field&lt; const Coord&lt; typename GridInterp::continuous\_dimension\_type &gt;, IdxRange&lt; Grid1D... &gt; &gt; coordinates, std::optional&lt; Field&lt; double const, [**batched\_derivs\_idx\_range\_type**](classIInterpolator.md#typedef-batched_derivs_idx_range_type) &gt; &gt; derivs\_xmin=std::nullopt, std::optional&lt; Field&lt; double const, [**batched\_derivs\_idx\_range\_type**](classIInterpolator.md#typedef-batched_derivs_idx_range_type) &gt; &gt; derivs\_xmax=std::nullopt) const = 0<br>_Approximate the value of a function at a set of coordinates using the current values at a known set of interpolation points._  |
| virtual  | [**~IInterpolator**](classIInterpolator.md#function-iinterpolator) () = default<br> |






















































## Detailed Description


The class is parametrised by multiple template parameters. Please note that CTAD will deduce all these template parameters from the Builder and Evaluator passed as constructor arguments.




**Template parameters:**


* `GridInterp` The dimension over which we interpolate. 
* `BSplines` The BSplines along the dimension of interest. 
* `BcMin` The boundary condition at the lower boundary. 
* `BcMax` The boundary condition at the upper boundary. 
* `Grid1D...` All the dimensions of the interpolation problem (batched + interpolated). 




    
## Public Functions Documentation




### function SplineInterpolator 

_Create a spline interpolator object._ 
```C++
inline SplineInterpolator::SplineInterpolator (
    BuilderType const & builder,
    EvaluatorType const & evaluator
) 
```





**Parameters:**


* `builder` An operator which builds spline coefficients from the values of a function at known interpolation points. 
* `evaluator` An operator which evaluates the value of a spline at requested coordinates. 




        

<hr>



### function batched\_derivs\_idx\_range\_xmax 

_Get the batched derivs index range on upper boundaries._ 
```C++
inline virtual batched_derivs_idx_range_type SplineInterpolator::batched_derivs_idx_range_xmax (
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
inline virtual batched_derivs_idx_range_type SplineInterpolator::batched_derivs_idx_range_xmin (
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
inline Field< double, IdxRange< Grid1D... > > SplineInterpolator::operator() (
    Field< double, IdxRange< Grid1D... > > const inout_data,
    ConstField< Coord< typename GridInterp::continuous_dimension_type >, IdxRange< Grid1D... > > const coordinates,
    std::optional< batched_deriv_field_type > derivs_xmin=std::nullopt,
    std::optional< batched_deriv_field_type > derivs_xmax=std::nullopt
) override const
```





**Parameters:**


* `inout_data` On input: an array containing the value of the function at the interpolation points. On output: an array containing the value of the function at the coordinates. 
* `coordinates` The coordinates where the function should be evaluated. 
* `derivs_xmin` The values of the derivatives at the lower boundary (used only with ddc::BoundCond::HERMITE lower boundary condition). 
* `derivs_xmax` The values of the derivatives at the upper boundary (used only with ddc::BoundCond::HERMITE upper boundary condition).



**Returns:**

A reference to the inout\_data array containing the value of the function at the coordinates. 





        

<hr>



### function ~SplineInterpolator 

```C++
SplineInterpolator::~SplineInterpolator () override
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/interpolation/spline_interpolator.hpp`

