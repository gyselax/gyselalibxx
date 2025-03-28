

# Class DiscreteToCartesianBuilder

**template &lt;class [**X**](structX.md), class [**Y**](structY.md), class SplineBuilder, class SplineEvaluator&gt;**



[**ClassList**](annotated.md) **>** [**DiscreteToCartesianBuilder**](classDiscreteToCartesianBuilder.md)



_A class to create a_ [_**DiscreteToCartesian**_](classDiscreteToCartesian.md) _instance from an analytical mapping. This class creates and stores splines memory spaces describing the analytical mapping. The discrete mapping is then created using the splines without copying data._[More...](#detailed-description)

* `#include <discrete_mapping_builder.hpp>`

















## Public Types

| Type | Name |
| ---: | :--- |
| typedef [**DiscreteToCartesian**](classDiscreteToCartesian.md)&lt; [**X**](structX.md), [**Y**](structY.md), SplineEvaluator &gt; | [**MappingType**](#typedef-mappingtype)  <br>_The type of the mapping that will be created._  |




















## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**DiscreteToCartesianBuilder**](#function-discretetocartesianbuilder) (ExecSpace exec\_space, Mapping const & analytical\_mapping, SplineBuilder const & builder, SplineEvaluator const & evaluator) <br>_Create an instance of the class capable of providing a_ [_**DiscreteToCartesian**_](classDiscreteToCartesian.md) _class instance._ |
|  [**DiscreteToCartesian**](classDiscreteToCartesian.md)&lt; [**X**](structX.md), [**Y**](structY.md), SplineEvaluator &gt; | [**operator()**](#function-operator) () const<br>_Get a_ [_**DiscreteToCartesian**_](classDiscreteToCartesian.md) _class instance._ |
|  void | [**set\_curvilinear\_to\_cartesian\_values**](#function-set_curvilinear_to_cartesian_values) (InterpolationField const & curvilinear\_to\_x\_vals, InterpolationField const & curvilinear\_to\_y\_vals, Mapping const & analytical\_mapping, IdxRangeInterpolationPoints const & interpolation\_idx\_range) <br>_Fill in the curvilinear fields with interpolation points mapped with the given analytical mapping._  |




























## Detailed Description




**Template parameters:**


* [**X**](structX.md) The first Cartesian dimension. 
* [**Y**](structY.md) The second Cartesian dimension. 
* `SplineBuilder` An operator for building spline coefficients. 
* `SplineEvaluator` An operator for evaluating a spline. 




    
## Public Types Documentation




### typedef MappingType 

_The type of the mapping that will be created._ 
```C++
using DiscreteToCartesianBuilder< X, Y, SplineBuilder, SplineEvaluator >::MappingType =  DiscreteToCartesian<X, Y, SplineEvaluator>;
```




<hr>
## Public Functions Documentation




### function DiscreteToCartesianBuilder 

_Create an instance of the class capable of providing a_ [_**DiscreteToCartesian**_](classDiscreteToCartesian.md) _class instance._
```C++
template<class Mapping>
inline DiscreteToCartesianBuilder::DiscreteToCartesianBuilder (
    ExecSpace exec_space,
    Mapping const & analytical_mapping,
    SplineBuilder const & builder,
    SplineEvaluator const & evaluator
) 
```





**Parameters:**


* `exec_space` The execution space where this class runs any for loops. 
* `analytical_mapping` The analytical mapping to be described by this discrete mapping. 
* `builder` A spline builder to be used to create a spline approximating the analytical mapping. 
* `evaluator` A spline evaluator to be used to evaluate a spline approximating the analytical mapping. 




        

<hr>



### function operator() 

_Get a_ [_**DiscreteToCartesian**_](classDiscreteToCartesian.md) _class instance._
```C++
inline DiscreteToCartesian < X , Y , SplineEvaluator > DiscreteToCartesianBuilder::operator() () const
```





**Returns:**

An instance of the mapping. 





        

<hr>



### function set\_curvilinear\_to\_cartesian\_values 

_Fill in the curvilinear fields with interpolation points mapped with the given analytical mapping._ 
```C++
template<class Mapping>
inline void DiscreteToCartesianBuilder::set_curvilinear_to_cartesian_values (
    InterpolationField const & curvilinear_to_x_vals,
    InterpolationField const & curvilinear_to_y_vals,
    Mapping const & analytical_mapping,
    IdxRangeInterpolationPoints const & interpolation_idx_range
) 
```



This function should be private. It is not due to the inclusion of a KOKKOS\_LAMBDA




**Template parameters:**


* `Mapping` Type of the analytical mapping. 



**Parameters:**


* `curvilinear_to_x_vals` Field of coordinate on [**X**](structX.md). 
* `curvilinear_to_y_vals` Field of coordinate on [**Y**](structY.md). 
* `analytical_mapping` Analytical mapping. 
* `interpolation_idx_range` Index range of an interpolation grid. 




        

<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/mapping/discrete_mapping_builder.hpp`

