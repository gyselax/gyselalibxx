

# Class RefinedDiscreteToCartesianBuilder

**template &lt;class [**X**](structX.md), class [**Y**](structY.md), class SplineBuilder, class SplineEvaluator, int ncells\_r, int ncells\_theta&gt;**



[**ClassList**](annotated.md) **>** [**RefinedDiscreteToCartesianBuilder**](classRefinedDiscreteToCartesianBuilder.md)



_A class to create a_ [_**DiscreteToCartesian**_](classDiscreteToCartesian.md) _instance from an analytical mapping. This class creates an instance which uses more refined splines than the provided builder and evaluator. This class creates and stores splines memory spaces describing the analytical mapping. The discrete mapping is then created using the splines without copying data._[More...](#detailed-description)

* `#include <discrete_mapping_builder.hpp>`















## Classes

| Type | Name |
| ---: | :--- |
| struct | [**BSplinesRRefined**](structRefinedDiscreteToCartesianBuilder_1_1BSplinesRRefined.md) <br>_The type of the radial B-splines on which the new mapping will be defined._  |
| struct | [**BSplinesThetaRefined**](structRefinedDiscreteToCartesianBuilder_1_1BSplinesThetaRefined.md) <br>_The type of the poloidal B-splines on which the new mapping will be defined._  |
| struct | [**GridRRefined**](structRefinedDiscreteToCartesianBuilder_1_1GridRRefined.md) <br>_The type of the grid of radial points on which the new mapping will be defined._  |
| struct | [**GridThetaRefined**](structRefinedDiscreteToCartesianBuilder_1_1GridThetaRefined.md) <br>_The type of the grid of poloidal points on which the new mapping will be defined._  |


## Public Types

| Type | Name |
| ---: | :--- |
| typedef [**DiscreteToCartesian**](classDiscreteToCartesian.md)&lt; [**X**](structX.md), [**Y**](structY.md), RefinedSplineEvaluator &gt; | [**MappingType**](#typedef-mappingtype)  <br>_The type of the mapping that will be created._  |




















## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**RefinedDiscreteToCartesianBuilder**](#function-refineddiscretetocartesianbuilder) (ExecSpace exec\_space, Mapping const & analytical\_mapping, SplineBuilder const & builder, SplineEvaluator const & evaluator) <br>_Create an instance of the class capable of providing a_ [_**DiscreteToCartesian**_](classDiscreteToCartesian.md) _class instance._ |
|  [**DiscreteToCartesian**](classDiscreteToCartesian.md)&lt; [**X**](structX.md), [**Y**](structY.md), RefinedSplineEvaluator &gt; | [**operator()**](#function-operator) () const<br>_Get a_ [_**DiscreteToCartesian**_](classDiscreteToCartesian.md) _class instance._ |
|  void | [**set\_curvilinear\_to\_cartesian\_values**](#function-set_curvilinear_to_cartesian_values) (InterpolationField curvilinear\_to\_x\_vals, InterpolationField curvilinear\_to\_y\_vals, Mapping const & analytical\_mapping, IdxRangeInterpolationPoints const & interpolation\_idx\_range) <br>_Fill in the curvilinear fields with interpolation points mapped with the given analytical mapping._  |




























## Detailed Description




**Template parameters:**


* [**X**](structX.md) The first Cartesian dimension. 
* [**Y**](structY.md) The second Cartesian dimension. 
* `SplineBuilder` An operator for building spline coefficients. 
* `SplineEvaluator` An operator for evaluating a spline. 
* `ncells_r` The number of cells in the refined spline in the radial direction. 
* `ncells_theta` The number of cells in the refined spline in the radial direction. 




    
## Public Types Documentation




### typedef MappingType 

_The type of the mapping that will be created._ 
```C++
using RefinedDiscreteToCartesianBuilder< X, Y, SplineBuilder, SplineEvaluator, ncells_r, ncells_theta >::MappingType =  DiscreteToCartesian<X, Y, RefinedSplineEvaluator>;
```




<hr>
## Public Functions Documentation




### function RefinedDiscreteToCartesianBuilder 

_Create an instance of the class capable of providing a_ [_**DiscreteToCartesian**_](classDiscreteToCartesian.md) _class instance._
```C++
template<class Mapping>
inline RefinedDiscreteToCartesianBuilder::RefinedDiscreteToCartesianBuilder (
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
inline DiscreteToCartesian < X , Y , RefinedSplineEvaluator > RefinedDiscreteToCartesianBuilder::operator() () const
```





**Returns:**

An instance of the mapping. 





        

<hr>



### function set\_curvilinear\_to\_cartesian\_values 

_Fill in the curvilinear fields with interpolation points mapped with the given analytical mapping._ 
```C++
template<class Mapping>
inline void RefinedDiscreteToCartesianBuilder::set_curvilinear_to_cartesian_values (
    InterpolationField curvilinear_to_x_vals,
    InterpolationField curvilinear_to_y_vals,
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

