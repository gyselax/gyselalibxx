

# Class RefinedDiscretePoloidalCSSplineMappingBuilder

**template &lt;class [**X**](structX.md), class [**Y**](structY.md), class SplineBuilder, class SplineEvaluator, int ncells\_r, int ncells\_theta&gt;**



[**ClassList**](annotated.md) **>** [**RefinedDiscretePoloidalCSSplineMappingBuilder**](classRefinedDiscretePoloidalCSSplineMappingBuilder.md)



_A class to create a_ [_**DiscretePoloidalCSSplineMapping**_](classDiscretePoloidalCSSplineMapping.md) _instance from an analytical mapping. This class creates an instance which uses more refined splines than the provided builder and evaluator. This class creates and stores splines memory spaces describing the analytical mapping. The discrete mapping is then created using the splines without copying data._[More...](#detailed-description)

* `#include <discrete_poloidal_cs_spline_mapping_builder.hpp>`















## Classes

| Type | Name |
| ---: | :--- |
| struct | [**BSplinesRRefined**](structRefinedDiscretePoloidalCSSplineMappingBuilder_1_1BSplinesRRefined.md) <br>_The type of the radial B-splines on which the new mapping will be defined._  |
| struct | [**BSplinesThetaRefined**](structRefinedDiscretePoloidalCSSplineMappingBuilder_1_1BSplinesThetaRefined.md) <br>_The type of the poloidal B-splines on which the new mapping will be defined._  |
| struct | [**GridRRefined**](structRefinedDiscretePoloidalCSSplineMappingBuilder_1_1GridRRefined.md) <br>_The type of the grid of radial points on which the new mapping will be defined._  |
| struct | [**GridThetaRefined**](structRefinedDiscretePoloidalCSSplineMappingBuilder_1_1GridThetaRefined.md) <br>_The type of the grid of poloidal points on which the new mapping will be defined._  |


## Public Types

| Type | Name |
| ---: | :--- |
| typedef [**DiscretePoloidalCSSplineMapping**](classDiscretePoloidalCSSplineMapping.md)&lt; [**X**](structX.md), [**Y**](structY.md), RefinedSplineEvaluator &gt; | [**MappingType**](#typedef-mappingtype)  <br>_The type of the mapping that will be created._  |




















## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**RefinedDiscretePoloidalCSSplineMappingBuilder**](#function-refineddiscretepoloidalcssplinemappingbuilder) (ExecSpace exec\_space, Mapping const & analytical\_mapping, SplineBuilder const & builder, SplineEvaluator const & evaluator) <br>_Create an instance of the class capable of providing a_ [_**DiscretePoloidalCSSplineMapping**_](classDiscretePoloidalCSSplineMapping.md) _class instance._ |
|  [**DiscretePoloidalCSSplineMapping**](classDiscretePoloidalCSSplineMapping.md)&lt; [**X**](structX.md), [**Y**](structY.md), RefinedSplineEvaluator &gt; | [**operator()**](#function-operator) () const<br>_Get a_ [_**DiscretePoloidalCSSplineMapping**_](classDiscretePoloidalCSSplineMapping.md) _class instance._ |
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
using RefinedDiscretePoloidalCSSplineMappingBuilder< X, Y, SplineBuilder, SplineEvaluator, ncells_r, ncells_theta >::MappingType =  DiscretePoloidalCSSplineMapping<X, Y, RefinedSplineEvaluator>;
```




<hr>
## Public Functions Documentation




### function RefinedDiscretePoloidalCSSplineMappingBuilder 

_Create an instance of the class capable of providing a_ [_**DiscretePoloidalCSSplineMapping**_](classDiscretePoloidalCSSplineMapping.md) _class instance._
```C++
template<class Mapping>
inline RefinedDiscretePoloidalCSSplineMappingBuilder::RefinedDiscretePoloidalCSSplineMappingBuilder (
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

_Get a_ [_**DiscretePoloidalCSSplineMapping**_](classDiscretePoloidalCSSplineMapping.md) _class instance._
```C++
inline DiscretePoloidalCSSplineMapping < X , Y , RefinedSplineEvaluator > RefinedDiscretePoloidalCSSplineMappingBuilder::operator() () const
```





**Returns:**

An instance of the mapping. 





        

<hr>



### function set\_curvilinear\_to\_cartesian\_values 

_Fill in the curvilinear fields with interpolation points mapped with the given analytical mapping._ 
```C++
template<class Mapping>
inline void RefinedDiscretePoloidalCSSplineMappingBuilder::set_curvilinear_to_cartesian_values (
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
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/coord_transformations/discrete_poloidal_cs_spline_mapping_builder.hpp`

