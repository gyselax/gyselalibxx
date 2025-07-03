

# Class GyroAverageOperator

**template &lt;class SplineRThetaBuilder, class SplineRThetaEvaluator, class IdxRangeRminorThetaBatch, class ToLogicalCoordTransform&gt;**



[**ClassList**](annotated.md) **>** [**GyroAverageOperator**](classGyroAverageOperator.md)



_Operator to compute the gyroaverage of a field in (r, theta) coordinates._ [More...](#detailed-description)

* `#include <gyroaverage_operator.hpp>`





































## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**GyroAverageOperator**](#function-gyroaverageoperator) (DConstFieldRminorTheta const & rho\_L, SplineRThetaBuilder const & spline\_builder, SplineRThetaEvaluator const & spline\_evaluator, ToLogicalCoordTransform coordinate\_transform, std::size\_t const nb\_gyro\_points=8) <br>_Constructor._  |
|  void | [**operator()**](#function-operator) (DFieldRminorThetaBatch const & A\_bar, DConstFieldRminorThetaBatch const & A) const<br>_Applies the gyroaverage operator to a batched field._  |




























## Detailed Description


This class performs the gyroaveraging operation on a batched field defined on a polar grid (r, theta). The gyroaverage is computed by integrating the field over a set of points along a circle (the Larmor orbit) centred at each grid point, with radius given by the local Larmor radius field (rho\_L).


The class uses 2D B-spline interpolation to evaluate the field at off-grid points along the orbit. The operation is performed in parallel over the (r, theta) grid.




**Template parameters:**


* `SplineRThetaBuilder` The type of the spline builder for the rtheta interpolation 
* `SplineRThetaEvaluator` The type of the spline evaluator for the rtheta interpolation 
* `IdxRangeRminorThetaBatch` The index range over [**R**](structR.md), [**Theta**](structTheta.md) and Batch directions. 
* `ToLogicalCoordTransform` Function to convert ([**R**](structR.md), Z) to (r, theta). 




    
## Public Functions Documentation




### function GyroAverageOperator 

_Constructor._ 
```C++
inline explicit GyroAverageOperator::GyroAverageOperator (
    DConstFieldRminorTheta const & rho_L,
    SplineRThetaBuilder const & spline_builder,
    SplineRThetaEvaluator const & spline_evaluator,
    ToLogicalCoordTransform coordinate_transform,
    std::size_t const nb_gyro_points=8
) 
```





**Parameters:**


* `rho_L` Field of Larmor radii on the (r, theta) grid. 
* `spline_builder` The spline builder for the rtheta interpolation 
* `spline_evaluator` The spline evaluator for the rtheta interpolation 
* `coordinate_transform` Function to convert ([**R**](structR.md), Z) to (r, theta). 
* `nb_gyro_points` Number of points to use in the gyroaverage integration (default: 8). 




        

<hr>



### function operator() 

_Applies the gyroaverage operator to a batched field._ 
```C++
inline void GyroAverageOperator::operator() (
    DFieldRminorThetaBatch const & A_bar,
    DConstFieldRminorThetaBatch const & A
) const
```



For each batch, and for each (r, theta) grid point, computes the gyroaverage by integrating the field along a circle of radius rho\_L centred at (r, theta). The field is interpolated at off-grid points using 2D B-splines.




**Parameters:**


* `A_bar` Output field to store the gyroaveraged result (batched). 
* `A` Input field to be gyroaveraged (batched). 




        

<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/gyroaverage/gyroaverage_operator.hpp`

