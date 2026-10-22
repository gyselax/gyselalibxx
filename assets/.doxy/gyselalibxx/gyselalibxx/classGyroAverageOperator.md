

# Class GyroAverageOperator

**template &lt;concepts::Interpolation InterpolatorType, class IdxRangeParams, class IdxRangeBatch&gt;**



[**ClassList**](annotated.md) **>** [**GyroAverageOperator**](classGyroAverageOperator.md)



_Operator to compute the gyroaverage of a field as a weighted sum of the its value at pre- determined sample points along the gyro-motion._ [More...](#detailed-description)

* `#include <gyroaverage_operator.hpp>`

















## Public Types

| Type | Name |
| ---: | :--- |
| typedef Field&lt; CoordsInterpSpace, IdxRangeGyroCoeffs &gt; | [**GyroCoeffsCoords**](#typedef-gyrocoeffscoords)  <br>_A Field of point coordinates for the sample points._  |
| typedef FieldMem&lt; CoordsInterpSpace, IdxRangeGyroCoeffs &gt; | [**GyroCoeffsCoordsMem**](#typedef-gyrocoeffscoordsmem)  <br>_A FieldMem of point coordinates for the sample points._  |
| typedef DField&lt; IdxRangeGyroCoeffs &gt; | [**GyroCoeffsWeights**](#typedef-gyrocoeffsweights)  <br>_A Field of doubles containing the sample points' weights._  |
| typedef DFieldMem&lt; IdxRangeGyroCoeffs &gt; | [**GyroCoeffsWeightsMem**](#typedef-gyrocoeffsweightsmem)  <br>_A FieldMem of doubles containing the sample points' weights._  |




















## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**GyroAverageOperator**](#function-gyroaverageoperator) (InterpolatorType const & interpolator, [**GyroCoeffsCoordsMem**](classGyroAverageOperator.md#typedef-gyrocoeffscoordsmem) & gyrocoeffs\_coords, [**GyroCoeffsWeightsMem**](classGyroAverageOperator.md#typedef-gyrocoeffsweightsmem) & gyrocoeffs\_weights) <br>_Constructor for the general case where coefficients are computed by the caller._  |
|  void | [**operator()**](#function-operator) (DField&lt; IdxRangeGyroBatchedSpace &gt; const & func) const<br>_Compute the gyroaverage of a field which depends both on space and parameter sets._  |
|  void | [**operator()**](#function-operator_1) (DField&lt; IdxRangeInterpBatchedSpace &gt; const & func, IdxParams idx\_param) const<br>_For a given parameter set, compute the gyroaverage of a field which only depends on space (and on the batch dimensions)._  |


## Public Static Functions

| Type | Name |
| ---: | :--- |
|  std::tuple&lt; [**GyroCoeffsCoordsMem**](classGyroAverageOperator.md#typedef-gyrocoeffscoordsmem), [**GyroCoeffsWeightsMem**](classGyroAverageOperator.md#typedef-gyrocoeffsweightsmem) &gt; | [**build\_circular\_coefficients**](#function-build_circular_coefficients) (TransformFromCartesian transform\_from\_cartesian, DField&lt; IdxRangeParams &gt; const larmor\_factor\_params, DField&lt; IdxRangeInterpSpace &gt; const larmor\_factor\_space, size\_t const nb\_gyro\_points=8) <br>_Helper function building coefficients for a given set of Larmor radii, in geometries where an inverse mapping to Cartesian coordinates is available._  |
|  std::tuple&lt; [**GyroCoeffsCoordsMem**](classGyroAverageOperator.md#typedef-gyrocoeffscoordsmem), [**GyroCoeffsWeightsMem**](classGyroAverageOperator.md#typedef-gyrocoeffsweightsmem) &gt; | [**build\_empty\_coefficients**](#function-build_empty_coefficients) (IdxRangeParams idx\_range\_params, IdxRangeInterpSpace idx\_range\_space, size\_t const nb\_gyro\_points=8) <br>_Helper function allocating empty coefficients of the right size and type._  |


























## Detailed Description


Interpolation is used to evaluate the field at points which are not on the grid. The choice of treatment for points outside the domain (Dirichlet, constant outside, etc) is determined by the interpolation class itself.


We distinguish 3 types of dimensions:
* "Params" refers to the dimensions which parametrise the gyroaverage of a given field; this usually corresponds to species and mu.
* "InterpSpace" refers to the space dimensions along which the interpolation is performed, e.g. poloidal cross-sections when the approximation B // e\_phi is assumed.
* "BatchSpace" refers to the remaining dimensions, if any. The same coefficients are used for every batch index.






**Template parameters:**


* `InterpolatorType` Type of the interpolator. 
* `IdxRangeParams` The additional dimensions parametrising the gyroaverage (usually species and mu). 
* `IdxRangeBatch` The batch dimensions on which the gyroaverage does not depend (default: none). 




    
## Public Types Documentation




### typedef GyroCoeffsCoords 

_A Field of point coordinates for the sample points._ 
```C++
using GyroAverageOperator< InterpolatorType, IdxRangeParams, IdxRangeBatch >::GyroCoeffsCoords =  Field<CoordsInterpSpace, IdxRangeGyroCoeffs>;
```




<hr>



### typedef GyroCoeffsCoordsMem 

_A FieldMem of point coordinates for the sample points._ 
```C++
using GyroAverageOperator< InterpolatorType, IdxRangeParams, IdxRangeBatch >::GyroCoeffsCoordsMem =  FieldMem<CoordsInterpSpace, IdxRangeGyroCoeffs>;
```




<hr>



### typedef GyroCoeffsWeights 

_A Field of doubles containing the sample points' weights._ 
```C++
using GyroAverageOperator< InterpolatorType, IdxRangeParams, IdxRangeBatch >::GyroCoeffsWeights =  DField<IdxRangeGyroCoeffs>;
```




<hr>



### typedef GyroCoeffsWeightsMem 

_A FieldMem of doubles containing the sample points' weights._ 
```C++
using GyroAverageOperator< InterpolatorType, IdxRangeParams, IdxRangeBatch >::GyroCoeffsWeightsMem =  DFieldMem<IdxRangeGyroCoeffs>;
```




<hr>
## Public Functions Documentation




### function GyroAverageOperator 

_Constructor for the general case where coefficients are computed by the caller._ 
```C++
inline explicit GyroAverageOperator::GyroAverageOperator (
    InterpolatorType const & interpolator,
    GyroCoeffsCoordsMem & gyrocoeffs_coords,
    GyroCoeffsWeightsMem & gyrocoeffs_weights
) 
```





**Parameters:**


* `interpolator` Interpolator for the evaluation at the sample points. 
* `gyrocoeffs_coords` Coordinates of the sample points. 
* `gyrocoeffs_weights` Weights associated with the sample points. 




        

<hr>



### function operator() 

_Compute the gyroaverage of a field which depends both on space and parameter sets._ 
```C++
inline void GyroAverageOperator::operator() (
    DField< IdxRangeGyroBatchedSpace > const & func
) const
```



The computation is done in-place for each parameter set.




**Parameters:**


* `func` Field to process on input and its gyroaverage on output. 




        

<hr>



### function operator() 

_For a given parameter set, compute the gyroaverage of a field which only depends on space (and on the batch dimensions)._ 
```C++
inline void GyroAverageOperator::operator() (
    DField< IdxRangeInterpBatchedSpace > const & func,
    IdxParams idx_param
) const
```





**Parameters:**


* `func` Field to process on input and its gyroaverage on output. 
* `idx_param` Parameter set whose coefficients will be used. 




        

<hr>
## Public Static Functions Documentation




### function build\_circular\_coefficients 

_Helper function building coefficients for a given set of Larmor radii, in geometries where an inverse mapping to Cartesian coordinates is available._ 
```C++
template<concepts::AnalyticalMapping TransformFromCartesian>
static inline std::tuple< GyroCoeffsCoordsMem , GyroCoeffsWeightsMem > GyroAverageOperator::build_circular_coefficients (
    TransformFromCartesian transform_from_cartesian,
    DField< IdxRangeParams > const larmor_factor_params,
    DField< IdxRangeInterpSpace > const larmor_factor_space,
    size_t const nb_gyro_points=8
) 
```



The Larmor radii are of the form `larmor_factor_params * larmor_factor_space`.




**Parameters:**


* `transform_from_cartesian` Transformation from Cartesian to mesh coordinates. 
* `larmor_factor_params` Factor of the Larmor radius depending on the parameter sets. 
* `larmor_factor_space` Factor of the Larmor radius depending on space. 
* `nb_gyro_points` Number of sample points on the gyro-trajectory (default: 8).



**Returns:**

A tuple containing the point coordinates and their weights that were computed. 





        

<hr>



### function build\_empty\_coefficients 

_Helper function allocating empty coefficients of the right size and type._ 
```C++
static inline std::tuple< GyroCoeffsCoordsMem , GyroCoeffsWeightsMem > GyroAverageOperator::build_empty_coefficients (
    IdxRangeParams idx_range_params,
    IdxRangeInterpSpace idx_range_space,
    size_t const nb_gyro_points=8
) 
```





**Parameters:**


* `idx_range_params` The index range for the parameter sets. 
* `idx_range_space` The index range for the interpolation space. 
* `nb_gyro_points` The number of sample points on the gyro-trajectory (default: 8).



**Returns:**

A tuple containing fields for the coordinates and their weights. 





        

<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/gyroaverage/gyroaverage_operator.hpp`

