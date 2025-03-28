

# Class ddcHelper::NonUniformInterpolationPoints

**template &lt;class BSplines, ddc::BoundCond BcXmin, ddc::BoundCond BcXmax&gt;**



[**ClassList**](annotated.md) **>** [**ddcHelper**](namespaceddcHelper.md) **>** [**NonUniformInterpolationPoints**](classddcHelper_1_1NonUniformInterpolationPoints.md)



_Helper class for the initialisation of the mesh of interpolation points._ [More...](#detailed-description)

* `#include <non_uniform_interpolation_points.hpp>`

















## Public Types

| Type | Name |
| ---: | :--- |
| typedef NonUniformGridBase&lt; Dim &gt; | [**interpolation\_discrete\_dimension\_type**](#typedef-interpolation_discrete_dimension_type)  <br> |






## Public Static Attributes

| Type | Name |
| ---: | :--- |
|  constexpr int | [**N\_BE\_MAX**](#variable-n_be_max)   = `n\_boundary\_equations(BcXmax, BSplines::degree())`<br>_The number of boundary equations at the upper bound in a spline interpolation._  |
|  constexpr int | [**N\_BE\_MIN**](#variable-n_be_min)   = `n\_boundary\_equations(BcXmin, BSplines::degree())`<br>_The number of boundary equations at the lower bound in a spline interpolation._  |
















## Public Static Functions

| Type | Name |
| ---: | :--- |
|  IdxRange&lt; Sampling &gt; | [**get\_domain**](#function-get_domain) () <br> |
|  auto | [**get\_sampling**](#function-get_sampling) (std::vector&lt; Coord&lt; Dim &gt; &gt; & interp\_points, double TOL=2e-14) <br> |


























## Detailed Description


A helper class for the initialisation of the mesh of interpolation points. This class should be used when the interpolation points are specified externally. This is possible with any kind of boundary condition. In the case of Greville boundary conditions, the second and second to last interpolation points should be the Greville points. 


    
## Public Types Documentation




### typedef interpolation\_discrete\_dimension\_type 

```C++
using ddcHelper::NonUniformInterpolationPoints< BSplines, BcXmin, BcXmax >::interpolation_discrete_dimension_type =  NonUniformGridBase<Dim>;
```



The type of the mesh.


This is always NonUniformPointSampling. 


        

<hr>
## Public Static Attributes Documentation




### variable N\_BE\_MAX 

_The number of boundary equations at the upper bound in a spline interpolation._ 
```C++
constexpr int ddcHelper::NonUniformInterpolationPoints< BSplines, BcXmin, BcXmax >::N_BE_MAX;
```




<hr>



### variable N\_BE\_MIN 

_The number of boundary equations at the lower bound in a spline interpolation._ 
```C++
constexpr int ddcHelper::NonUniformInterpolationPoints< BSplines, BcXmin, BcXmax >::N_BE_MIN;
```




<hr>
## Public Static Functions Documentation




### function get\_domain 

```C++
template<typename Sampling>
static inline IdxRange< Sampling > ddcHelper::NonUniformInterpolationPoints::get_domain () 
```



Get the domain which can be used to access the interpolation points in the sampling.




**Returns:**

domain The discrete domain which maps to the sampling of interpolation points. 





        

<hr>



### function get\_sampling 

```C++
template<typename Sampling, typename U>
static inline auto ddcHelper::NonUniformInterpolationPoints::get_sampling (
    std::vector< Coord< Dim > > & interp_points,
    double TOL=2e-14
) 
```



Get the sampling of interpolation points.




**Parameters:**


* `interp_points` The grid points on which the simulation will evolve. 
* `TOL` The tolerance above which points are considered to be non-uniform.



**Returns:**

sampling The DDC point sampling of the interpolation points. 





        

<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/utils/non_uniform_interpolation_points.hpp`

