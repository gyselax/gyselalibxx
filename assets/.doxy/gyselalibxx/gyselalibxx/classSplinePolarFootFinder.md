

# Class SplinePolarFootFinder

**template &lt;class TimeStepper, class LogicalToPhysicalMapping, class LogicalToPseudoPhysicalMapping, class SplineRThetaBuilderAdvection, class SplineRThetaEvaluatorAdvection&gt;**



[**ClassList**](annotated.md) **>** [**SplinePolarFootFinder**](classSplinePolarFootFinder.md)



_A class to find the foot of the characteristics on the_  _plane._[More...](#detailed-description)

* `#include <spline_polar_foot_finder.hpp>`



Inherits the following classes: [IPolarFootFinder](classIPolarFootFinder.md)














## Public Types

| Type | Name |
| ---: | :--- |
| typedef ConstField&lt; CoordRTheta, [**IdxRangeOperator**](classSplinePolarFootFinder.md#typedef-idxrangeoperator), [**memory\_space**](classSplinePolarFootFinder.md#typedef-memory_space) &gt; | [**CConstFieldFeet**](#typedef-cconstfieldfeet)  <br>_The type of a constant field of (r, theta) coordinates at every grid point, saved on a compatible memory space._  |
| typedef Field&lt; CoordRTheta, [**IdxRangeOperator**](classSplinePolarFootFinder.md#typedef-idxrangeoperator), [**memory\_space**](classSplinePolarFootFinder.md#typedef-memory_space) &gt; | [**CFieldFeet**](#typedef-cfieldfeet)  <br>_The type of a field of (r, theta) coordinates at every grid point, saved on a compatible memory space._  |
| typedef [**DVectorConstField**](classVectorField.md)&lt; [**IdxRangeOperator**](classSplinePolarFootFinder.md#typedef-idxrangeoperator), PseudoCartesianBasis, [**memory\_space**](classSplinePolarFootFinder.md#typedef-memory_space) &gt; | [**DVectorConstFieldAdvection**](#typedef-dvectorconstfieldadvection)  <br>_The type of a constant vector field defined on the pseudo-Cartesian basis at every grid point, saved on a compatible memory space._  |
| typedef [**DVectorField**](classVectorField.md)&lt; [**IdxRangeOperator**](classSplinePolarFootFinder.md#typedef-idxrangeoperator), PseudoCartesianBasis, [**memory\_space**](classSplinePolarFootFinder.md#typedef-memory_space) &gt; | [**DVectorFieldAdvection**](#typedef-dvectorfieldadvection)  <br>_The type of a vector field defined on the pseudo-Cartesian basis at every grid point, saved on a compatible memory space._  |
| typedef GridRadial | [**GridR**](#typedef-gridr)  <br>_The continuous radial dimension._  |
| typedef GridPoloidal | [**GridTheta**](#typedef-gridtheta)  <br>_The continuous poloidal dimension._  |
| typedef IdxRangeBatched | [**IdxRangeOperator**](#typedef-idxrangeoperator)  <br>_The type of the index range over which the operator works._  |
| typedef typename GridR::continuous\_dimension\_type | [**R**](#typedef-r)  <br>_The continuous radial dimension._  |
| typedef typename GridTheta::continuous\_dimension\_type | [**Theta**](#typedef-theta)  <br>_The continuous poloidal dimension._  |
| typedef VectorIndexSetAdvDims | [**VectorIndexSetAdvectionDims**](#typedef-vectorindexsetadvectiondims)  <br>_The continuous radial dimension._  |
| typedef [**DVectorFieldMem**](classVectorFieldMem.md)&lt; IdxRangeSplineBatched, PseudoCartesianBasis, [**memory\_space**](classSplinePolarFootFinder.md#typedef-memory_space) &gt; | [**VectorSplineCoeffsMem**](#typedef-vectorsplinecoeffsmem)  <br>_The type of 2 batched splines representing the x and y components of a vector on the polar plane on a compatible memory space._  |
| typedef MemorySpace | [**memory\_space**](#typedef-memory_space)  <br>_The type of the memory space where the field is saved (CPU vs GPU)._  |


## Public Types inherited from IPolarFootFinder

See [IPolarFootFinder](classIPolarFootFinder.md)

| Type | Name |
| ---: | :--- |
| typedef IdxRangeBatched | [**IdxRangeOperator**](classIPolarFootFinder.md#typedef-idxrangeoperator)  <br>_The type of the index range over which the operator works._  |
| typedef MemorySpace | [**memory\_space**](classIPolarFootFinder.md#typedef-memory_space)  <br>_The type of the memory space where the field is saved (CPU vs GPU)._  |






































## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**SplinePolarFootFinder**](#function-splinepolarfootfinder) (TimeStepper const & time\_stepper, LogicalToPhysicalMapping const & logical\_to\_physical\_mapping, LogicalToPseudoPhysicalMapping const & logical\_to\_pseudo\_physical\_mapping, SplineRThetaBuilderAdvection const & builder\_advection\_field, SplineRThetaEvaluatorAdvection const & evaluator\_advection\_field, double epsilon=1e-12) <br>_Instantiate a time integration method for the advection operator._  |
|  void | [**is\_unified**](#function-is_unified) (Field&lt; [**T**](structT.md), IdxRangeRTheta, [**memory\_space**](classSplinePolarFootFinder.md#typedef-memory_space) &gt; const & values) const<br>_Check if the values at the centre point are the same._  |
|  void | [**operator()**](#function-operator) ([**CFieldFeet**](classSplinePolarFootFinder.md#typedef-cfieldfeet) feet, [**DVectorConstField**](classVectorField.md)&lt; [**IdxRangeOperator**](classSplinePolarFootFinder.md#typedef-idxrangeoperator), [**VectorIndexSetAdvectionDims**](classSplinePolarFootFinder.md#typedef-vectorindexsetadvectiondims), [**memory\_space**](classSplinePolarFootFinder.md#typedef-memory_space) &gt; advection\_field, double dt) const<br>_Advect the feet over_  _._ |
|  void | [**unify\_value\_at\_centre\_pt**](#function-unify_value_at_centre_pt) (Field&lt; [**T**](structT.md), IdxRangeRTheta, [**memory\_space**](classSplinePolarFootFinder.md#typedef-memory_space) &gt; values) const<br>_Replace the value at_  _point by the value at_ _for all_ _._ |


## Public Functions inherited from IPolarFootFinder

See [IPolarFootFinder](classIPolarFootFinder.md)

| Type | Name |
| ---: | :--- |
| virtual void | [**operator()**](classIPolarFootFinder.md#function-operator) (Field&lt; Coord&lt; [**R**](classIPolarFootFinder.md#typedef-r), [**Theta**](classIPolarFootFinder.md#typedef-theta) &gt;, [**IdxRangeOperator**](classIPolarFootFinder.md#typedef-idxrangeoperator), [**memory\_space**](classIPolarFootFinder.md#typedef-memory_space) &gt; feet, [**DVectorConstField**](classVectorField.md)&lt; [**IdxRangeOperator**](classIPolarFootFinder.md#typedef-idxrangeoperator), [**VectorIndexSetAdvectionDims**](classIPolarFootFinder.md#typedef-vectorindexsetadvectiondims), [**memory\_space**](classIPolarFootFinder.md#typedef-memory_space) &gt; advection\_field, double dt) const = 0<br>_Advect the feet over_  _._ |
| virtual  | [**~IPolarFootFinder**](classIPolarFootFinder.md#function-ipolarfootfinder) () = default<br> |








## Protected Types inherited from IPolarFootFinder

See [IPolarFootFinder](classIPolarFootFinder.md)

| Type | Name |
| ---: | :--- |
| typedef GridRadial | [**GridR**](classIPolarFootFinder.md#typedef-gridr)  <br>_The continuous radial dimension._  |
| typedef GridPoloidal | [**GridTheta**](classIPolarFootFinder.md#typedef-gridtheta)  <br>_The continuous poloidal dimension._  |
| typedef typename GridR::continuous\_dimension\_type | [**R**](classIPolarFootFinder.md#typedef-r)  <br>_The continuous radial dimension._  |
| typedef typename GridTheta::continuous\_dimension\_type | [**Theta**](classIPolarFootFinder.md#typedef-theta)  <br>_The continuous poloidal dimension._  |
| typedef VectorIndexSetAdvDims | [**VectorIndexSetAdvectionDims**](classIPolarFootFinder.md#typedef-vectorindexsetadvectiondims)  <br>_The continuous radial dimension._  |














































## Detailed Description


The natural advection domain is the physical domain, where the studied equation is given. However, not all the mappings used are analytically invertible and inverting the Jacobian matrix of the mapping could be costly and could introduce numerical errors. That is why, we also introduce a pseudo-Cartesian domain.


More details can be found in Edoardo Zoni's article ([https://doi.org/10.1016/j.jcp.2019.108889](https://doi.org/10.1016/j.jcp.2019.108889)).




**Template parameters:**


* `TimeStepper` A child class of [**ITimeStepper**](classITimeStepper.md) providing a time integration method. 
* `LogicalToPhysicalMapping` A mapping from the logical domain to the physical domain. 
* `LogicalToPseudoPhysicalMapping` A mapping from the logical domain to the domain where the advection is carried out. This may be a pseudo-physical domain or the physical domain itself. 
* `SplineRThetaBuilderAdvection` A 2D SplineBuilder to construct a spline on a polar domain. 
* `SplineRThetaEvaluatorAdvection` A 2D SplineEvaluator to evaluate a spline on a polar domain. A boundary condition must be provided in case the foot of the characteristic is found outside the domain.



**See also:** [**BslAdvectionPolar**](classBslAdvectionPolar.md) 



    
## Public Types Documentation




### typedef CConstFieldFeet 

_The type of a constant field of (r, theta) coordinates at every grid point, saved on a compatible memory space._ 
```C++
using SplinePolarFootFinder< TimeStepper, LogicalToPhysicalMapping, LogicalToPseudoPhysicalMapping, SplineRThetaBuilderAdvection, SplineRThetaEvaluatorAdvection >::CConstFieldFeet =  ConstField<CoordRTheta, IdxRangeOperator, memory_space>;
```




<hr>



### typedef CFieldFeet 

_The type of a field of (r, theta) coordinates at every grid point, saved on a compatible memory space._ 
```C++
using SplinePolarFootFinder< TimeStepper, LogicalToPhysicalMapping, LogicalToPseudoPhysicalMapping, SplineRThetaBuilderAdvection, SplineRThetaEvaluatorAdvection >::CFieldFeet =  Field<CoordRTheta, IdxRangeOperator, memory_space>;
```




<hr>



### typedef DVectorConstFieldAdvection 

_The type of a constant vector field defined on the pseudo-Cartesian basis at every grid point, saved on a compatible memory space._ 
```C++
using SplinePolarFootFinder< TimeStepper, LogicalToPhysicalMapping, LogicalToPseudoPhysicalMapping, SplineRThetaBuilderAdvection, SplineRThetaEvaluatorAdvection >::DVectorConstFieldAdvection =  DVectorConstField<IdxRangeOperator, PseudoCartesianBasis, memory_space>;
```




<hr>



### typedef DVectorFieldAdvection 

_The type of a vector field defined on the pseudo-Cartesian basis at every grid point, saved on a compatible memory space._ 
```C++
using SplinePolarFootFinder< TimeStepper, LogicalToPhysicalMapping, LogicalToPseudoPhysicalMapping, SplineRThetaBuilderAdvection, SplineRThetaEvaluatorAdvection >::DVectorFieldAdvection =  DVectorField<IdxRangeOperator, PseudoCartesianBasis, memory_space>;
```




<hr>



### typedef GridR 

_The continuous radial dimension._ 
```C++
using IPolarFootFinder< GridRadial, GridPoloidal, VectorIndexSetAdvDims, IdxRangeBatched, MemorySpace >::GridR =  GridRadial;
```




<hr>



### typedef GridTheta 

_The continuous poloidal dimension._ 
```C++
using IPolarFootFinder< GridRadial, GridPoloidal, VectorIndexSetAdvDims, IdxRangeBatched, MemorySpace >::GridTheta =  GridPoloidal;
```




<hr>



### typedef IdxRangeOperator 

_The type of the index range over which the operator works._ 
```C++
using IPolarFootFinder< GridRadial, GridPoloidal, VectorIndexSetAdvDims, IdxRangeBatched, MemorySpace >::IdxRangeOperator =  IdxRangeBatched;
```




<hr>



### typedef R 

_The continuous radial dimension._ 
```C++
using IPolarFootFinder< GridRadial, GridPoloidal, VectorIndexSetAdvDims, IdxRangeBatched, MemorySpace >::R =  typename GridR::continuous_dimension_type;
```




<hr>



### typedef Theta 

_The continuous poloidal dimension._ 
```C++
using IPolarFootFinder< GridRadial, GridPoloidal, VectorIndexSetAdvDims, IdxRangeBatched, MemorySpace >::Theta =  typename GridTheta::continuous_dimension_type;
```




<hr>



### typedef VectorIndexSetAdvectionDims 

_The continuous radial dimension._ 
```C++
using IPolarFootFinder< GridRadial, GridPoloidal, VectorIndexSetAdvDims, IdxRangeBatched, MemorySpace >::VectorIndexSetAdvectionDims =  VectorIndexSetAdvDims;
```




<hr>



### typedef VectorSplineCoeffsMem 

_The type of 2 batched splines representing the x and y components of a vector on the polar plane on a compatible memory space._ 
```C++
using SplinePolarFootFinder< TimeStepper, LogicalToPhysicalMapping, LogicalToPseudoPhysicalMapping, SplineRThetaBuilderAdvection, SplineRThetaEvaluatorAdvection >::VectorSplineCoeffsMem =  DVectorFieldMem<IdxRangeSplineBatched, PseudoCartesianBasis, memory_space>;
```




<hr>



### typedef memory\_space 

_The type of the memory space where the field is saved (CPU vs GPU)._ 
```C++
using IPolarFootFinder< GridRadial, GridPoloidal, VectorIndexSetAdvDims, IdxRangeBatched, MemorySpace >::memory_space =  MemorySpace;
```




<hr>
## Public Functions Documentation




### function SplinePolarFootFinder 

_Instantiate a time integration method for the advection operator._ 
```C++
inline SplinePolarFootFinder::SplinePolarFootFinder (
    TimeStepper const & time_stepper,
    LogicalToPhysicalMapping const & logical_to_physical_mapping,
    LogicalToPseudoPhysicalMapping const & logical_to_pseudo_physical_mapping,
    SplineRThetaBuilderAdvection const & builder_advection_field,
    SplineRThetaEvaluatorAdvection const & evaluator_advection_field,
    double epsilon=1e-12
) 
```





**Parameters:**


* `time_stepper` The time integration method used to solve the characteristic equation ([**ITimeStepper**](classITimeStepper.md)). 
* `logical_to_physical_mapping` The mapping from the logical domain to the physical domain. 
* `logical_to_pseudo_physical_mapping` The mapping from the logical domain to the pseudo-physical domain. 
* `builder_advection_field` The spline builder which computes the spline representation of the advection field. 
* `evaluator_advection_field` The B-splines evaluator to evaluate the advection field. 
* `epsilon`  parameter used for the linearisation of the advection field around the central point.



**See also:** [**ITimeStepper**](classITimeStepper.md) 



        

<hr>



### function is\_unified 

_Check if the values at the centre point are the same._ 
```C++
template<class T>
inline void SplinePolarFootFinder::is_unified (
    Field< T , IdxRangeRTheta, memory_space > const & values
) const
```



For polar geometry, to ensure continuity at the centre point, we have to be sure that all the points for  have the same value. This function check if for , the values  are the same.




**Parameters:**


* `values` A table of values we want to check if the centre point has an unique value. 




        

<hr>



### function operator() 

_Advect the feet over_  _._
```C++
inline void SplinePolarFootFinder::operator() (
    CFieldFeet feet,
    DVectorConstField < IdxRangeOperator , VectorIndexSetAdvectionDims , memory_space > advection_field,
    double dt
) const
```



From the advection field in the physical index range, compute the advection field in the right index range an compute its B-splines coefficients. Then, use the given time integration method (time\_stepper) to solve the characteristic equation over .




**Parameters:**


* `feet` On input: the mesh points. On output: the characteristic feet. 
* `advection_field` The advection field in the physical index range. 
* `dt` The time step. 




        

<hr>



### function unify\_value\_at\_centre\_pt 

_Replace the value at_  _point by the value at_ _for all_ _._
```C++
template<class T>
inline void SplinePolarFootFinder::unify_value_at_centre_pt (
    Field< T , IdxRangeRTheta, memory_space > values
) const
```



For polar geometry, to ensure continuity at the centre point, we have to be sure that all the points for  have the same value. As the computation of the values of a table can induces machine errors, this function is useful to reset the values at the central point at the same value.




**Parameters:**


* `values` The table of values we want to unify at the central point. 




        

<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/advection/spline_polar_foot_finder.hpp`

