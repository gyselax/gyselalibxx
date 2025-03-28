

# Class SplinePolarFootFinder

**template &lt;class TimeStepper, class LogicalToPhysicalMapping, class LogicalToPseudoPhysicalMapping, class SplineRThetaBuilder, class SplineRThetaEvaluatorConstBound&gt;**



[**ClassList**](annotated.md) **>** [**SplinePolarFootFinder**](classSplinePolarFootFinder.md)



_A class to find the foot of the characteristics on the_  _plane._[More...](#detailed-description)

* `#include <spline_polar_foot_finder.hpp>`



Inherits the following classes: [IPolarFootFinder](classIPolarFootFinder.md)














## Public Types

| Type | Name |
| ---: | :--- |
| typedef ConstField&lt; ElementType, IdxRangeRTheta, memory\_space &gt; | [**ConstFieldRTheta**](#typedef-constfieldrtheta)  <br>_The type of a constant field on the polar plane on a compatible memory space._  |
| typedef [**VectorConstField**](classVectorField.md)&lt; double, IdxRangeRTheta, VectorIndexSet&lt; Dim1, Dim2 &gt;, memory\_space &gt; | [**DVectorConstFieldRTheta**](#typedef-dvectorconstfieldrtheta)  <br>_The type of a constant vector (x,y) field on the polar plane on a compatible memory space._  |
| typedef [**VectorField**](classVectorField.md)&lt; double, IdxRangeRTheta, VectorIndexSet&lt; Dim1, Dim2 &gt;, memory\_space &gt; | [**DVectorFieldRTheta**](#typedef-dvectorfieldrtheta)  <br>_The type of a vector (x,y) field on the polar plane on a compatible memory space._  |
| typedef Field&lt; ElementType, IdxRangeRTheta, memory\_space &gt; | [**FieldRTheta**](#typedef-fieldrtheta)  <br>_The type of a field on the polar plane on a compatible memory space._  |
| typedef [**VectorFieldMem**](classVectorFieldMem.md)&lt; double, IdxRange&lt; [**BSplinesR**](structBSplinesR.md), [**BSplinesTheta**](structBSplinesTheta.md) &gt;, VectorIndexSet&lt; Dim1, Dim2 &gt;, memory\_space &gt; | [**VectorSplineCoeffsMem2D**](#typedef-vectorsplinecoeffsmem2d)  <br>_The type of 2 splines representing the x and y components of a vector on the polar plane on a compatible memory space._  |








































## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**SplinePolarFootFinder**](#function-splinepolarfootfinder) (TimeStepper const & time\_stepper, LogicalToPhysicalMapping const & logical\_to\_physical\_mapping, LogicalToPseudoPhysicalMapping const & logical\_to\_pseudo\_physical\_mapping, SplineRThetaBuilder const & builder\_advection\_field, SplineRThetaEvaluatorConstBound const & evaluator\_advection\_field, double epsilon=1e-12) <br>_Instantiate a time integration method for the advection operator._  |
|  void | [**is\_unified**](#function-is_unified) (Field&lt; [**T**](structT.md), IdxRangeRTheta, memory\_space &gt; const & values) const<br>_Check if the values at the centre point are the same._  |
|  void | [**operator()**](#function-operator) ([**FieldRTheta**](classSplinePolarFootFinder.md#typedef-fieldrtheta)&lt; CoordRTheta &gt; feet, [**DVectorConstFieldRTheta**](classSplinePolarFootFinder.md#typedef-dvectorconstfieldrtheta)&lt; [**X**](structX.md), [**Y**](structY.md) &gt; advection\_field, double dt) const<br>_Advect the feet over_  _._ |
|  void | [**unify\_value\_at\_centre\_pt**](#function-unify_value_at_centre_pt) (Field&lt; [**T**](structT.md), IdxRangeRTheta, memory\_space &gt; values) const<br>_Replace the value at_  _point by the value at_ _for all_ _._ |


## Public Functions inherited from IPolarFootFinder

See [IPolarFootFinder](classIPolarFootFinder.md)

| Type | Name |
| ---: | :--- |
| virtual void | [**operator()**](classIPolarFootFinder.md#function-operator) (Field&lt; Coord&lt; [**R**](classIPolarFootFinder.md#typedef-r), [**Theta**](classIPolarFootFinder.md#typedef-theta) &gt;, [**IdxRangeRTheta**](classIPolarFootFinder.md#typedef-idxrangertheta), [**memory\_space**](classIPolarFootFinder.md#typedef-memory_space) &gt; feet, [**DVectorConstField**](classVectorField.md)&lt; [**IdxRangeRTheta**](classIPolarFootFinder.md#typedef-idxrangertheta), VectorIndexSet&lt; [**X**](classIPolarFootFinder.md#typedef-x), [**Y**](classIPolarFootFinder.md#typedef-y) &gt;, [**memory\_space**](classIPolarFootFinder.md#typedef-memory_space) &gt; advection\_field, double dt) const = 0<br>_Advect the feet over_  _._ |
| virtual  | [**~IPolarFootFinder**](classIPolarFootFinder.md#function-ipolarfootfinder) () = default<br> |








## Protected Types inherited from IPolarFootFinder

See [IPolarFootFinder](classIPolarFootFinder.md)

| Type | Name |
| ---: | :--- |
| typedef GridRadial | [**GridR**](classIPolarFootFinder.md#typedef-gridr)  <br>_The continuous radial dimension._  |
| typedef GridPoloidal | [**GridTheta**](classIPolarFootFinder.md#typedef-gridtheta)  <br>_The continuous poloidal dimension._  |
| typedef IdxRange&lt; [**GridR**](classIPolarFootFinder.md#typedef-gridr), [**GridTheta**](classIPolarFootFinder.md#typedef-gridtheta) &gt; | [**IdxRangeRTheta**](classIPolarFootFinder.md#typedef-idxrangertheta)  <br>_The type of the index range over which the operator works._  |
| typedef typename GridR::continuous\_dimension\_type | [**R**](classIPolarFootFinder.md#typedef-r)  <br>_The continuous radial dimension._  |
| typedef typename GridTheta::continuous\_dimension\_type | [**Theta**](classIPolarFootFinder.md#typedef-theta)  <br>_The continuous poloidal dimension._  |
| typedef AdvectionDim1 | [**X**](classIPolarFootFinder.md#typedef-x)  <br>_The continuous radial dimension._  |
| typedef AdvectionDim2 | [**Y**](classIPolarFootFinder.md#typedef-y)  <br>_The continuous poloidal dimension._  |
| typedef MemorySpace | [**memory\_space**](classIPolarFootFinder.md#typedef-memory_space)  <br>_The type of the memory space where the field is saved (CPU vs GPU)._  |














































## Detailed Description


The natural advection domain is the physical domain, where the studied equation is given. However, not all the mappings used are analytically invertible and inverting the Jacobian matrix of the mapping could be costly and could introduce numerical errors. That is why, we also introduce a pseudo-Cartesian domain.


More details can be found in Edoardo Zoni's article ([https://doi.org/10.1016/j.jcp.2019.108889](https://doi.org/10.1016/j.jcp.2019.108889)).




**Template parameters:**


* `TimeStepper` A child class of [**ITimeStepper**](classITimeStepper.md) providing a time integration method. 
* `LogicalToPhysicalMapping` A mapping from the logical domain to the physical domain. 
* `LogicalToPseudoPhysicalMapping` A mapping from the logical domain to the domain where the advection is carried out. This may be a pseudo-physical domain or the physical domain itself. 
* `SplineRThetaBuilder` A 2D SplineBuilder to construct a spline on a polar domain. 
* `SplineRThetaEvaluatorConstBound` A 2D SplineEvaluator to evaluate a spline on a polar domain. A boundary condition must be provided in case the foot of the characteristic is found outside the domain.



**See also:** [**BslAdvectionRTheta**](classBslAdvectionRTheta.md) 



    
## Public Types Documentation




### typedef ConstFieldRTheta 

_The type of a constant field on the polar plane on a compatible memory space._ 
```C++
using SplinePolarFootFinder< TimeStepper, LogicalToPhysicalMapping, LogicalToPseudoPhysicalMapping, SplineRThetaBuilder, SplineRThetaEvaluatorConstBound >::ConstFieldRTheta =  ConstField<ElementType, IdxRangeRTheta, memory_space>;
```




<hr>



### typedef DVectorConstFieldRTheta 

_The type of a constant vector (x,y) field on the polar plane on a compatible memory space._ 
```C++
using SplinePolarFootFinder< TimeStepper, LogicalToPhysicalMapping, LogicalToPseudoPhysicalMapping, SplineRThetaBuilder, SplineRThetaEvaluatorConstBound >::DVectorConstFieldRTheta =  VectorConstField<double, IdxRangeRTheta, VectorIndexSet<Dim1, Dim2>, memory_space>;
```




<hr>



### typedef DVectorFieldRTheta 

_The type of a vector (x,y) field on the polar plane on a compatible memory space._ 
```C++
using SplinePolarFootFinder< TimeStepper, LogicalToPhysicalMapping, LogicalToPseudoPhysicalMapping, SplineRThetaBuilder, SplineRThetaEvaluatorConstBound >::DVectorFieldRTheta =  VectorField<double, IdxRangeRTheta, VectorIndexSet<Dim1, Dim2>, memory_space>;
```




<hr>



### typedef FieldRTheta 

_The type of a field on the polar plane on a compatible memory space._ 
```C++
using SplinePolarFootFinder< TimeStepper, LogicalToPhysicalMapping, LogicalToPseudoPhysicalMapping, SplineRThetaBuilder, SplineRThetaEvaluatorConstBound >::FieldRTheta =  Field<ElementType, IdxRangeRTheta, memory_space>;
```




<hr>



### typedef VectorSplineCoeffsMem2D 

_The type of 2 splines representing the x and y components of a vector on the polar plane on a compatible memory space._ 
```C++
using SplinePolarFootFinder< TimeStepper, LogicalToPhysicalMapping, LogicalToPseudoPhysicalMapping, SplineRThetaBuilder, SplineRThetaEvaluatorConstBound >::VectorSplineCoeffsMem2D =  VectorFieldMem< double, IdxRange<BSplinesR, BSplinesTheta>, VectorIndexSet<Dim1, Dim2>, memory_space>;
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
    SplineRThetaBuilder const & builder_advection_field,
    SplineRThetaEvaluatorConstBound const & evaluator_advection_field,
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
    FieldRTheta < CoordRTheta > feet,
    DVectorConstFieldRTheta < X , Y > advection_field,
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

