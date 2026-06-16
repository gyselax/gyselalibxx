

# Class SplinePolarFootFinder

**template &lt;class IdxRangeBatched, class TimeStepperBuilder, concepts::Mapping LogicalToPhysicalMapping, concepts::AnalyticalMapping LogicalToPseudoPhysicalMapping, class SplineRThetaBuilderAdvection, class SplineRThetaEvaluatorAdvection&gt;**



[**ClassList**](annotated.md) **>** [**SplinePolarFootFinder**](classSplinePolarFootFinder.md)



_A class to find the foot of the characteristics on the_ \((r,\theta)\) _plane._[More...](#detailed-description)

* `#include <spline_polar_foot_finder.hpp>`

















## Public Types

| Type | Name |
| ---: | :--- |
| typedef ddc::type\_seq\_element\_t&lt; 0, [**VectorIndexSetAdvectionDims**](classSplinePolarFootFinder.md#typedef-vectorindexsetadvectiondims) &gt; | [**AdvDim1**](#typedef-advdim1)  <br>_The first dimension of the advection field._  |
| typedef ddc::type\_seq\_element\_t&lt; 1, [**VectorIndexSetAdvectionDims**](classSplinePolarFootFinder.md#typedef-vectorindexsetadvectiondims) &gt; | [**AdvDim2**](#typedef-advdim2)  <br>_The second dimension of the advection field._  |
| typedef ConstField&lt; CoordRTheta, [**IdxRangeOperator**](classSplinePolarFootFinder.md#typedef-idxrangeoperator), [**memory\_space**](classSplinePolarFootFinder.md#typedef-memory_space) &gt; | [**CConstFieldFeet**](#typedef-cconstfieldfeet)  <br>_The type of a constant field of (r, theta) coordinates at every grid point, saved on a compatible memory space._  |
| typedef Field&lt; CoordRTheta, [**IdxRangeOperator**](classSplinePolarFootFinder.md#typedef-idxrangeoperator), [**memory\_space**](classSplinePolarFootFinder.md#typedef-memory_space) &gt; | [**CFieldFeet**](#typedef-cfieldfeet)  <br>_The type of a field of (r, theta) coordinates at every grid point, saved on a compatible memory space._  |
| typedef [**DVectorConstField**](classVectorField.md)&lt; [**IdxRangeOperator**](classSplinePolarFootFinder.md#typedef-idxrangeoperator), PseudoCartesianBasis, [**memory\_space**](classSplinePolarFootFinder.md#typedef-memory_space) &gt; | [**DVectorConstFieldAdvection**](#typedef-dvectorconstfieldadvection)  <br>_The type of a constant vector field defined on the pseudo-Cartesian basis at every grid point, saved on a compatible memory space._  |
| typedef [**DVectorField**](classVectorField.md)&lt; [**IdxRangeOperator**](classSplinePolarFootFinder.md#typedef-idxrangeoperator), PseudoCartesianBasis, [**memory\_space**](classSplinePolarFootFinder.md#typedef-memory_space) &gt; | [**DVectorFieldAdvection**](#typedef-dvectorfieldadvection)  <br>_The type of a vector field defined on the pseudo-Cartesian basis at every grid point, saved on a compatible memory space._  |
| typedef [**ElementwiseSplinePolarFootFinderMem**](classElementwiseSplinePolarFootFinderMem.md)&lt; [**GridR**](classSplinePolarFootFinder.md#typedef-gridr), [**GridTheta**](classSplinePolarFootFinder.md#typedef-gridtheta), [**IdxRangeOperator**](classSplinePolarFootFinder.md#typedef-idxrangeoperator), SplineRThetaEvaluatorAdvection, [**PseudoPhysicalToAdvectionMapping**](classCombinedMapping.md), PseudoPhysicalToLogicalMapping, LogicalToPseudoPhysicalMapping, [**DVectorFieldMem**](classVectorFieldMem.md)&lt; IdxRangeSplineBatched, [**VectorIndexSetAdvectionDims**](classSplinePolarFootFinder.md#typedef-vectorindexsetadvectiondims), [**memory\_space**](classSplinePolarFootFinder.md#typedef-memory_space) &gt;, TimeStepper &gt; | [**ElementwiseOperator**](#typedef-elementwiseoperator)  <br>_The operator returned by operator() which calculates the feet elementwise._  |
| typedef typename SplineRThetaBuilderAdvection::exec\_space | [**ExecSpace**](#typedef-execspace)  <br>_Execution space._  |
| typedef typename SplineRThetaBuilderAdvection::interpolation\_discrete\_dimension\_type1 | [**GridR**](#typedef-gridr)  <br>_The continuous radial dimension._  |
| typedef typename SplineRThetaBuilderAdvection::interpolation\_discrete\_dimension\_type2 | [**GridTheta**](#typedef-gridtheta)  <br>_The continuous poloidal dimension._  |
| typedef IdxRangeBatched | [**IdxRangeOperator**](#typedef-idxrangeoperator)  <br>_The type of the index range over which the operator works._  |
| typedef typename GridR::continuous\_dimension\_type | [**R**](#typedef-r)  <br>_The continuous radial dimension._  |
| typedef typename GridTheta::continuous\_dimension\_type | [**Theta**](#typedef-theta)  <br>_The continuous poloidal dimension._  |
| typedef ddc::to\_type\_seq\_t&lt; typename LogicalToPhysicalMapping::CoordResult &gt; | [**VectorIndexSetAdvectionDims**](#typedef-vectorindexsetadvectiondims)  <br>_The continuous radial dimension._  |
| typedef [**DVectorFieldMem**](classVectorFieldMem.md)&lt; IdxRangeSplineBatched, PseudoCartesianBasis, [**memory\_space**](classSplinePolarFootFinder.md#typedef-memory_space) &gt; | [**VectorSplineCoeffsMem**](#typedef-vectorsplinecoeffsmem)  <br>_The type of 2 batched splines representing the x and y components of a vector on the polar plane on a compatible memory space._  |
| typedef typename SplineRThetaBuilderAdvection::memory\_space | [**memory\_space**](#typedef-memory_space)  <br>_The type of the memory space where the field is saved (CPU vs GPU)._  |




















## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**SplinePolarFootFinder**](#function-splinepolarfootfinder) (IdxRangeBatched const & idx\_range\_operator, TimeStepperBuilder const & time\_stepper\_builder, LogicalToPhysicalMapping const & logical\_to\_physical\_mapping, LogicalToPseudoPhysicalMapping const & logical\_to\_pseudo\_physical\_mapping, SplineRThetaBuilderAdvection const & builder\_advection\_field, SplineRThetaEvaluatorAdvection const & evaluator\_advection\_field, double epsilon=1e-12) <br>_Instantiate a time integration method for the advection operator._  |
|  [**ElementwiseOperator**](classSplinePolarFootFinder.md#typedef-elementwiseoperator) | [**operator()**](#function-operator) ([**DVectorConstField**](classVectorField.md)&lt; [**IdxRangeOperator**](classSplinePolarFootFinder.md#typedef-idxrangeoperator), [**VectorIndexSetAdvectionDims**](classSplinePolarFootFinder.md#typedef-vectorindexsetadvectiondims), [**memory\_space**](classSplinePolarFootFinder.md#typedef-memory_space) &gt; advection\_field) const<br>_Get an elementwise operator providing a GPU copyable functor capable of calculating the feet of the characteristics._  |
|  void | [**operator()**](#function-operator_1) ([**CFieldFeet**](classSplinePolarFootFinder.md#typedef-cfieldfeet) feet, [**DVectorConstField**](classVectorField.md)&lt; [**IdxRangeOperator**](classSplinePolarFootFinder.md#typedef-idxrangeoperator), [**VectorIndexSetAdvectionDims**](classSplinePolarFootFinder.md#typedef-vectorindexsetadvectiondims), [**memory\_space**](classSplinePolarFootFinder.md#typedef-memory_space) &gt; advection\_field, double dt) const<br>_Advect the feet over_ \(dt\) _._ |


## Public Static Functions

| Type | Name |
| ---: | :--- |
|  void | [**is\_unified**](#function-is_unified) (Field&lt; T, [**IdxRangeOperator**](classSplinePolarFootFinder.md#typedef-idxrangeoperator), [**memory\_space**](classSplinePolarFootFinder.md#typedef-memory_space) &gt; const & values) <br>_Check if the values at the centre point are the same._  |
|  void | [**unify\_value\_at\_centre\_pt**](#function-unify_value_at_centre_pt) (Field&lt; T, [**IdxRangeOperator**](classSplinePolarFootFinder.md#typedef-idxrangeoperator), [**memory\_space**](classSplinePolarFootFinder.md#typedef-memory_space) &gt; values) <br>_Replace the value at_ \((r=0, \theta)\) _point by the value at_\((r=0,0)\) _for all_\(\theta\) _._ |


























## Detailed Description


The natural advection domain is the physical domain, where the studied equation is given. However, not all the mappings used are analytically invertible and inverting the Jacobian matrix of the mapping could be costly and could introduce numerical errors. That is why, we also introduce a pseudo-Cartesian domain.


More details can be found in Edoardo Zoni's article ([https://doi.org/10.1016/j.jcp.2019.108889](https://doi.org/10.1016/j.jcp.2019.108889)).




**Template parameters:**


* `TimeStepperBuilder` A time stepper builder indicating which time integration method should be applied to solve the characteristic equation. 
* `LogicalToPhysicalMapping` A mapping from the logical domain to the physical domain. 
* `LogicalToPseudoPhysicalMapping` A mapping from the logical domain to the domain where the advection is carried out. This may be a pseudo-physical domain or the physical domain itself. 
* `SplineRThetaBuilderAdvection` A 2D SplineBuilder to construct a spline on a polar domain. 
* `SplineRThetaEvaluatorAdvection` A 2D SplineEvaluator to evaluate a spline on a polar domain. A boundary condition must be provided in case the foot of the characteristic is found outside the domain.



**See also:** [**BslAdvectionPolar**](classBslAdvectionPolar.md) 



    
## Public Types Documentation




### typedef AdvDim1 

_The first dimension of the advection field._ 
```C++
using SplinePolarFootFinder< IdxRangeBatched, TimeStepperBuilder, LogicalToPhysicalMapping, LogicalToPseudoPhysicalMapping, SplineRThetaBuilderAdvection, SplineRThetaEvaluatorAdvection >::AdvDim1 =  ddc::type_seq_element_t<0, VectorIndexSetAdvectionDims>;
```




<hr>



### typedef AdvDim2 

_The second dimension of the advection field._ 
```C++
using SplinePolarFootFinder< IdxRangeBatched, TimeStepperBuilder, LogicalToPhysicalMapping, LogicalToPseudoPhysicalMapping, SplineRThetaBuilderAdvection, SplineRThetaEvaluatorAdvection >::AdvDim2 =  ddc::type_seq_element_t<1, VectorIndexSetAdvectionDims>;
```




<hr>



### typedef CConstFieldFeet 

_The type of a constant field of (r, theta) coordinates at every grid point, saved on a compatible memory space._ 
```C++
using SplinePolarFootFinder< IdxRangeBatched, TimeStepperBuilder, LogicalToPhysicalMapping, LogicalToPseudoPhysicalMapping, SplineRThetaBuilderAdvection, SplineRThetaEvaluatorAdvection >::CConstFieldFeet =  ConstField<CoordRTheta, IdxRangeOperator, memory_space>;
```




<hr>



### typedef CFieldFeet 

_The type of a field of (r, theta) coordinates at every grid point, saved on a compatible memory space._ 
```C++
using SplinePolarFootFinder< IdxRangeBatched, TimeStepperBuilder, LogicalToPhysicalMapping, LogicalToPseudoPhysicalMapping, SplineRThetaBuilderAdvection, SplineRThetaEvaluatorAdvection >::CFieldFeet =  Field<CoordRTheta, IdxRangeOperator, memory_space>;
```




<hr>



### typedef DVectorConstFieldAdvection 

_The type of a constant vector field defined on the pseudo-Cartesian basis at every grid point, saved on a compatible memory space._ 
```C++
using SplinePolarFootFinder< IdxRangeBatched, TimeStepperBuilder, LogicalToPhysicalMapping, LogicalToPseudoPhysicalMapping, SplineRThetaBuilderAdvection, SplineRThetaEvaluatorAdvection >::DVectorConstFieldAdvection =  DVectorConstField<IdxRangeOperator, PseudoCartesianBasis, memory_space>;
```




<hr>



### typedef DVectorFieldAdvection 

_The type of a vector field defined on the pseudo-Cartesian basis at every grid point, saved on a compatible memory space._ 
```C++
using SplinePolarFootFinder< IdxRangeBatched, TimeStepperBuilder, LogicalToPhysicalMapping, LogicalToPseudoPhysicalMapping, SplineRThetaBuilderAdvection, SplineRThetaEvaluatorAdvection >::DVectorFieldAdvection =  DVectorField<IdxRangeOperator, PseudoCartesianBasis, memory_space>;
```




<hr>



### typedef ElementwiseOperator 

_The operator returned by operator() which calculates the feet elementwise._ 
```C++
using SplinePolarFootFinder< IdxRangeBatched, TimeStepperBuilder, LogicalToPhysicalMapping, LogicalToPseudoPhysicalMapping, SplineRThetaBuilderAdvection, SplineRThetaEvaluatorAdvection >::ElementwiseOperator =  ElementwiseSplinePolarFootFinderMem< GridR, GridTheta, IdxRangeOperator, SplineRThetaEvaluatorAdvection, PseudoPhysicalToAdvectionMapping, PseudoPhysicalToLogicalMapping, LogicalToPseudoPhysicalMapping, DVectorFieldMem<IdxRangeSplineBatched, VectorIndexSetAdvectionDims, memory_space>, TimeStepper>;
```




<hr>



### typedef ExecSpace 

_Execution space._ 
```C++
using SplinePolarFootFinder< IdxRangeBatched, TimeStepperBuilder, LogicalToPhysicalMapping, LogicalToPseudoPhysicalMapping, SplineRThetaBuilderAdvection, SplineRThetaEvaluatorAdvection >::ExecSpace =  typename SplineRThetaBuilderAdvection::exec_space;
```




<hr>



### typedef GridR 

_The continuous radial dimension._ 
```C++
using SplinePolarFootFinder< IdxRangeBatched, TimeStepperBuilder, LogicalToPhysicalMapping, LogicalToPseudoPhysicalMapping, SplineRThetaBuilderAdvection, SplineRThetaEvaluatorAdvection >::GridR =  typename SplineRThetaBuilderAdvection::interpolation_discrete_dimension_type1;
```




<hr>



### typedef GridTheta 

_The continuous poloidal dimension._ 
```C++
using SplinePolarFootFinder< IdxRangeBatched, TimeStepperBuilder, LogicalToPhysicalMapping, LogicalToPseudoPhysicalMapping, SplineRThetaBuilderAdvection, SplineRThetaEvaluatorAdvection >::GridTheta =  typename SplineRThetaBuilderAdvection::interpolation_discrete_dimension_type2;
```




<hr>



### typedef IdxRangeOperator 

_The type of the index range over which the operator works._ 
```C++
using SplinePolarFootFinder< IdxRangeBatched, TimeStepperBuilder, LogicalToPhysicalMapping, LogicalToPseudoPhysicalMapping, SplineRThetaBuilderAdvection, SplineRThetaEvaluatorAdvection >::IdxRangeOperator =  IdxRangeBatched;
```




<hr>



### typedef R 

_The continuous radial dimension._ 
```C++
using SplinePolarFootFinder< IdxRangeBatched, TimeStepperBuilder, LogicalToPhysicalMapping, LogicalToPseudoPhysicalMapping, SplineRThetaBuilderAdvection, SplineRThetaEvaluatorAdvection >::R =  typename GridR::continuous_dimension_type;
```




<hr>



### typedef Theta 

_The continuous poloidal dimension._ 
```C++
using SplinePolarFootFinder< IdxRangeBatched, TimeStepperBuilder, LogicalToPhysicalMapping, LogicalToPseudoPhysicalMapping, SplineRThetaBuilderAdvection, SplineRThetaEvaluatorAdvection >::Theta =  typename GridTheta::continuous_dimension_type;
```




<hr>



### typedef VectorIndexSetAdvectionDims 

_The continuous radial dimension._ 
```C++
using SplinePolarFootFinder< IdxRangeBatched, TimeStepperBuilder, LogicalToPhysicalMapping, LogicalToPseudoPhysicalMapping, SplineRThetaBuilderAdvection, SplineRThetaEvaluatorAdvection >::VectorIndexSetAdvectionDims =  ddc::to_type_seq_t<typename LogicalToPhysicalMapping::CoordResult>;
```




<hr>



### typedef VectorSplineCoeffsMem 

_The type of 2 batched splines representing the x and y components of a vector on the polar plane on a compatible memory space._ 
```C++
using SplinePolarFootFinder< IdxRangeBatched, TimeStepperBuilder, LogicalToPhysicalMapping, LogicalToPseudoPhysicalMapping, SplineRThetaBuilderAdvection, SplineRThetaEvaluatorAdvection >::VectorSplineCoeffsMem =  DVectorFieldMem<IdxRangeSplineBatched, PseudoCartesianBasis, memory_space>;
```




<hr>



### typedef memory\_space 

_The type of the memory space where the field is saved (CPU vs GPU)._ 
```C++
using SplinePolarFootFinder< IdxRangeBatched, TimeStepperBuilder, LogicalToPhysicalMapping, LogicalToPseudoPhysicalMapping, SplineRThetaBuilderAdvection, SplineRThetaEvaluatorAdvection >::memory_space =  typename SplineRThetaBuilderAdvection::memory_space;
```




<hr>
## Public Functions Documentation




### function SplinePolarFootFinder 

_Instantiate a time integration method for the advection operator._ 
```C++
inline SplinePolarFootFinder::SplinePolarFootFinder (
    IdxRangeBatched const & idx_range_operator,
    TimeStepperBuilder const & time_stepper_builder,
    LogicalToPhysicalMapping const & logical_to_physical_mapping,
    LogicalToPseudoPhysicalMapping const & logical_to_pseudo_physical_mapping,
    SplineRThetaBuilderAdvection const & builder_advection_field,
    SplineRThetaEvaluatorAdvection const & evaluator_advection_field,
    double epsilon=1e-12
) 
```





**Parameters:**


* `idx_range_operator` The index range on which the operator should act. 
* `time_stepper_builder` A builder for the time integration method used for the characteristic equation. 
* `logical_to_physical_mapping` The mapping from the logical domain to the physical domain. 
* `logical_to_pseudo_physical_mapping` The mapping from the logical domain to the pseudo-physical domain. 
* `builder_advection_field` The spline builder which computes the spline representation of the advection field. 
* `evaluator_advection_field` The B-splines evaluator to evaluate the advection field. 
* `epsilon` \(\varepsilon\) parameter used for the linearisation of the advection field around the central point.



**See also:** [**ITimeStepper**](classITimeStepper.md) 



        

<hr>



### function operator() 

_Get an elementwise operator providing a GPU copyable functor capable of calculating the feet of the characteristics._ 
```C++
inline ElementwiseOperator SplinePolarFootFinder::operator() (
    DVectorConstField < IdxRangeOperator , VectorIndexSetAdvectionDims , memory_space > advection_field
) const
```



From the advection field in the physical domain, compute the advection field in the right domain an compute its B-splines coefficients. Then, use the given time integration method (time\_stepper) to solve the characteristic equation over \(dt\).




**Parameters:**


* `advection_field` The advection field in the chosen domain.



**Returns:**

An elementwise operator providing a GPU copyable functor capable of calculating the feet of the characteristics. 





        

<hr>



### function operator() 

_Advect the feet over_ \(dt\) _._
```C++
inline void SplinePolarFootFinder::operator() (
    CFieldFeet feet,
    DVectorConstField < IdxRangeOperator , VectorIndexSetAdvectionDims , memory_space > advection_field,
    double dt
) const
```



From the advection field in the physical domain, compute the advection field in the right domain an compute its B-splines coefficients. Then, use the given time integration method (time\_stepper) to solve the characteristic equation over \(dt\).




**Parameters:**


* `feet` On input: the mesh points. On output: the characteristic feet. 
* `advection_field` The advection field in the physical domain. 
* `dt` The time step. 




        

<hr>
## Public Static Functions Documentation




### function is\_unified 

_Check if the values at the centre point are the same._ 
```C++
template<class T>
static inline void SplinePolarFootFinder::is_unified (
    Field< T, IdxRangeOperator , memory_space > const & values
) 
```



For polar geometry, to ensure continuity at the centre point, we have to be sure that all the points for \(r = 0\) have the same value. This function check if for \(r= 0\), the values \(\forall \theta\) are the same.




**Parameters:**


* `values` A table of values we want to check if the centre point has an unique value. 




        

<hr>



### function unify\_value\_at\_centre\_pt 

_Replace the value at_ \((r=0, \theta)\) _point by the value at_\((r=0,0)\) _for all_\(\theta\) _._
```C++
template<class T>
static inline void SplinePolarFootFinder::unify_value_at_centre_pt (
    Field< T, IdxRangeOperator , memory_space > values
) 
```



For polar geometry, to ensure continuity at the centre point, we have to be sure that all the points for \(r = 0\) have the same value. As the computation of the values of a table can induces machine errors, this function is useful to reset the values at the central point at the same value.




**Parameters:**


* `values` The table of values we want to unify at the central point. 




        

<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/advection/spline_polar_foot_finder.hpp`

