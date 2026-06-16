

# Class PolarFootFinder

**template &lt;FootFindingSpace FFSpace, AdvectionFieldSpace AFSpace, concepts::Mapping LogicalToPhysicalMapping, class IdxRangeBatched, class TimeStepperBuilder, class RThetaAdvectionBuilder, class RThetaAdvectionEvaluator&gt;**



[**ClassList**](annotated.md) **>** [**PolarFootFinder**](classPolarFootFinder.md)



_Operator for finding the feet of the characteristics on a polar slice._ [More...](#detailed-description)

* `#include <polar_foot_finder.hpp>`

















## Public Types

| Type | Name |
| ---: | :--- |
| typedef ConstField&lt; CoordRTheta, IdxRangeBatched, [**memory\_space**](classPolarFootFinder.md#typedef-memory_space) &gt; | [**CConstFieldFeet**](#typedef-cconstfieldfeet)  <br>_A read-only field of_ \((r, \theta)\) _coordinates over the batched operator domain._ |
| typedef Field&lt; CoordRTheta, IdxRangeBatched, [**memory\_space**](classPolarFootFinder.md#typedef-memory_space) &gt; | [**CFieldFeet**](#typedef-cfieldfeet)  <br> |
| typedef typename [**polar\_foot\_finder\_details::ElementwiseChoice**](classpolar__foot__finder__details_1_1ElementwiseChoice.md)&lt; FFSpace, AFSpace, [**GridR**](structGridR.md), [**GridTheta**](structGridTheta.md), IdxRangeBatched, RThetaAdvectionEvaluator, [**AdvecCoefField**](classVectorFieldMem.md), TimeStepperBuilder, LogicalToPhysicalMapping &gt;::type | [**ElementwiseOperator**](#typedef-elementwiseoperator)  <br>_The operator returned by operator() which calculates the feet elementwise._  |
| typedef typename RThetaAdvectionBuilder::exec\_space | [**ExecSpace**](#typedef-execspace)  <br>_The execution space where kernels are launched._  |
| typedef IdxRangeBatched | [**IdxRangeOperator**](#typedef-idxrangeoperator)  <br>_The type of the index range over which the operator works._  |
| typedef typename LogicalToPhysicalMapping::curvilinear\_tag\_r | [**R**](#typedef-r)  <br>_The continuous radial dimension._  |
| typedef typename LogicalToPhysicalMapping::curvilinear\_tag\_theta | [**Theta**](#typedef-theta)  <br>_The continuous poloidal dimension._  |
| typedef typename RThetaAdvectionBuilder::memory\_space | [**memory\_space**](#typedef-memory_space)  <br>_The memory space where fields are stored (e.g._ `Kokkos::HostSpace` _or a GPU space)._ |




















## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**PolarFootFinder**](#function-polarfootfinder) (TimeStepperBuilder const & time\_stepper\_builder, LogicalToPhysicalMapping const & logical\_to\_physical, RThetaAdvectionBuilder const & builder\_advection\_field, RThetaAdvectionEvaluator const & evaluator\_advection\_field, Coord&lt; [**X\_pC**](structX__pC.md), [**Y\_pC**](structY__pC.md) &gt; coord\_centre\_pc=Coord&lt; [**X\_pC**](structX__pC.md), [**Y\_pC**](structY__pC.md) &gt;(0, 0), double epsilon=1e-12) <br>_Construct a_ [_**PolarFootFinder**_](classPolarFootFinder.md) _._ |
|  [**ElementwiseOperator**](classPolarFootFinder.md#typedef-elementwiseoperator) | [**operator()**](#function-operator) ([**DVectorConstField**](classVectorField.md)&lt; IdxRangeBatched, VectorIndexSet&lt; AdvDim1, AdvDim2 &gt;, [**memory\_space**](classPolarFootFinder.md#typedef-memory_space) &gt; advection\_field) const<br>_Get an elementwise operator capable of calculating the feet of the characteristics._  |
|  void | [**operator()**](#function-operator_1) ([**CFieldFeet**](classPolarFootFinder.md#typedef-cfieldfeet) feet, [**DVectorConstField**](classVectorField.md)&lt; IdxRangeBatched, VectorIndexSet&lt; AdvDim1, AdvDim2 &gt;, [**memory\_space**](classPolarFootFinder.md#typedef-memory_space) &gt; advection\_field, double dt) const<br>_Advect the feet over_ \(dt\) _._ |




























## Detailed Description


Calculates the spline representation of the advection field and uses it together with a time-stepping method to solve the characteristic equation. The space in which the advection field is expressed and the space in which foot-finding is performed are selected at compile time via `FFSpace` and `AFSpace`.




**Template parameters:**


* `FFSpace` The space in which the foot of the characteristic is found. 
* `AFSpace` The space in which the advection field is expressed. 
* `LogicalToPhysicalMapping` The mapping from the logical \((r, \theta)\) domain to the physical \((x, y)\) domain. 
* `IdxRangeBatched` The batched index range over which the operator works. 
* `TimeStepperBuilder` The factory type for the time integration method. 
* `RThetaAdvectionBuilder` The spline builder for the advection field. 
* `RThetaAdvectionEvaluator` The spline evaluator for the advection field. 




    
## Public Types Documentation




### typedef CConstFieldFeet 

_A read-only field of_ \((r, \theta)\) _coordinates over the batched operator domain._
```C++
using PolarFootFinder< FFSpace, AFSpace, LogicalToPhysicalMapping, IdxRangeBatched, TimeStepperBuilder, RThetaAdvectionBuilder, RThetaAdvectionEvaluator >::CConstFieldFeet =  ConstField<CoordRTheta, IdxRangeBatched, memory_space>;
```




<hr>



### typedef CFieldFeet 

```C++
using PolarFootFinder< FFSpace, AFSpace, LogicalToPhysicalMapping, IdxRangeBatched, TimeStepperBuilder, RThetaAdvectionBuilder, RThetaAdvectionEvaluator >::CFieldFeet =  Field<CoordRTheta, IdxRangeBatched, memory_space>;
```



The type of the memory space where the field is saved (CPU vs GPU). A mutable field of \((r, \theta)\) coordinates over the batched operator domain. 


        

<hr>



### typedef ElementwiseOperator 

_The operator returned by operator() which calculates the feet elementwise._ 
```C++
using PolarFootFinder< FFSpace, AFSpace, LogicalToPhysicalMapping, IdxRangeBatched, TimeStepperBuilder, RThetaAdvectionBuilder, RThetaAdvectionEvaluator >::ElementwiseOperator =  typename polar_foot_finder_details::ElementwiseChoice< FFSpace, AFSpace, GridR, GridTheta, IdxRangeBatched, RThetaAdvectionEvaluator, AdvecCoefField, TimeStepperBuilder, LogicalToPhysicalMapping>::type;
```




<hr>



### typedef ExecSpace 

_The execution space where kernels are launched._ 
```C++
using PolarFootFinder< FFSpace, AFSpace, LogicalToPhysicalMapping, IdxRangeBatched, TimeStepperBuilder, RThetaAdvectionBuilder, RThetaAdvectionEvaluator >::ExecSpace =  typename RThetaAdvectionBuilder::exec_space;
```




<hr>



### typedef IdxRangeOperator 

_The type of the index range over which the operator works._ 
```C++
using PolarFootFinder< FFSpace, AFSpace, LogicalToPhysicalMapping, IdxRangeBatched, TimeStepperBuilder, RThetaAdvectionBuilder, RThetaAdvectionEvaluator >::IdxRangeOperator =  IdxRangeBatched;
```




<hr>



### typedef R 

_The continuous radial dimension._ 
```C++
using PolarFootFinder< FFSpace, AFSpace, LogicalToPhysicalMapping, IdxRangeBatched, TimeStepperBuilder, RThetaAdvectionBuilder, RThetaAdvectionEvaluator >::R =  typename LogicalToPhysicalMapping::curvilinear_tag_r;
```




<hr>



### typedef Theta 

_The continuous poloidal dimension._ 
```C++
using PolarFootFinder< FFSpace, AFSpace, LogicalToPhysicalMapping, IdxRangeBatched, TimeStepperBuilder, RThetaAdvectionBuilder, RThetaAdvectionEvaluator >::Theta =  typename LogicalToPhysicalMapping::curvilinear_tag_theta;
```




<hr>



### typedef memory\_space 

_The memory space where fields are stored (e.g._ `Kokkos::HostSpace` _or a GPU space)._
```C++
using PolarFootFinder< FFSpace, AFSpace, LogicalToPhysicalMapping, IdxRangeBatched, TimeStepperBuilder, RThetaAdvectionBuilder, RThetaAdvectionEvaluator >::memory_space =  typename RThetaAdvectionBuilder::memory_space;
```




<hr>
## Public Functions Documentation




### function PolarFootFinder 

_Construct a_ [_**PolarFootFinder**_](classPolarFootFinder.md) _._
```C++
inline PolarFootFinder::PolarFootFinder (
    TimeStepperBuilder const & time_stepper_builder,
    LogicalToPhysicalMapping const & logical_to_physical,
    RThetaAdvectionBuilder const & builder_advection_field,
    RThetaAdvectionEvaluator const & evaluator_advection_field,
    Coord< X_pC , Y_pC > coord_centre_pc=Coord< X_pC , Y_pC >(0, 0),
    double epsilon=1e-12
) 
```





**Parameters:**


* `time_stepper_builder` A builder for the time integration method. 
* `logical_to_physical` The mapping from the logical domain to the physical domain. 
* `builder_advection_field` The spline builder for the advection field coefficients. 
* `evaluator_advection_field` The spline evaluator for the advection field. 
* `coord_centre_pc` The coordinate of the polar centre in the pseudo-Cartesian domain \((X_{pC}, Y_{pC})\). Ignored for `LOGICAL` foot finding. 
* `epsilon` \(\varepsilon\) parameter for the [**CombinedMapping**](classCombinedMapping.md) linearisation near the O-point. Only used for `PHYSICAL` advection field. 




        

<hr>



### function operator() 

_Get an elementwise operator capable of calculating the feet of the characteristics._ 
```C++
inline ElementwiseOperator PolarFootFinder::operator() (
    DVectorConstField < IdxRangeBatched, VectorIndexSet< AdvDim1, AdvDim2 >, memory_space > advection_field
) const
```



Computes the spline coefficients of the advection field, then packages them together with the mappings and time stepper into an [**ElementwiseOperator**](classPolarFootFinder.md#typedef-elementwiseoperator). Calling `operator()(dt)` on the returned object yields a GPU-copyable functor.




**Parameters:**


* `advection_field` The advection field in the domain selected by `AFSpace`. 



**Returns:**

An [**ElementwiseOperator**](classPolarFootFinder.md#typedef-elementwiseoperator) for the given advection field. 





        

<hr>



### function operator() 

_Advect the feet over_ \(dt\) _._
```C++
inline void PolarFootFinder::operator() (
    CFieldFeet feet,
    DVectorConstField < IdxRangeBatched, VectorIndexSet< AdvDim1, AdvDim2 >, memory_space > advection_field,
    double dt
) const
```



Computes the spline coefficients of the advection field, solves the characteristic equation over \(dt\) at every grid point, and writes the resulting feet in-place.




**Parameters:**


* `feet` On input: the mesh points. On output: the characteristic feet. 
* `advection_field` The advection field in the domain selected by `AFSpace`. 
* `dt` The time step. 




        

<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/advection/polar_foot_finder.hpp`

