

# Class PolarFootFinder

**template &lt;FootFindingSpace FFSpace, AdvectionFieldSpace AFSpace, concepts::Mapping LogicalToPhysicalMapping, class IdxRangeBatched, class TimeStepperBuilder, concepts::Interpolation RThetaAdvectionInterpolator&gt;**



[**ClassList**](annotated.md) **>** [**PolarFootFinder**](classPolarFootFinder.md)



_Operator for finding the feet of the characteristics on a polar slice._ [More...](#detailed-description)

* `#include <polar_foot_finder.hpp>`

















## Public Types

| Type | Name |
| ---: | :--- |
| typedef std::conditional\_t&lt; AFSpace==AdvectionFieldSpace::PHYSICAL, [**X**](structX.md), [**R**](classPolarFootFinder.md#typedef-r) &gt; | [**AdvDim1**](#typedef-advdim1)  <br>_The first dimension on which the advection field is defined._  |
| typedef std::conditional\_t&lt; AFSpace==AdvectionFieldSpace::PHYSICAL, [**Y**](structY.md), [**Theta**](classPolarFootFinder.md#typedef-theta) &gt; | [**AdvDim2**](#typedef-advdim2)  <br>_The second dimension on which the advection field is defined._  |
| typedef ConstField&lt; CoordRTheta, IdxRangeBatched, [**memory\_space**](classPolarFootFinder.md#typedef-memory_space) &gt; | [**CConstFieldFeet**](#typedef-cconstfieldfeet)  <br>_A read-only field of_ \((r, \theta)\) _coordinates over the batched operator domain._ |
| typedef Field&lt; CoordRTheta, IdxRangeBatched, [**memory\_space**](classPolarFootFinder.md#typedef-memory_space) &gt; | [**CFieldFeet**](#typedef-cfieldfeet)  <br> |
| typedef typename [**polar\_foot\_finder\_details::ElementwiseChoice**](structpolar__foot__finder__details_1_1ElementwiseChoice.md)&lt; FFSpace, AFSpace, [**GridR**](structGridR.md), [**GridTheta**](structGridTheta.md), IdxRangeBatched, RThetaAdvectionEvaluator, [**AdvecCoefField**](classVectorFieldMem.md), TimeStepperBuilder, LogicalToPhysicalMapping &gt;::type | [**ElementwiseOperator**](#typedef-elementwiseoperator)  <br>_The operator returned by operator() which calculates the feet elementwise._  |
| typedef typename RThetaAdvectionBuilder::exec\_space | [**ExecSpace**](#typedef-execspace)  <br>_The execution space where kernels are launched._  |
| typedef IdxRangeBatched | [**IdxRangeOperator**](#typedef-idxrangeoperator)  <br>_The type of the index range over which the operator works._  |
| typedef typename [**CoordWithOPoint**](structCoordWithOPoint.md)&lt; typename LogicalToPhysicalMapping::CoordArg &gt;::curvilinear\_tag\_r | [**R**](#typedef-r)  <br>_The continuous radial dimension._  |
| typedef typename [**CoordWithOPoint**](structCoordWithOPoint.md)&lt; typename LogicalToPhysicalMapping::CoordArg &gt;::curvilinear\_tag\_theta | [**Theta**](#typedef-theta)  <br>_The continuous poloidal dimension._  |
| typedef typename RThetaAdvectionBuilder::memory\_space | [**memory\_space**](#typedef-memory_space)  <br>_The memory space where fields are stored (e.g._ `Kokkos::HostSpace` _or a GPU space)._ |




















## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**PolarFootFinder**](#function-polarfootfinder) (TimeStepperBuilder const & time\_stepper\_builder, LogicalToPhysicalMapping const & logical\_to\_physical, RThetaAdvectionInterpolator const & interpolator\_advection\_field, Coord&lt; [**X\_pC**](structX__pC.md), [**Y\_pC**](structY__pC.md) &gt; coord\_centre\_pc=Coord&lt; [**X\_pC**](structX__pC.md), [**Y\_pC**](structY__pC.md) &gt;(0, 0), double epsilon=1e-12) <br>_Construct a_ [_**PolarFootFinder**_](classPolarFootFinder.md) _._ |
|  [**ElementwiseOperator**](classPolarFootFinder.md#typedef-elementwiseoperator) | [**operator()**](#function-operator) ([**DVectorConstField**](classVectorField.md)&lt; IdxRangeBatched, VectorIndexSet&lt; [**AdvDim1**](classPolarFootFinder.md#typedef-advdim1), [**AdvDim2**](classPolarFootFinder.md#typedef-advdim2) &gt;, [**memory\_space**](classPolarFootFinder.md#typedef-memory_space) &gt; advection\_field) const<br>_Get an elementwise operator capable of calculating the feet of the characteristics._  |
|  void | [**operator()**](#function-operator_1) ([**CFieldFeet**](classPolarFootFinder.md#typedef-cfieldfeet) feet, [**DVectorConstField**](classVectorField.md)&lt; IdxRangeBatched, VectorIndexSet&lt; [**AdvDim1**](classPolarFootFinder.md#typedef-advdim1), [**AdvDim2**](classPolarFootFinder.md#typedef-advdim2) &gt;, [**memory\_space**](classPolarFootFinder.md#typedef-memory_space) &gt; advection\_field, double dt) const<br>_Advect the feet over_ \(dt\) _._ |




























## Detailed Description


Calculates the interpolation representation of the advection field and uses it together with a time-stepping method to solve the characteristic equation. The space in which the advection field is expressed and the space in which foot-finding is performed are selected at compile time via `FFSpace` and `AFSpace`.




**Template parameters:**


* `FFSpace` The space in which the foot of the characteristic is found. 
* `AFSpace` The space in which the advection field is expressed. 
* `LogicalToPhysicalMapping` The mapping from the logical \((r, \theta)\) domain to the physical \((x, y)\) domain. 
* `IdxRangeBatched` The batched index range over which the operator works. 
* `TimeStepperBuilder` The factory type for the time integration method. 
* `RThetaAdvectionInterpolator` The interpolator for the advection field. 




    
## Public Types Documentation




### typedef AdvDim1 

_The first dimension on which the advection field is defined._ 
```C++
using PolarFootFinder< FFSpace, AFSpace, LogicalToPhysicalMapping, IdxRangeBatched, TimeStepperBuilder, RThetaAdvectionInterpolator >::AdvDim1 =  std::conditional_t<AFSpace == AdvectionFieldSpace::PHYSICAL, X, R>;
```




<hr>



### typedef AdvDim2 

_The second dimension on which the advection field is defined._ 
```C++
using PolarFootFinder< FFSpace, AFSpace, LogicalToPhysicalMapping, IdxRangeBatched, TimeStepperBuilder, RThetaAdvectionInterpolator >::AdvDim2 =  std::conditional_t<AFSpace == AdvectionFieldSpace::PHYSICAL, Y, Theta>;
```




<hr>



### typedef CConstFieldFeet 

_A read-only field of_ \((r, \theta)\) _coordinates over the batched operator domain._
```C++
using PolarFootFinder< FFSpace, AFSpace, LogicalToPhysicalMapping, IdxRangeBatched, TimeStepperBuilder, RThetaAdvectionInterpolator >::CConstFieldFeet =  ConstField<CoordRTheta, IdxRangeBatched, memory_space>;
```




<hr>



### typedef CFieldFeet 

```C++
using PolarFootFinder< FFSpace, AFSpace, LogicalToPhysicalMapping, IdxRangeBatched, TimeStepperBuilder, RThetaAdvectionInterpolator >::CFieldFeet =  Field<CoordRTheta, IdxRangeBatched, memory_space>;
```



The type of the memory space where the field is saved (CPU vs GPU). A mutable field of \((r, \theta)\) coordinates over the batched operator domain. 


        

<hr>



### typedef ElementwiseOperator 

_The operator returned by operator() which calculates the feet elementwise._ 
```C++
using PolarFootFinder< FFSpace, AFSpace, LogicalToPhysicalMapping, IdxRangeBatched, TimeStepperBuilder, RThetaAdvectionInterpolator >::ElementwiseOperator =  typename polar_foot_finder_details::ElementwiseChoice< FFSpace, AFSpace, GridR, GridTheta, IdxRangeBatched, RThetaAdvectionEvaluator, AdvecCoefField, TimeStepperBuilder, LogicalToPhysicalMapping>::type;
```




<hr>



### typedef ExecSpace 

_The execution space where kernels are launched._ 
```C++
using PolarFootFinder< FFSpace, AFSpace, LogicalToPhysicalMapping, IdxRangeBatched, TimeStepperBuilder, RThetaAdvectionInterpolator >::ExecSpace =  typename RThetaAdvectionBuilder::exec_space;
```




<hr>



### typedef IdxRangeOperator 

_The type of the index range over which the operator works._ 
```C++
using PolarFootFinder< FFSpace, AFSpace, LogicalToPhysicalMapping, IdxRangeBatched, TimeStepperBuilder, RThetaAdvectionInterpolator >::IdxRangeOperator =  IdxRangeBatched;
```




<hr>



### typedef R 

_The continuous radial dimension._ 
```C++
using PolarFootFinder< FFSpace, AFSpace, LogicalToPhysicalMapping, IdxRangeBatched, TimeStepperBuilder, RThetaAdvectionInterpolator >::R =  typename CoordWithOPoint< typename LogicalToPhysicalMapping::CoordArg>::curvilinear_tag_r;
```




<hr>



### typedef Theta 

_The continuous poloidal dimension._ 
```C++
using PolarFootFinder< FFSpace, AFSpace, LogicalToPhysicalMapping, IdxRangeBatched, TimeStepperBuilder, RThetaAdvectionInterpolator >::Theta =  typename CoordWithOPoint< typename LogicalToPhysicalMapping::CoordArg>::curvilinear_tag_theta;
```




<hr>



### typedef memory\_space 

_The memory space where fields are stored (e.g._ `Kokkos::HostSpace` _or a GPU space)._
```C++
using PolarFootFinder< FFSpace, AFSpace, LogicalToPhysicalMapping, IdxRangeBatched, TimeStepperBuilder, RThetaAdvectionInterpolator >::memory_space =  typename RThetaAdvectionBuilder::memory_space;
```




<hr>
## Public Functions Documentation




### function PolarFootFinder 

_Construct a_ [_**PolarFootFinder**_](classPolarFootFinder.md) _._
```C++
inline PolarFootFinder::PolarFootFinder (
    TimeStepperBuilder const & time_stepper_builder,
    LogicalToPhysicalMapping const & logical_to_physical,
    RThetaAdvectionInterpolator const & interpolator_advection_field,
    Coord< X_pC , Y_pC > coord_centre_pc=Coord< X_pC , Y_pC >(0, 0),
    double epsilon=1e-12
) 
```





**Parameters:**


* `time_stepper_builder` A builder for the time integration method. 
* `logical_to_physical` The mapping from the logical domain to the physical domain. 
* `interpolator_advection_field` An interpolator to build and evaluate an approximation of the advection field. 
* `coord_centre_pc` The coordinate of the polar centre in the pseudo-Cartesian domain \((X_{pC}, Y_{pC})\). Ignored for `LOGICAL` foot finding. 
* `epsilon` \(\varepsilon\) parameter for the [**CombinedMapping**](classCombinedMapping.md) linearisation near the O-point. Only used for `PHYSICAL` advection field. 




        

<hr>



### function operator() 

_Get an elementwise operator capable of calculating the feet of the characteristics._ 
```C++
inline ElementwiseOperator PolarFootFinder::operator() (
    DVectorConstField < IdxRangeBatched, VectorIndexSet< AdvDim1 , AdvDim2 >, memory_space > advection_field
) const
```



Computes the interpolation coefficients of the advection field, then packages them together with the mappings and time stepper into an [**ElementwiseOperator**](classPolarFootFinder.md#typedef-elementwiseoperator). Calling `operator()(dt)` on the returned object yields a GPU-copyable functor.




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
    DVectorConstField < IdxRangeBatched, VectorIndexSet< AdvDim1 , AdvDim2 >, memory_space > advection_field,
    double dt
) const
```



Computes the coefficients of the advection field, solves the characteristic equation over \(dt\) at every grid point, and writes the resulting feet in-place.




**Parameters:**


* `feet` On input: the mesh points. On output: the characteristic feet. 
* `advection_field` The advection field in the domain selected by `AFSpace`. 
* `dt` The time step. 




        

<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/advection/polar_foot_finder.hpp`

