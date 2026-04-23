

# Class RK2&lt; FieldMem, DerivFieldMem, ExecSpace &gt;

**template &lt;timestepper\_detail::FieldLike FieldMem, timestepper\_detail::FieldLike DerivFieldMem, class ExecSpace&gt;**



[**ClassList**](annotated.md) **>** [**RK2&lt; FieldMem, DerivFieldMem, ExecSpace &gt;**](classRK2_3_01FieldMem_00_01DerivFieldMem_00_01ExecSpace_01_4.md)



_A class which provides an implementation of a second-order Runge-Kutta method._ [More...](#detailed-description)

* `#include <rk2.hpp>`



Inherits the following classes: [ITimeStepper](classITimeStepper.md)
















## Public Types inherited from ITimeStepper

See [ITimeStepper](classITimeStepper.md)

| Type | Name |
| ---: | :--- |
| typedef timestepper\_detail::const\_reference\_t&lt; [**DerivFieldMem**](classITimeStepper.md#typedef-derivfieldmem) &gt; | [**DerivConstField**](classITimeStepper.md#typedef-derivconstfield)  <br>_The constant type of the derivatives values of the function being evolved._  |
| typedef timestepper\_detail::reference\_t&lt; [**DerivFieldMem**](classITimeStepper.md#typedef-derivfieldmem) &gt; | [**DerivField**](classITimeStepper.md#typedef-derivfield)  <br>_The type of the derivatives of the function being evolved._  |
| typedef DerivFieldMemType | [**DerivFieldMem**](classITimeStepper.md#typedef-derivfieldmem)  <br>_The type of the memory allocation for the derivatives of the function being evolved._  |
| typedef typename timestepper\_detail::IdxRangeType&lt; FieldMem &gt;::type | [**IdxRange**](classITimeStepper.md#typedef-idxrange)  <br>_The type of the index range on which the values of the function are defined._  |
| typedef timestepper\_detail::const\_reference\_t&lt; FieldMem &gt; | [**ValConstField**](classITimeStepper.md#typedef-valconstfield)  <br>_The constant type of the values of the function being evolved._  |
| typedef timestepper\_detail::reference\_t&lt; FieldMem &gt; | [**ValField**](classITimeStepper.md#typedef-valfield)  <br>_The type of the values of the function being evolved._  |
| typedef FieldMem | [**ValFieldMem**](classITimeStepper.md#typedef-valfieldmem)  <br>_The type of the memory allocation for the values of the function being evolved._  |
| typedef ExecSpace | [**exec\_space**](classITimeStepper.md#typedef-exec_space)  <br>_The space (CPU/GPU) where the calculations are carried out._  |






































## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**RK2**](#function-rk2) ([**IdxRange**](classITimeStepper.md#typedef-idxrange) idx\_range) <br>_Create a_ [_**RK2**_](classRK2.md) _object._ |
| virtual void | [**update**](#function-update) (ExecSpace const & exec\_space, [**ValField**](classITimeStepper.md#typedef-valfield) y, double dt, std::function&lt; void([**DerivField**](classITimeStepper.md#typedef-derivfield), [**ValConstField**](classITimeStepper.md#typedef-valconstfield))&gt; dy\_calculator, std::function&lt; void([**ValField**](classITimeStepper.md#typedef-valfield), [**DerivConstField**](classITimeStepper.md#typedef-derivconstfield), double)&gt; y\_update) const<br>_Carry out one step of the Runge-Kutta scheme._  |


## Public Functions inherited from ITimeStepper

See [ITimeStepper](classITimeStepper.md)

| Type | Name |
| ---: | :--- |
|  void | [**update**](classITimeStepper.md#function-update-13) ([**ValField**](classITimeStepper.md#typedef-valfield) y, double dt, std::function&lt; void([**DerivField**](classITimeStepper.md#typedef-derivfield), [**ValConstField**](classITimeStepper.md#typedef-valconstfield))&gt; dy\_calculator) const<br>_Carry out one step of the timestepping scheme._  |
|  void | [**update**](classITimeStepper.md#function-update-23) (ExecSpace const & exec\_space, [**ValField**](classITimeStepper.md#typedef-valfield) y, double dt, std::function&lt; void([**DerivField**](classITimeStepper.md#typedef-derivfield), [**ValConstField**](classITimeStepper.md#typedef-valconstfield))&gt; dy\_calculator) const<br>_Carry out one step of the timestepping scheme._  |
| virtual void | [**update**](classITimeStepper.md#function-update-33) (ExecSpace const & exec\_space, [**ValField**](classITimeStepper.md#typedef-valfield) y, double dt, std::function&lt; void([**DerivField**](classITimeStepper.md#typedef-derivfield), [**ValConstField**](classITimeStepper.md#typedef-valconstfield))&gt; dy\_calculator, std::function&lt; void([**ValField**](classITimeStepper.md#typedef-valfield), [**DerivConstField**](classITimeStepper.md#typedef-derivconstfield), double)&gt; y\_update) const = 0<br>_Carry out one step of the timestepping scheme._  |






















































## Detailed Description


A class which provides an implementation of a second-order Runge-Kutta method in order to evolve values over time. The values may be either scalars or vectors. In the case of vectors the appropriate dimensions must be passed as template parameters. This specialisation handles Field-like objects (Field, [**VectorField**](classVectorField.md), [**MultipatchField**](classMultipatchField.md)). The values which evolve are defined on an index range.


For the following ODE : \(\partial_t y(t) = f(t, y(t))\),


the Runge-Kutta 2 method is given by : \(y^{n+1} =  y^{n} + dt k_2\),


with



* \(k_1 = f(t^{n}, y^{n})\),
* \(k_2 = f(t^{n+1/2}, y^{n} + \frac{dt}{2} k_1 )\), 




    
## Public Functions Documentation




### function RK2 

_Create a_ [_**RK2**_](classRK2.md) _object._
```C++
inline explicit RK2< FieldMem, DerivFieldMem, ExecSpace >::RK2 (
    IdxRange idx_range
) 
```





**Parameters:**


* `idx_range` The index range on which the points which evolve over time are defined. 




        

<hr>



### function update 

_Carry out one step of the Runge-Kutta scheme._ 
```C++
inline virtual void RK2< FieldMem, DerivFieldMem, ExecSpace >::update (
    ExecSpace const & exec_space,
    ValField y,
    double dt,
    std::function< void( DerivField , ValConstField )> dy_calculator,
    std::function< void( ValField , DerivConstField , double)> y_update
) const
```





**Parameters:**


* `exec_space` The space on which the function is executed (CPU/GPU). 
* `y` The value(s) which should be evolved over time defined on each of the dimensions at each point of the index range. 
* `dt` The time step over which the values should be evolved. 
* `dy_calculator` The function describing how the derivative of the evolve function is calculated. 
* `y_update` The function describing how the value(s) are updated using the derivative. 




        
Implements [*ITimeStepper::update*](classITimeStepper.md#function-update-33)


<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/timestepper/rk2.hpp`

