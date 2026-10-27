

# Class RK3&lt; FieldMem, DerivFieldMem, ExecSpace &gt;

**template &lt;timestepper\_detail::FieldLike FieldMem, timestepper\_detail::FieldLike DerivFieldMem, class ExecSpace&gt;**



[**ClassList**](annotated.md) **>** [**RK3&lt; FieldMem, DerivFieldMem, ExecSpace &gt;**](classRK3_3_01FieldMem_00_01DerivFieldMem_00_01ExecSpace_01_4.md)



_A class which provides an implementation of a third-order Runge-Kutta method._ [More...](#detailed-description)

* `#include <rk3.hpp>`



Inherits the following classes: [ITimeStepper](classITimeStepper.md)














## Public Types

| Type | Name |
| ---: | :--- |
| typedef timestepper\_detail::const\_reference\_t&lt; [**DerivFieldMem**](classRK3_3_01FieldMem_00_01DerivFieldMem_00_01ExecSpace_01_4.md#typedef-derivfieldmem)[**DerivFieldMem**](classDerivFieldMem.md) &gt; | [**DerivConstField**](#typedef-derivconstfield)  <br>_The constant type of the derivatives values of the function being evolved._  |
| typedef timestepper\_detail::reference\_t&lt; [**DerivFieldMem**](classRK3_3_01FieldMem_00_01DerivFieldMem_00_01ExecSpace_01_4.md#typedef-derivfieldmem)[**DerivFieldMem**](classDerivFieldMem.md) &gt; | [**DerivField**](#typedef-derivfield)  <br>_The type of the derivatives of the function being evolved._  |
| typedef [**DerivFieldMem**](classRK3_3_01FieldMem_00_01DerivFieldMem_00_01ExecSpace_01_4.md#typedef-derivfieldmem)[**DerivFieldMem**](classDerivFieldMem.md) | [**DerivFieldMem**](#typedef-derivfieldmem)  <br>_The type of the memory allocation for the derivatives of the function being evolved._  |
| typedef typename timestepper\_detail::IdxRangeType&lt; FieldMem &gt;::type | [**IdxRange**](#typedef-idxrange)  <br>_The type of the index range on which the values of the function are defined._  |
| typedef timestepper\_detail::const\_reference\_t&lt; FieldMem &gt; | [**ValConstField**](#typedef-valconstfield)  <br>_The constant type of the values of the function being evolved._  |
| typedef timestepper\_detail::reference\_t&lt; FieldMem &gt; | [**ValField**](#typedef-valfield)  <br>_The type of the values of the function being evolved._  |
| typedef FieldMem | [**ValFieldMem**](#typedef-valfieldmem)  <br>_The type of the memory allocation for the values of the function being evolved._  |
| typedef ExecSpace | [**exec\_space**](#typedef-exec_space)  <br>_The space (CPU/GPU) where the calculations are carried out._  |


## Public Types inherited from ITimeStepper

See [ITimeStepper](classITimeStepper.md)

| Type | Name |
| ---: | :--- |
| typedef timestepper\_detail::const\_reference\_t&lt; [**DerivFieldMem**](classITimeStepper.md#typedef-derivfieldmem)[**DerivFieldMem**](classDerivFieldMem.md) &gt; | [**DerivConstField**](classITimeStepper.md#typedef-derivconstfield)  <br>_The constant type of the derivatives values of the function being evolved._  |
| typedef timestepper\_detail::reference\_t&lt; [**DerivFieldMem**](classITimeStepper.md#typedef-derivfieldmem)[**DerivFieldMem**](classDerivFieldMem.md) &gt; | [**DerivField**](classITimeStepper.md#typedef-derivfield)  <br>_The type of the derivatives of the function being evolved._  |
| typedef DerivFieldMemType | [**DerivFieldMem**](classITimeStepper.md#typedef-derivfieldmem)  <br>_The type of the memory allocation for the derivatives of the function being evolved._  |
| typedef typename timestepper\_detail::IdxRangeType&lt; FieldMem &gt;::type | [**IdxRange**](classITimeStepper.md#typedef-idxrange)  <br>_The type of the index range on which the values of the function are defined._  |
| typedef timestepper\_detail::const\_reference\_t&lt; FieldMem &gt; | [**ValConstField**](classITimeStepper.md#typedef-valconstfield)  <br>_The constant type of the values of the function being evolved._  |
| typedef timestepper\_detail::reference\_t&lt; FieldMem &gt; | [**ValField**](classITimeStepper.md#typedef-valfield)  <br>_The type of the values of the function being evolved._  |
| typedef FieldMem | [**ValFieldMem**](classITimeStepper.md#typedef-valfieldmem)  <br>_The type of the memory allocation for the values of the function being evolved._  |
| typedef ExecSpace | [**exec\_space**](classITimeStepper.md#typedef-exec_space)  <br>_The space (CPU/GPU) where the calculations are carried out._  |






































## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**RK3**](#function-rk3) ([**IdxRange**](classRK3_3_01FieldMem_00_01DerivFieldMem_00_01ExecSpace_01_4.md#typedef-idxrange) idx\_range) <br>_Create a_ [_**RK3**_](classRK3.md) _object._ |
| virtual void | [**update**](#function-update-13) (ExecSpace const & exec\_space, [**ValField**](classRK3_3_01FieldMem_00_01DerivFieldMem_00_01ExecSpace_01_4.md#typedef-valfield) y, double dt, std::function&lt; void([**DerivField**](classRK3_3_01FieldMem_00_01DerivFieldMem_00_01ExecSpace_01_4.md#typedef-derivfield)[**DerivField**](classDerivField.md), [**ValConstField**](classRK3_3_01FieldMem_00_01DerivFieldMem_00_01ExecSpace_01_4.md#typedef-valconstfield))&gt; dy\_calculator, std::function&lt; void([**ValField**](classRK3_3_01FieldMem_00_01DerivFieldMem_00_01ExecSpace_01_4.md#typedef-valfield), [**DerivConstField**](classRK3_3_01FieldMem_00_01DerivFieldMem_00_01ExecSpace_01_4.md#typedef-derivconstfield), double)&gt; y\_update) const<br>_Carry out one step of the Runge-Kutta scheme._  |
|  void | [**update**](#function-update-23) ([**ValField**](classRK3_3_01FieldMem_00_01DerivFieldMem_00_01ExecSpace_01_4.md#typedef-valfield) y, double dt, std::function&lt; void([**DerivField**](classRK3_3_01FieldMem_00_01DerivFieldMem_00_01ExecSpace_01_4.md#typedef-derivfield)[**DerivField**](classDerivField.md), [**ValConstField**](classRK3_3_01FieldMem_00_01DerivFieldMem_00_01ExecSpace_01_4.md#typedef-valconstfield))&gt; dy\_calculator) const<br>_Carry out one step of the timestepping scheme._  |
|  void | [**update**](#function-update-33) (ExecSpace const & exec\_space, [**ValField**](classRK3_3_01FieldMem_00_01DerivFieldMem_00_01ExecSpace_01_4.md#typedef-valfield) y, double dt, std::function&lt; void([**DerivField**](classRK3_3_01FieldMem_00_01DerivFieldMem_00_01ExecSpace_01_4.md#typedef-derivfield)[**DerivField**](classDerivField.md), [**ValConstField**](classRK3_3_01FieldMem_00_01DerivFieldMem_00_01ExecSpace_01_4.md#typedef-valconstfield))&gt; dy\_calculator) const<br>_Carry out one step of the timestepping scheme._  |


## Public Functions inherited from ITimeStepper

See [ITimeStepper](classITimeStepper.md)

| Type | Name |
| ---: | :--- |
|  void | [**update**](classITimeStepper.md#function-update-13) ([**ValField**](classITimeStepper.md#typedef-valfield) y, double dt, std::function&lt; void([**DerivField**](classITimeStepper.md#typedef-derivfield)[**DerivField**](classDerivField.md), [**ValConstField**](classITimeStepper.md#typedef-valconstfield))&gt; dy\_calculator) const<br>_Carry out one step of the timestepping scheme._  |
|  void | [**update**](classITimeStepper.md#function-update-23) (ExecSpace const & exec\_space, [**ValField**](classITimeStepper.md#typedef-valfield) y, double dt, std::function&lt; void([**DerivField**](classITimeStepper.md#typedef-derivfield)[**DerivField**](classDerivField.md), [**ValConstField**](classITimeStepper.md#typedef-valconstfield))&gt; dy\_calculator) const<br>_Carry out one step of the timestepping scheme._  |
| virtual void | [**update**](classITimeStepper.md#function-update-33) (ExecSpace const & exec\_space, [**ValField**](classITimeStepper.md#typedef-valfield) y, double dt, std::function&lt; void([**DerivField**](classITimeStepper.md#typedef-derivfield)[**DerivField**](classDerivField.md), [**ValConstField**](classITimeStepper.md#typedef-valconstfield))&gt; dy\_calculator, std::function&lt; void([**ValField**](classITimeStepper.md#typedef-valfield), [**DerivConstField**](classITimeStepper.md#typedef-derivconstfield), double)&gt; y\_update) const = 0<br>_Carry out one step of the timestepping scheme._  |






















































## Detailed Description


A class which provides an implementation of a third-order Runge-Kutta method in order to evolve values over time. The values may be either scalars or vectors. In the case of vectors the appropriate dimensions must be passed as template parameters. This specialisation handles Field-like objects (Field, [**VectorField**](classVectorField.md), MultipatchField). The values which evolve are defined on an index range.


For the following ODE : \(\partial_t y(t) = f(t, y(t))\),


the Runge-Kutta 3 method is given by : \(y^{n+1} =  y^{n} + \frac{dt}{6} \left(k_1 + 4 k_2 + k_3 \right)\),


with



* \(k_1 = f(t^{n}, y^{n})\),
* \(k_2 = f(t^{n+1/2}, y^{n} + \frac{dt}{2} k_1 )\),
* \(k_3 = f(t^{n+1/2}, y^{n} + dt ( 2 k_2 - k_1) )\). 




    
## Public Types Documentation




### typedef DerivConstField 

_The constant type of the derivatives values of the function being evolved._ 
```C++
using ITimeStepper< FieldMem, DerivFieldMem, ExecSpace >::DerivConstField = ;
```




<hr>



### typedef DerivField 

_The type of the derivatives of the function being evolved._ 
```C++
using ITimeStepper< FieldMem, DerivFieldMem, ExecSpace >::DerivField = ;
```




<hr>



### typedef DerivFieldMem 

_The type of the memory allocation for the derivatives of the function being evolved._ 
```C++
using RK3< FieldMem, DerivFieldMem, ExecSpace >::DerivFieldMem = ;
```




<hr>



### typedef IdxRange 

_The type of the index range on which the values of the function are defined._ 
```C++
using ITimeStepper< FieldMem, DerivFieldMem, ExecSpace >::IdxRange = ;
```




<hr>



### typedef ValConstField 

_The constant type of the values of the function being evolved._ 
```C++
using ITimeStepper< FieldMem, DerivFieldMem, ExecSpace >::ValConstField = ;
```




<hr>



### typedef ValField 

_The type of the values of the function being evolved._ 
```C++
using ITimeStepper< FieldMem, DerivFieldMem, ExecSpace >::ValField = ;
```




<hr>



### typedef ValFieldMem 

_The type of the memory allocation for the values of the function being evolved._ 
```C++
using RK3< FieldMem, DerivFieldMem, ExecSpace >::ValFieldMem = ;
```




<hr>



### typedef exec\_space 

_The space (CPU/GPU) where the calculations are carried out._ 
```C++
using RK3< FieldMem, DerivFieldMem, ExecSpace >::exec_space = ;
```




<hr>
## Public Functions Documentation




### function RK3 

_Create a_ [_**RK3**_](classRK3.md) _object._
```C++
inline explicit RK3< FieldMem, DerivFieldMem, ExecSpace >::RK3 (
    IdxRange idx_range
) 
```





**Parameters:**


* `idx_range` The index range on which the points which evolve over time are defined. 




        

<hr>



### function update [1/3]

_Carry out one step of the Runge-Kutta scheme._ 
```C++
inline virtual void RK3< FieldMem, DerivFieldMem, ExecSpace >::update (
    ExecSpace const & exec_space,
    ValField y,
    double dt,
    std::function< void( DerivField DerivField , ValConstField )> dy_calculator,
    std::function< void( ValField , DerivConstField , double)> y_update
) const
```





**Parameters:**


* [**exec\_space**](classRK3_3_01FieldMem_00_01DerivFieldMem_00_01ExecSpace_01_4.md#typedef-exec_space) The space on which the function is executed (CPU/GPU). 
* `y` The value(s) which should be evolved over time defined on each of the dimensions at each point of the index range. 
* `dt` The time step over which the values should be evolved. 
* `dy_calculator` The function describing how the derivative of the evolve function is calculated. 
* `y_update` The function describing how the value(s) are updated using the derivative. 




        
Implements [*ITimeStepper::update*](classITimeStepper.md#function-update-33)


<hr>



### function update [2/3]

_Carry out one step of the timestepping scheme._ 
```C++
inline void RK3< FieldMem, DerivFieldMem, ExecSpace >::update (
    ValField y,
    double dt,
    std::function< void( DerivField DerivField , ValConstField )> dy_calculator
) const
```



This function is a wrapper around the update function below. The values of the function are updated using the trivial method $f += df \* dt$. This is the standard method however some cases may need a more complex update function which is why the more explicit method is also provided.




**Parameters:**


* `y` The value(s) which should be evolved over time defined on each of the dimensions at each point of the index range. 
* `dt` The time step over which the values should be evolved. 
* `dy_calculator` The function describing how the derivative of the evolve function is calculated. 




        

<hr>



### function update [3/3]

_Carry out one step of the timestepping scheme._ 
```C++
inline void RK3< FieldMem, DerivFieldMem, ExecSpace >::update (
    ExecSpace const & exec_space,
    ValField y,
    double dt,
    std::function< void( DerivField DerivField , ValConstField )> dy_calculator
) const
```



This function is a wrapper around the update function below. The values of the function are updated using the trivial method $f += df \* dt$. This is the standard method however some cases may need a more complex update function which is why the more explicit method is also provided.




**Parameters:**


* [**exec\_space**](classRK3_3_01FieldMem_00_01DerivFieldMem_00_01ExecSpace_01_4.md#typedef-exec_space) The space on which the function is executed (CPU/GPU). 
* `y` The value(s) which should be evolved over time defined on each of the dimensions at each point of the index range. 
* `dt` The time step over which the values should be evolved. 
* `dy_calculator` The function describing how the derivative of the evolve function is calculated. 




        

<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/timestepper/rk3.hpp`

