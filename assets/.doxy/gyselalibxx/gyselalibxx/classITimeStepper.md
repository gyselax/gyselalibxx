

# Class ITimeStepper

**template &lt;class FieldMem, class DerivFieldMemType, class ExecSpace&gt;**



[**ClassList**](annotated.md) **>** [**ITimeStepper**](classITimeStepper.md)



_The superclass from which all timestepping methods inherit._ [More...](#detailed-description)

* `#include <itimestepper.hpp>`

















## Public Types

| Type | Name |
| ---: | :--- |
| typedef timestepper\_detail::const\_reference\_t&lt; [**DerivFieldMem**](classITimeStepper.md#typedef-derivfieldmem) &gt; | [**DerivConstField**](#typedef-derivconstfield)  <br>_The constant type of the derivatives values of the function being evolved._  |
| typedef timestepper\_detail::reference\_t&lt; [**DerivFieldMem**](classITimeStepper.md#typedef-derivfieldmem) &gt; | [**DerivField**](#typedef-derivfield)  <br>_The type of the derivatives of the function being evolved._  |
| typedef DerivFieldMemType | [**DerivFieldMem**](#typedef-derivfieldmem)  <br>_The type of the memory allocation for the derivatives of the function being evolved._  |
| typedef typename timestepper\_detail::IdxRangeType&lt; FieldMem &gt;::type | [**IdxRange**](#typedef-idxrange)  <br>_The type of the index range on which the values of the function are defined._  |
| typedef timestepper\_detail::const\_reference\_t&lt; FieldMem &gt; | [**ValConstField**](#typedef-valconstfield)  <br>_The constant type of the values of the function being evolved._  |
| typedef timestepper\_detail::reference\_t&lt; FieldMem &gt; | [**ValField**](#typedef-valfield)  <br>_The type of the values of the function being evolved._  |
| typedef FieldMem | [**ValFieldMem**](#typedef-valfieldmem)  <br>_The type of the memory allocation for the values of the function being evolved._  |
| typedef ExecSpace | [**exec\_space**](#typedef-exec_space)  <br>_The space (CPU/GPU) where the calculations are carried out._  |




















## Public Functions

| Type | Name |
| ---: | :--- |
|  void | [**update**](#function-update-13) ([**ValField**](classITimeStepper.md#typedef-valfield) y, double dt, std::function&lt; void([**DerivField**](classITimeStepper.md#typedef-derivfield), [**ValConstField**](classITimeStepper.md#typedef-valconstfield))&gt; dy\_calculator) const<br>_Carry out one step of the timestepping scheme._  |
|  void | [**update**](#function-update-23) (ExecSpace const & exec\_space, [**ValField**](classITimeStepper.md#typedef-valfield) y, double dt, std::function&lt; void([**DerivField**](classITimeStepper.md#typedef-derivfield), [**ValConstField**](classITimeStepper.md#typedef-valconstfield))&gt; dy\_calculator) const<br>_Carry out one step of the timestepping scheme._  |
| virtual void | [**update**](#function-update-33) (ExecSpace const & exec\_space, [**ValField**](classITimeStepper.md#typedef-valfield) y, double dt, std::function&lt; void([**DerivField**](classITimeStepper.md#typedef-derivfield), [**ValConstField**](classITimeStepper.md#typedef-valconstfield))&gt; dy\_calculator, std::function&lt; void([**ValField**](classITimeStepper.md#typedef-valfield), [**DerivConstField**](classITimeStepper.md#typedef-derivconstfield), double)&gt; y\_update) const = 0<br>_Carry out one step of the timestepping scheme._  |




























## Detailed Description


The class exposes three update functions which are used to carry out one step of the chosen timestepping method to solve an ODE of the form: \(\partial_t y(t) = f(t, y(t))\), 


    
## Public Types Documentation




### typedef DerivConstField 

_The constant type of the derivatives values of the function being evolved._ 
```C++
using ITimeStepper< FieldMem, DerivFieldMemType, ExecSpace >::DerivConstField =  timestepper_detail::const_reference_t<DerivFieldMem>;
```




<hr>



### typedef DerivField 

_The type of the derivatives of the function being evolved._ 
```C++
using ITimeStepper< FieldMem, DerivFieldMemType, ExecSpace >::DerivField =  timestepper_detail::reference_t<DerivFieldMem>;
```




<hr>



### typedef DerivFieldMem 

_The type of the memory allocation for the derivatives of the function being evolved._ 
```C++
using ITimeStepper< FieldMem, DerivFieldMemType, ExecSpace >::DerivFieldMem =  DerivFieldMemType;
```




<hr>



### typedef IdxRange 

_The type of the index range on which the values of the function are defined._ 
```C++
using ITimeStepper< FieldMem, DerivFieldMemType, ExecSpace >::IdxRange =  typename timestepper_detail::IdxRangeType<FieldMem>::type;
```




<hr>



### typedef ValConstField 

_The constant type of the values of the function being evolved._ 
```C++
using ITimeStepper< FieldMem, DerivFieldMemType, ExecSpace >::ValConstField =  timestepper_detail::const_reference_t<FieldMem>;
```




<hr>



### typedef ValField 

_The type of the values of the function being evolved._ 
```C++
using ITimeStepper< FieldMem, DerivFieldMemType, ExecSpace >::ValField =  timestepper_detail::reference_t<FieldMem>;
```




<hr>



### typedef ValFieldMem 

_The type of the memory allocation for the values of the function being evolved._ 
```C++
using ITimeStepper< FieldMem, DerivFieldMemType, ExecSpace >::ValFieldMem =  FieldMem;
```




<hr>



### typedef exec\_space 

_The space (CPU/GPU) where the calculations are carried out._ 
```C++
using ITimeStepper< FieldMem, DerivFieldMemType, ExecSpace >::exec_space =  ExecSpace;
```




<hr>
## Public Functions Documentation




### function update [1/3]

_Carry out one step of the timestepping scheme._ 
```C++
inline void ITimeStepper::update (
    ValField y,
    double dt,
    std::function< void( DerivField , ValConstField )> dy_calculator
) const
```



This function is a wrapper around the update function below. The values of the function are updated using the trivial method $f += df \* dt$. This is the standard method however some cases may need a more complex update function which is why the more explicit method is also provided.




**Parameters:**


* `y` The value(s) which should be evolved over time defined on each of the dimensions at each point of the index range. 
* `dt` The time step over which the values should be evolved. 
* `dy_calculator` The function describing how the derivative of the evolve function is calculated. 




        

<hr>



### function update [2/3]

_Carry out one step of the timestepping scheme._ 
```C++
inline void ITimeStepper::update (
    ExecSpace const & exec_space,
    ValField y,
    double dt,
    std::function< void( DerivField , ValConstField )> dy_calculator
) const
```



This function is a wrapper around the update function below. The values of the function are updated using the trivial method $f += df \* dt$. This is the standard method however some cases may need a more complex update function which is why the more explicit method is also provided.




**Parameters:**


* `exec_space` The space on which the function is executed (CPU/GPU). 
* `y` The value(s) which should be evolved over time defined on each of the dimensions at each point of the index range. 
* `dt` The time step over which the values should be evolved. 
* `dy_calculator` The function describing how the derivative of the evolve function is calculated. 




        

<hr>



### function update [3/3]

_Carry out one step of the timestepping scheme._ 
```C++
virtual void ITimeStepper::update (
    ExecSpace const & exec_space,
    ValField y,
    double dt,
    std::function< void( DerivField , ValConstField )> dy_calculator,
    std::function< void( ValField , DerivConstField , double)> y_update
) const = 0
```





**Parameters:**


* `exec_space` The space on which the function is executed (CPU/GPU). 
* `y` The value(s) which should be evolved over time defined on each of the dimensions at each point of the index range. 
* `dt` The time step over which the values should be evolved. 
* `dy_calculator` The function describing how the derivative of the evolve function is calculated. 
* `y_update` The function describing how the value(s) are updated using the derivative. 




        

<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/timestepper/itimestepper.hpp`

