

# Class Euler

**template &lt;class FieldMem, class [**DerivFieldMem**](classDerivFieldMem.md), class ExecSpace&gt;**



[**ClassList**](annotated.md) **>** [**Euler**](classEuler.md)



_A class which provides an implementation of an explicit_ [_**Euler**_](classEuler.md) _method._[More...](#detailed-description)

* `#include <euler.hpp>`



Inherits the following classes: [ITimeStepper](classITimeStepper.md)














## Public Types

| Type | Name |
| ---: | :--- |
| typedef typename DerivFieldMem::view\_type | [**DerivConstField**](#typedef-derivconstfield)  <br>_The constant type of the derivatives values of the function being evolved._  |
| typedef typename DerivFieldMem::span\_type | [**DerivField**](#typedef-derivfield)  <br>_The type of the derivatives of the function being evolved._  |
| typedef typename FieldMem::discrete\_domain\_type | [**IdxRange**](#typedef-idxrange)  <br>_The type of the index range on which the values of the function are defined._  |
| typedef typename FieldMem::view\_type | [**ValConstField**](#typedef-valconstfield)  <br>_The constant type of the values of the function being evolved._  |
| typedef typename FieldMem::span\_type | [**ValField**](#typedef-valfield)  <br>_The type of the values of the function being evolved._  |


## Public Types inherited from ITimeStepper

See [ITimeStepper](classITimeStepper.md)

| Type | Name |
| ---: | :--- |
| typedef typename DerivFieldMem::view\_type | [**DerivConstField**](classITimeStepper.md#typedef-derivconstfield)  <br>_The constant type of the derivatives values of the function being evolved._  |
| typedef typename DerivFieldMem::span\_type | [**DerivField**](classITimeStepper.md#typedef-derivfield)  <br>_The type of the derivatives of the function being evolved._  |
| typedef typename FieldMem::discrete\_domain\_type | [**IdxRange**](classITimeStepper.md#typedef-idxrange)  <br>_The type of the index range on which the values of the function are defined._  |
| typedef typename FieldMem::view\_type | [**ValConstField**](classITimeStepper.md#typedef-valconstfield)  <br>_The constant type of the values of the function being evolved._  |
| typedef typename FieldMem::span\_type | [**ValField**](classITimeStepper.md#typedef-valfield)  <br>_The type of the values of the function being evolved._  |






































## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**Euler**](#function-euler) ([**IdxRange**](classEuler.md#typedef-idxrange) idx\_range) <br>_Create a_ [_**Euler**_](classEuler.md) _object._ |
| virtual void | [**update**](#function-update-14) (ExecSpace const & exec\_space, [**ValField**](classEuler.md#typedef-valfield) y, double dt, std::function&lt; void([**DerivField**](classEuler.md#typedef-derivfield), [**ValConstField**](classEuler.md#typedef-valconstfield))&gt; dy\_calculator, std::function&lt; void([**ValField**](classEuler.md#typedef-valfield), [**DerivConstField**](classEuler.md#typedef-derivconstfield), double)&gt; y\_update) const<br>_Carry out one step of the explicit_ [_**Euler**_](classEuler.md) _scheme._ |
|  void | [**update**](#function-update-24) ([**ValField**](classEuler.md#typedef-valfield) y, double dt, std::function&lt; void([**DerivField**](classEuler.md#typedef-derivfield), [**ValConstField**](classEuler.md#typedef-valconstfield))&gt; dy\_calculator) const<br>_Carry out one step of the timestepping scheme._  |
|  void | [**update**](#function-update-34) (ExecSpace const & exec\_space, [**ValField**](classEuler.md#typedef-valfield) y, double dt, std::function&lt; void([**DerivField**](classEuler.md#typedef-derivfield), [**ValConstField**](classEuler.md#typedef-valconstfield))&gt; dy\_calculator) const<br>_Carry out one step of the timestepping scheme._  |
| virtual void | [**update**](#function-update-44) (ExecSpace const & exec\_space, [**ValField**](classEuler.md#typedef-valfield) y, double dt, std::function&lt; void([**DerivField**](classEuler.md#typedef-derivfield), [**ValConstField**](classEuler.md#typedef-valconstfield))&gt; dy\_calculator, std::function&lt; void([**ValField**](classEuler.md#typedef-valfield), [**DerivConstField**](classEuler.md#typedef-derivconstfield), double)&gt; y\_update) const<br>_Carry out one step of the timestepping scheme._  |


## Public Functions inherited from ITimeStepper

See [ITimeStepper](classITimeStepper.md)

| Type | Name |
| ---: | :--- |
|  void | [**assemble\_field\_k\_total**](classITimeStepper.md#function-assemble_field_k_total) (ExecSpace const & exec\_space, FieldType k\_total, FuncType func, std::array&lt; FieldType, n\_args &gt; k\_arr) const<br> |
|  void | [**assemble\_vector\_field\_k\_total**](classITimeStepper.md#function-assemble_vector_field_k_total) (ExecSpace const & exec\_space, FieldType k\_total, FuncType func, std::array&lt; FieldType, n\_args &gt; k\_arr) const<br> |
|  void | [**update**](classITimeStepper.md#function-update-13) ([**ValField**](classITimeStepper.md#typedef-valfield) y, double dt, std::function&lt; void([**DerivField**](classITimeStepper.md#typedef-derivfield), [**ValConstField**](classITimeStepper.md#typedef-valconstfield))&gt; dy\_calculator) const<br>_Carry out one step of the timestepping scheme._  |
|  void | [**update**](classITimeStepper.md#function-update-23) (ExecSpace const & exec\_space, [**ValField**](classITimeStepper.md#typedef-valfield) y, double dt, std::function&lt; void([**DerivField**](classITimeStepper.md#typedef-derivfield), [**ValConstField**](classITimeStepper.md#typedef-valconstfield))&gt; dy\_calculator) const<br>_Carry out one step of the timestepping scheme._  |
| virtual void | [**update**](classITimeStepper.md#function-update-33) (ExecSpace const & exec\_space, [**ValField**](classITimeStepper.md#typedef-valfield) y, double dt, std::function&lt; void([**DerivField**](classITimeStepper.md#typedef-derivfield), [**ValConstField**](classITimeStepper.md#typedef-valconstfield))&gt; dy\_calculator, std::function&lt; void([**ValField**](classITimeStepper.md#typedef-valfield), [**DerivConstField**](classITimeStepper.md#typedef-derivconstfield), double)&gt; y\_update) const = 0<br>_Carry out one step of the timestepping scheme._  |
















































## Protected Functions inherited from ITimeStepper

See [ITimeStepper](classITimeStepper.md)

| Type | Name |
| ---: | :--- |
|  void | [**assemble\_k\_total**](classITimeStepper.md#function-assemble_k_total) (ExecSpace const & exec\_space, [**DerivField**](classITimeStepper.md#typedef-derivfield) k\_total, FuncType func, T... k) const<br>_A method to assemble multiple derivative fields into one. This method is responsible for choosing how this is done depending on the type of the derivative field._  |
|  void | [**copy**](classITimeStepper.md#function-copy) ([**ValField**](classITimeStepper.md#typedef-valfield) copy\_to, [**ValConstField**](classITimeStepper.md#typedef-valconstfield) copy\_from) const<br>_Make a copy of the values of the function being evolved._  |




## Protected Static Functions inherited from ITimeStepper

See [ITimeStepper](classITimeStepper.md)

| Type | Name |
| ---: | :--- |
|  KOKKOS\_FUNCTION void | [**fill\_k\_total**](classITimeStepper.md#function-fill_k_total) (DerivFieldType k\_total, Idx i, [**DVector**](classTensor.md)&lt; DDims... &gt; new\_val) <br>_A method to fill an element of a vector field._  |


## Detailed Description


A class which provides an implementation of an explicit [**Euler**](classEuler.md) method in order to evolve values over time. The values may be either scalars or vectors. In the case of vectors the appropriate dimensions must be passed as template parameters. The values which evolve are defined on an index range.


For the following ODE : ,


the explicit [**Euler**](classEuler.md) method is given by : .


The method is order 1. 


    
## Public Types Documentation




### typedef DerivConstField 

_The constant type of the derivatives values of the function being evolved._ 
```C++
using ITimeStepper< FieldMem, DerivFieldMem, ExecSpace >::DerivConstField =  typename DerivFieldMem::view_type;
```




<hr>



### typedef DerivField 

_The type of the derivatives of the function being evolved._ 
```C++
using ITimeStepper< FieldMem, DerivFieldMem, ExecSpace >::DerivField =  typename DerivFieldMem::span_type;
```




<hr>



### typedef IdxRange 

_The type of the index range on which the values of the function are defined._ 
```C++
using ITimeStepper< FieldMem, DerivFieldMem, ExecSpace >::IdxRange =  typename FieldMem::discrete_domain_type;
```




<hr>



### typedef ValConstField 

_The constant type of the values of the function being evolved._ 
```C++
using ITimeStepper< FieldMem, DerivFieldMem, ExecSpace >::ValConstField =  typename FieldMem::view_type;
```




<hr>



### typedef ValField 

_The type of the values of the function being evolved._ 
```C++
using ITimeStepper< FieldMem, DerivFieldMem, ExecSpace >::ValField =  typename FieldMem::span_type;
```




<hr>
## Public Functions Documentation




### function Euler 

_Create a_ [_**Euler**_](classEuler.md) _object._
```C++
inline explicit Euler::Euler (
    IdxRange idx_range
) 
```





**Parameters:**


* `idx_range` The index range on which the points which evolve over time are defined. 




        

<hr>



### function update [1/4]

_Carry out one step of the explicit_ [_**Euler**_](classEuler.md) _scheme._
```C++
inline virtual void Euler::update (
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



### function update [2/4]

_Carry out one step of the timestepping scheme._ 
```C++
inline void Euler::update (
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



### function update [3/4]

_Carry out one step of the timestepping scheme._ 
```C++
inline void Euler::update (
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



### function update [4/4]

_Carry out one step of the timestepping scheme._ 
```C++
virtual void Euler::update (
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
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/timestepper/euler.hpp`

