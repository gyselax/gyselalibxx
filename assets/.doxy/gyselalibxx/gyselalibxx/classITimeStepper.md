

# Class ITimeStepper

**template &lt;class FieldMem, class [**DerivFieldMem**](classDerivFieldMem.md), class ExecSpace&gt;**



[**ClassList**](annotated.md) **>** [**ITimeStepper**](classITimeStepper.md)



_The superclass from which all timestepping methods inherit._ [More...](#detailed-description)

* `#include <itimestepper.hpp>`





Inherited by the following classes: [CrankNicolson](classCrankNicolson.md),  [Euler](classEuler.md),  [RK2](classRK2.md),  [RK3](classRK3.md),  [RK4](classRK4.md)












## Public Types

| Type | Name |
| ---: | :--- |
| typedef typename DerivFieldMem::view\_type | [**DerivConstField**](#typedef-derivconstfield)  <br>_The constant type of the derivatives values of the function being evolved._  |
| typedef typename DerivFieldMem::span\_type | [**DerivField**](#typedef-derivfield)  <br>_The type of the derivatives of the function being evolved._  |
| typedef typename FieldMem::discrete\_domain\_type | [**IdxRange**](#typedef-idxrange)  <br>_The type of the index range on which the values of the function are defined._  |
| typedef typename FieldMem::view\_type | [**ValConstField**](#typedef-valconstfield)  <br>_The constant type of the values of the function being evolved._  |
| typedef typename FieldMem::span\_type | [**ValField**](#typedef-valfield)  <br>_The type of the values of the function being evolved._  |




















## Public Functions

| Type | Name |
| ---: | :--- |
|  void | [**assemble\_field\_k\_total**](#function-assemble_field_k_total) (ExecSpace const & exec\_space, FieldType k\_total, FuncType func, std::array&lt; FieldType, n\_args &gt; k\_arr) const<br> |
|  void | [**assemble\_vector\_field\_k\_total**](#function-assemble_vector_field_k_total) (ExecSpace const & exec\_space, FieldType k\_total, FuncType func, std::array&lt; FieldType, n\_args &gt; k\_arr) const<br> |
|  void | [**update**](#function-update-13) ([**ValField**](classITimeStepper.md#typedef-valfield) y, double dt, std::function&lt; void([**DerivField**](classITimeStepper.md#typedef-derivfield), [**ValConstField**](classITimeStepper.md#typedef-valconstfield))&gt; dy\_calculator) const<br>_Carry out one step of the timestepping scheme._  |
|  void | [**update**](#function-update-23) (ExecSpace const & exec\_space, [**ValField**](classITimeStepper.md#typedef-valfield) y, double dt, std::function&lt; void([**DerivField**](classITimeStepper.md#typedef-derivfield), [**ValConstField**](classITimeStepper.md#typedef-valconstfield))&gt; dy\_calculator) const<br>_Carry out one step of the timestepping scheme._  |
| virtual void | [**update**](#function-update-33) (ExecSpace const & exec\_space, [**ValField**](classITimeStepper.md#typedef-valfield) y, double dt, std::function&lt; void([**DerivField**](classITimeStepper.md#typedef-derivfield), [**ValConstField**](classITimeStepper.md#typedef-valconstfield))&gt; dy\_calculator, std::function&lt; void([**ValField**](classITimeStepper.md#typedef-valfield), [**DerivConstField**](classITimeStepper.md#typedef-derivconstfield), double)&gt; y\_update) const = 0<br>_Carry out one step of the timestepping scheme._  |
























## Protected Functions

| Type | Name |
| ---: | :--- |
|  void | [**assemble\_k\_total**](#function-assemble_k_total) (ExecSpace const & exec\_space, [**DerivField**](classITimeStepper.md#typedef-derivfield) k\_total, FuncType func, T... k) const<br>_A method to assemble multiple derivative fields into one. This method is responsible for choosing how this is done depending on the type of the derivative field._  |
|  void | [**copy**](#function-copy) ([**ValField**](classITimeStepper.md#typedef-valfield) copy\_to, [**ValConstField**](classITimeStepper.md#typedef-valconstfield) copy\_from) const<br>_Make a copy of the values of the function being evolved._  |


## Protected Static Functions

| Type | Name |
| ---: | :--- |
|  KOKKOS\_FUNCTION void | [**fill\_k\_total**](#function-fill_k_total) (DerivFieldType k\_total, Idx i, [**DVector**](classTensor.md)&lt; DDims... &gt; new\_val) <br>_A method to fill an element of a vector field._  |


## Detailed Description


The class exposes three update functions which are used to carry out one step of the chosen timestepping method to solve an ODE of the form: , 


    
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




### function assemble\_field\_k\_total 

```C++
template<class FieldType, class FuncType, std::size_t n_args>
inline void ITimeStepper::assemble_field_k_total (
    ExecSpace const & exec_space,
    FieldType k_total,
    FuncType func,
    std::array< FieldType, n_args > k_arr
) const
```



Calculate func(k\_arr[0], k\_arr[1], ...) when FieldType is a Field (ddc::ChunkSpan). This function should be private but is public due to Cuda restrictions.




**Parameters:**


* `exec_space` The space (CPU/GPU) where the calculation should be executed. 
* `k_total` The field to be filled with the combined derivative fields. 
* `func` A function which combines an element from each of the derivative fields. 
* `k_arr` The derivative fields being combined. 




        

<hr>



### function assemble\_vector\_field\_k\_total 

```C++
template<class FieldType, class FuncType, std::size_t n_args>
inline void ITimeStepper::assemble_vector_field_k_total (
    ExecSpace const & exec_space,
    FieldType k_total,
    FuncType func,
    std::array< FieldType, n_args > k_arr
) const
```



Calculate func(k\_arr[0], k\_arr[1], ...) when FieldType is a [**VectorField**](classVectorField.md). This function should be private but is public due to Cuda restrictions.




**Parameters:**


* `exec_space` The space (CPU/GPU) where the calculation should be executed. 
* `k_total` The field to be filled with the combined derivative fields. 
* `func` A function which combines an element from each of the derivative fields. 
* `k_arr` The derivative fields being combined. 




        

<hr>



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
## Protected Functions Documentation




### function assemble\_k\_total 

_A method to assemble multiple derivative fields into one. This method is responsible for choosing how this is done depending on the type of the derivative field._ 
```C++
template<class FuncType, class... T>
inline void ITimeStepper::assemble_k_total (
    ExecSpace const & exec_space,
    DerivField k_total,
    FuncType func,
    T... k
) const
```





**Parameters:**


* `exec_space` The space (CPU/GPU) where the calculation should be executed. 
* `k_total` The field to be filled with the combined derivative fields. 
* `func` A function which combines an element from each of the derivative fields. 
* `k` The derivative fields being combined. 




        

<hr>



### function copy 

_Make a copy of the values of the function being evolved._ 
```C++
inline void ITimeStepper::copy (
    ValField copy_to,
    ValConstField copy_from
) const
```





**Parameters:**


* `copy_to` the field that the values should be copied to. 
* `copy_from` The field that the values should be copied from. 




        

<hr>
## Protected Static Functions Documentation




### function fill\_k\_total 

_A method to fill an element of a vector field._ 
```C++
template<class DerivFieldType, class Idx, class... DDims>
static inline KOKKOS_FUNCTION void ITimeStepper::fill_k_total (
    DerivFieldType k_total,
    Idx i,
    DVector < DDims... > new_val
) 
```





**Parameters:**


* `k_total` The vector field that will be filled. 
* `i` The index where the vector field should be filled. 
* `new_val` The coordinate that should be saved to the vector field. 




        

<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/timestepper/itimestepper.hpp`

