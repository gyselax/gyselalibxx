

# Class RK4

**template &lt;class ValType, class DerivType, class ExecSpace&gt;**



[**ClassList**](annotated.md) **>** [**RK4**](classRK4.md)



_A class which provides an implementation of a fourth-order Runge-Kutta method._ [More...](#detailed-description)

* `#include <rk4.hpp>`

















## Public Types

| Type | Name |
| ---: | :--- |
| typedef DerivType | [**DerivFieldMem**](#typedef-derivfieldmem)  <br>_The type of the memory allocation for the derivatives of the function being evolved._  |
| typedef ValType | [**ValFieldMem**](#typedef-valfieldmem)  <br>_The type of the memory allocation for the values of the function being evolved._  |
| typedef ExecSpace | [**exec\_space**](#typedef-exec_space)  <br>_The space (CPU/GPU) where the calculations are carried out._  |




















## Public Functions

| Type | Name |
| ---: | :--- |
|  KOKKOS\_DEFAULTED\_FUNCTION | [**RK4**](#function-rk4) () = default<br>_Create a_ [_**RK4**_](classRK4.md) _object to operate on scalars._ |
|  KOKKOS\_FUNCTION void | [**update**](#function-update) (ValType & y, double dt, DYFunctor dy\_calculator, YFunctor y\_update=timestepper\_detail::serial\_y\_update&lt; ValType &, DerivType const & &gt;) const<br>_Carry out one step of the Runge-Kutta scheme on a scalar._  |




























## Detailed Description


A class which provides an implementation of a fourth-order Runge-Kutta method in order to evolve values over time. This specialisation handles elementwise operations and can be called from GPU.


For the following ODE : \(\partial_t y(t) = f(t, y(t))\),


the Runge-Kutta 3 method is given by : \(y^{n+1} =  y^{n} + \frac{dt}{6} \left(k_1 + 2k_2 + 2k_3 + k_4 \right)\),


with



* \(k_1 = f(t^{n}, y^{n})\),
* \(k_2 = f(t^{n+1/2}, y^{n} + \frac{dt}{2} k_1 )\),
* \(k_3 = f(t^{n+1/2}, y^{n} + \frac{dt}{2} k_2 )\),
* \(k_3 = f(t^{n}, y^{n} + dt k_3 )\). 




    
## Public Types Documentation




### typedef DerivFieldMem 

_The type of the memory allocation for the derivatives of the function being evolved._ 
```C++
using RK4< ValType, DerivType, ExecSpace >::DerivFieldMem =  DerivType;
```




<hr>



### typedef ValFieldMem 

_The type of the memory allocation for the values of the function being evolved._ 
```C++
using RK4< ValType, DerivType, ExecSpace >::ValFieldMem =  ValType;
```




<hr>



### typedef exec\_space 

_The space (CPU/GPU) where the calculations are carried out._ 
```C++
using RK4< ValType, DerivType, ExecSpace >::exec_space =  ExecSpace;
```




<hr>
## Public Functions Documentation




### function RK4 

_Create a_ [_**RK4**_](classRK4.md) _object to operate on scalars._
```C++
explicit KOKKOS_DEFAULTED_FUNCTION RK4::RK4 () = default
```




<hr>



### function update 

_Carry out one step of the Runge-Kutta scheme on a scalar._ 
```C++
template<class DYFunctor, class YFunctor>
inline KOKKOS_FUNCTION void RK4::update (
    ValType & y,
    double dt,
    DYFunctor dy_calculator,
    YFunctor y_update=timestepper_detail::serial_y_update< ValType &, DerivType const & >
) const
```





**Parameters:**


* `y` The value(s) which should be evolved over time defined on each of the dimensions at each point of the index range. 
* `dt` The time step over which the values should be evolved. 
* `dy_calculator` The function describing how the derivative of the evolve function is calculated. 
* `y_update` The function describing how the value(s) are updated using the derivative. 




        

<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/timestepper/rk4.hpp`

