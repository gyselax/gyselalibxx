

# Class CrankNicolson

**template &lt;class ValType, class DerivType, class ExecSpace&gt;**



[**ClassList**](annotated.md) **>** [**CrankNicolson**](classCrankNicolson.md)



_A class which provides an implementation of a Crank-Nicolson method._ [More...](#detailed-description)

* `#include <crank_nicolson.hpp>`

















## Public Types

| Type | Name |
| ---: | :--- |
| typedef DerivType | [**DerivFieldMem**](#typedef-derivfieldmem)  <br>_The type of the memory allocation for the derivatives of the function being evolved._  |
| typedef ValType | [**ValFieldMem**](#typedef-valfieldmem)  <br>_The type of the memory allocation for the values of the function being evolved._  |
| typedef ExecSpace | [**exec\_space**](#typedef-exec_space)  <br>_The space (CPU/GPU) where the calculations are carried out._  |




















## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**CrankNicolson**](#function-cranknicolson) (int const counter=int(20), double const epsilon=1e-12) <br>_Create a_ [_**CrankNicolson**_](classCrankNicolson.md) _object to operate on scalars._ |
|  KOKKOS\_FUNCTION void | [**update**](#function-update) (ValType & y, double dt, DYFunctor dy\_calculator, YFunctor y\_update=timestepper\_detail::serial\_y\_update&lt; ValType &, DerivType const & &gt;) const<br>_Carry out one step of the Crank-Nicolson scheme on a scalar._  |




























## Detailed Description


A class which provides an implementation of a Crank-Nicolson method in order to evolve values over time. This specialisation handles elementwise operations and can be called from GPU.


For the following ODE : \(\partial_t y(t) = f(t, y(t))\),


the Crank-Nicolson method is given by : \(y^{k} =  y^{n} + \frac{dt}{2} \left(f(t^{n}, y^{n}) + f(t^{k}, y^{k}) \right)\).


The method is an implicit method. If \(|y^{k+1} -  y^{k}| < \varepsilon\), then we set \(y^{n+1} = y^{k+1}\).


The method is order 2. 


    
## Public Types Documentation




### typedef DerivFieldMem 

_The type of the memory allocation for the derivatives of the function being evolved._ 
```C++
using CrankNicolson< ValType, DerivType, ExecSpace >::DerivFieldMem =  DerivType;
```




<hr>



### typedef ValFieldMem 

_The type of the memory allocation for the values of the function being evolved._ 
```C++
using CrankNicolson< ValType, DerivType, ExecSpace >::ValFieldMem =  ValType;
```




<hr>



### typedef exec\_space 

_The space (CPU/GPU) where the calculations are carried out._ 
```C++
using CrankNicolson< ValType, DerivType, ExecSpace >::exec_space =  ExecSpace;
```




<hr>
## Public Functions Documentation




### function CrankNicolson 

_Create a_ [_**CrankNicolson**_](classCrankNicolson.md) _object to operate on scalars._
```C++
inline explicit CrankNicolson::CrankNicolson (
    int const counter=int(20),
    double const epsilon=1e-12
) 
```





**Parameters:**


* `counter` The maximal number of loops for the implicit method. 
* `epsilon` The \(\varepsilon\) upperbound of the difference of two steps in the implicit method: \(|y^{k+1} -  y^{k}| < \varepsilon\). 




        

<hr>



### function update 

_Carry out one step of the Crank-Nicolson scheme on a scalar._ 
```C++
template<class DYFunctor, class YFunctor>
inline KOKKOS_FUNCTION void CrankNicolson::update (
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
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/timestepper/crank_nicolson.hpp`

