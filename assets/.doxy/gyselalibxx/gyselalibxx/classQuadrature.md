

# Class Quadrature

**template &lt;class IdxRangeQuadrature, class IdxRangeTotal, class MemorySpace&gt;**



[**ClassList**](annotated.md) **>** [**Quadrature**](classQuadrature.md)



_A class providing an operator for integrating functions defined on known grid points._ [More...](#detailed-description)

* `#include <quadrature.hpp>`





































## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**Quadrature**](#function-quadrature) (QuadConstField coeffs) <br>_Create a_ [_**Quadrature**_](classQuadrature.md) _object._ |
|  double | [**operator()**](#function-operator) (ExecutionSpace exec\_space, IntegratorFunction integrated\_function) const<br>_An operator for calculating the integral of a function defined on known grid points._  |
|  void | [**operator()**](#function-operator_1) (ExecutionSpace exec\_space, Field&lt; double, BatchIdxRange, MemorySpace &gt; const result, IntegratorFunction integrated\_function) const<br>_An operator for calculating the integral of a function defined on a discrete domain by cycling over batch dimensions._  |




























## Detailed Description




**Template parameters:**


* `IdxRangeQuadrature` The index range over which the function is integrated. 
* `IdxRangeTotal` The index range of the chunk which can be passed to the operator(). This is the IdxRangeQuadrature combined with any batch dimensions. If there are no batch dimensions then this argument does not need to be provided as by default it is equal to the IdxRangeQuadrature. 
* `MemorySpace` The memory space (cpu/gpu) where the quadrature coefficients are saved. 




    
## Public Functions Documentation




### function Quadrature 

_Create a_ [_**Quadrature**_](classQuadrature.md) _object._
```C++
inline explicit Quadrature::Quadrature (
    QuadConstField coeffs
) 
```





**Parameters:**


* `coeffs` The coefficients of the quadrature. 




        

<hr>



### function operator() 

_An operator for calculating the integral of a function defined on known grid points._ 
```C++
template<class ExecutionSpace, class IntegratorFunction>
inline double Quadrature::operator() (
    ExecutionSpace exec_space,
    IntegratorFunction integrated_function
) const
```





**Parameters:**


* `exec_space` The space on which the function is executed (CPU/GPU). 
* `integrated_function` A function taking an index of a position in the index range over which the quadrature is calculated and returning the value of the function to be integrated at that point. It should be noted that a Field fulfils these criteria and can be passed as the function to be integrated. If the exec\_space is a GPU the function that is passed must be accessible from GPU.



**Returns:**

The integral of the function over the domain. 





        

<hr>



### function operator() 

_An operator for calculating the integral of a function defined on a discrete domain by cycling over batch dimensions._ 
```C++
template<class ExecutionSpace, class BatchIdxRange, class IntegratorFunction>
inline void Quadrature::operator() (
    ExecutionSpace exec_space,
    Field< double, BatchIdxRange, MemorySpace > const result,
    IntegratorFunction integrated_function
) const
```





**Parameters:**


* `exec_space` The space on which the function is executed (CPU/GPU). 
* `result` The result of the quadrature calculation. 
* `integrated_function` A function taking an index of a position in the index range over which the quadrature is calculated (including the batch index range) and returning the value of the function to be integrated at that point. Please note that a Field fulfils the described criteria. If the exec\_space is a GPU the function that is passed must be accessible from GPU. 




        

<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/quadrature/quadrature.hpp`

