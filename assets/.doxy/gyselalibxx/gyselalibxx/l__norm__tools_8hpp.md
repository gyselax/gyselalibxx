

# File l\_norm\_tools.hpp



[**FileList**](files.md) **>** [**math\_tools**](dir_3ced5d1c6eac490d7704c2e023d148d8.md) **>** [**l\_norm\_tools.hpp**](l__norm__tools_8hpp.md)

[Go to the source code of this file](l__norm__tools_8hpp_source.md)

[More...](#detailed-description)

* `#include <ddc/ddc.hpp>`
* `#include "quadrature.hpp"`
* `#include "vector_field.hpp"`
* `#include "vector_index_tools.hpp"`





































## Public Functions

| Type | Name |
| ---: | :--- |
|  double | [**error\_norm\_L1**](#function-error_norm_l1) (ExecSpace exec\_space, [**Quadrature**](classQuadrature.md)&lt; IdxRangeQuad, IdxRangeQuad, typename ExecSpace::memory\_space &gt; quadrature, DField&lt; IdxRangeQuad, typename ExecSpace::memory\_space &gt; function, DField&lt; IdxRangeQuad, typename ExecSpace::memory\_space &gt; exact\_function) <br>_Compute the L1 norm of the error between 2 Fields._  |
|  double | [**error\_norm\_L2**](#function-error_norm_l2) (ExecSpace exec\_space, [**Quadrature**](classQuadrature.md)&lt; IdxRangeQuad, IdxRangeQuad, typename ExecSpace::memory\_space &gt; quadrature, DField&lt; IdxRangeQuad, typename ExecSpace::memory\_space &gt; function, DField&lt; IdxRangeQuad, typename ExecSpace::memory\_space &gt; exact\_function) <br>_Compute the L2 norm of the error between 2 Fields._  |
|  double | [**error\_norm\_inf**](#function-error_norm_inf) (ExecSpace exec\_space, ConstField&lt; ElementType, IdxRange, typename ExecSpace::memory\_space &gt; function, ConstField&lt; ElementType, IdxRange, typename ExecSpace::memory\_space &gt; exact\_function) <br>_Compute the infinity norm of the error between 2 Fields._  |
|  double | [**error\_norm\_inf**](#function-error_norm_inf) (ExecSpace exec\_space, [**VectorConstField**](classVectorField.md)&lt; ElementType, IdxRange, VectorIndexSetType, typename ExecSpace::memory\_space &gt; function, [**VectorConstField**](classVectorField.md)&lt; ElementType, IdxRange, VectorIndexSetType, typename ExecSpace::memory\_space &gt; exact\_function) <br>_Compute the infinity norm of the error between 2 VectorFields._  |
|  double | [**norm\_L1**](#function-norm_l1) (ExecSpace exec\_space, [**Quadrature**](classQuadrature.md)&lt; IdxRangeQuad, IdxRangeQuad, typename ExecSpace::memory\_space &gt; quadrature, DField&lt; IdxRangeQuad, typename ExecSpace::memory\_space &gt; function) <br>_Compute L1 norm of a function with a given quadrature._  |
|  double | [**norm\_L2**](#function-norm_l2) (ExecSpace exec\_space, [**Quadrature**](classQuadrature.md)&lt; IdxRangeQuad, IdxRangeQuad, typename ExecSpace::memory\_space &gt; quadrature, DField&lt; IdxRangeQuad, typename ExecSpace::memory\_space &gt; function) <br>_Compute L2 norm of a function with a given quadrature._  |
|  KOKKOS\_FUNCTION double | [**norm\_inf**](#function-norm_inf) (ddc::Coordinate&lt; Tags... &gt; coord) <br>_Compute the infinity norm._  |
|  KOKKOS\_FUNCTION double | [**norm\_inf**](#function-norm_inf) ([**DVector**](classTensor.md)&lt; Tags... &gt; vec) <br>_Compute the infinity norm of a vector on an orthonormal coordinate system._  |
|  KOKKOS\_INLINE\_FUNCTION double | [**norm\_inf**](#function-norm_inf) (double const coord) <br>_Compute the infinity norm._  |
|  double | [**norm\_inf**](#function-norm_inf) (ExecSpace exec\_space, ConstField&lt; ElementType, IdxRange, typename ExecSpace::memory\_space &gt; function) <br>_Compute the infinity norm for a Field._  |
|  double | [**norm\_inf**](#function-norm_inf) (ExecSpace exec\_space, [**VectorConstField**](classVectorField.md)&lt; ElementType, IdxRange, VectorIndexSetType, typename ExecSpace::memory\_space &gt; function) <br>_Compute the infinity norm for a_ [_**VectorField**_](classVectorField.md) _._ |




























## Detailed Description


File Describing useful mathematical functions to compute Lnorms 


    
## Public Functions Documentation




### function error\_norm\_L1 

_Compute the L1 norm of the error between 2 Fields._ 
```C++
template<class IdxRangeQuad, class ExecSpace>
double error_norm_L1 (
    ExecSpace exec_space,
    Quadrature < IdxRangeQuad, IdxRangeQuad, typename ExecSpace::memory_space > quadrature,
    DField< IdxRangeQuad, typename ExecSpace::memory_space > function,
    DField< IdxRangeQuad, typename ExecSpace::memory_space > exact_function
) 
```





**Parameters:**


* `exec_space` The space on which the function is executed (CPU/GPU). 
* `quadrature` The quadrature used to compute the integral. 
* `function` The calculated function. 
* `exact_function` The exact function with which the calculated function is compared. 



**Returns:**

A double containing the value of the infinity norm. 





        

<hr>



### function error\_norm\_L2 

_Compute the L2 norm of the error between 2 Fields._ 
```C++
template<class IdxRangeQuad, class ExecSpace>
double error_norm_L2 (
    ExecSpace exec_space,
    Quadrature < IdxRangeQuad, IdxRangeQuad, typename ExecSpace::memory_space > quadrature,
    DField< IdxRangeQuad, typename ExecSpace::memory_space > function,
    DField< IdxRangeQuad, typename ExecSpace::memory_space > exact_function
) 
```





**Parameters:**


* `exec_space` The space on which the function is executed (CPU/GPU). 
* `quadrature` The quadrature used to compute the integral. 
* `function` The calculated function. 
* `exact_function` The exact function with which the calculated function is compared. 



**Returns:**

A double containing the value of the infinity norm. 





        

<hr>



### function error\_norm\_inf 

_Compute the infinity norm of the error between 2 Fields._ 
```C++
template<class ExecSpace, class ElementType, class IdxRange>
inline double error_norm_inf (
    ExecSpace exec_space,
    ConstField< ElementType, IdxRange, typename ExecSpace::memory_space > function,
    ConstField< ElementType, IdxRange, typename ExecSpace::memory_space > exact_function
) 
```





**Parameters:**


* `exec_space` The space on which the function is executed (CPU/GPU). 
* `function` The calculated function. 
* `exact_function` The exact function with which the calculated function is compared. 



**Returns:**

A double containing the value of the infinity norm. 





        

<hr>



### function error\_norm\_inf 

_Compute the infinity norm of the error between 2 VectorFields._ 
```C++
template<class ExecSpace, class ElementType, class IdxRange, class VectorIndexSetType>
inline double error_norm_inf (
    ExecSpace exec_space,
    VectorConstField < ElementType, IdxRange, VectorIndexSetType, typename ExecSpace::memory_space > function,
    VectorConstField < ElementType, IdxRange, VectorIndexSetType, typename ExecSpace::memory_space > exact_function
) 
```





**Parameters:**


* `exec_space` The space on which the function is executed (CPU/GPU). 
* `function` The calculated function. 
* `exact_function` The exact function with which the calculated function is compared. 



**Returns:**

A double containing the value of the infinity norm. 





        

<hr>



### function norm\_L1 

_Compute L1 norm of a function with a given quadrature._ 
```C++
template<class IdxRangeQuad, class ExecSpace>
double norm_L1 (
    ExecSpace exec_space,
    Quadrature < IdxRangeQuad, IdxRangeQuad, typename ExecSpace::memory_space > quadrature,
    DField< IdxRangeQuad, typename ExecSpace::memory_space > function
) 
```








**Parameters:**


* `exec_space` The space on which the function is executed (CPU/GPU). 
* `quadrature` The quadrature used to compute the integral. 
* `function` A Field to the value of the function on the quadrature grid.



**Returns:**

A double containing the L1 norm of the function. 





        

<hr>



### function norm\_L2 

_Compute L2 norm of a function with a given quadrature._ 
```C++
template<class IdxRangeQuad, class ExecSpace>
double norm_L2 (
    ExecSpace exec_space,
    Quadrature < IdxRangeQuad, IdxRangeQuad, typename ExecSpace::memory_space > quadrature,
    DField< IdxRangeQuad, typename ExecSpace::memory_space > function
) 
```








**Parameters:**


* `exec_space` The space on which the function is executed (CPU/GPU). 
* `quadrature` The quadrature used to compute the integral. 
* `function` A Field to the value of the function on the quadrature grid.



**Returns:**

A double containing the L2 norm of the function. 





        

<hr>



### function norm\_inf 

_Compute the infinity norm._ 
```C++
template<class... Tags>
KOKKOS_FUNCTION double norm_inf (
    ddc::Coordinate< Tags... > coord
) 
```



For a given vector  , compute .




**Parameters:**


* `coord` The given vector.



**Returns:**

A double containing the value of the infinity norm. 





        

<hr>



### function norm\_inf 

_Compute the infinity norm of a vector on an orthonormal coordinate system._ 
```C++
template<class... Tags>
KOKKOS_FUNCTION double norm_inf (
    DVector < Tags... > vec
) 
```



For a given vector  , compute .




**Parameters:**


* `vec` The given vector.



**Returns:**

A double containing the value of the infinity norm. 





        

<hr>



### function norm\_inf 

_Compute the infinity norm._ 
```C++
KOKKOS_INLINE_FUNCTION double norm_inf (
    double const coord
) 
```



In case of scalar, the infinity norm returns the scalar.




**Parameters:**


* `coord` The given double.



**Returns:**

A double containing the value of the infinity norm. 





        

<hr>



### function norm\_inf 

_Compute the infinity norm for a Field._ 
```C++
template<class ExecSpace, class ElementType, class IdxRange>
inline double norm_inf (
    ExecSpace exec_space,
    ConstField< ElementType, IdxRange, typename ExecSpace::memory_space > function
) 
```





**Parameters:**


* `exec_space` The space on which the function is executed (CPU/GPU). 
* `function` The function whose norm is calculated. 



**Returns:**

A double containing the value of the infinity norm. 





        

<hr>



### function norm\_inf 

_Compute the infinity norm for a_ [_**VectorField**_](classVectorField.md) _._
```C++
template<class ExecSpace, class ElementType, class IdxRange, class VectorIndexSetType>
inline double norm_inf (
    ExecSpace exec_space,
    VectorConstField < ElementType, IdxRange, VectorIndexSetType, typename ExecSpace::memory_space > function
) 
```





**Parameters:**


* `exec_space` The space on which the function is executed (CPU/GPU). 
* `function` The function whose norm is calculated. 



**Returns:**

A double containing the value of the infinity norm. 





        

<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/math_tools/l_norm_tools.hpp`

