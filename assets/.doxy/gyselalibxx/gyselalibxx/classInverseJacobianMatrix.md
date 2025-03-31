

# Class InverseJacobianMatrix

**template &lt;class Mapping, class PositionCoordinate&gt;**



[**ClassList**](annotated.md) **>** [**InverseJacobianMatrix**](classInverseJacobianMatrix.md)



[More...](#detailed-description)

* `#include <inverse_jacobian_matrix.hpp>`

















## Public Types

| Type | Name |
| ---: | :--- |
| typedef [**DTensor**](classTensor.md)&lt; ValidArgIndices, vector\_index\_set\_dual\_t&lt; ValidResultIndices &gt; &gt; | [**InverseJacobianTensor**](#typedef-inversejacobiantensor)  <br>_The type of the tensor representing the inverse Jacobian._  |




















## Public Functions

| Type | Name |
| ---: | :--- |
|  KOKKOS\_FUNCTION | [**InverseJacobianMatrix**](#function-inversejacobianmatrix) (Mapping const & mapping) <br>_A constructor for the_ [_**InverseJacobianMatrix**_](classInverseJacobianMatrix.md) _._ |
|  KOKKOS\_INLINE\_FUNCTION double | [**inv\_jacobian\_component**](#function-inv_jacobian_component) (PositionCoordinate const & coord) const<br>_Compute the (i,j) coefficient of the inverse Jacobian matrix._  |
|  KOKKOS\_INLINE\_FUNCTION [**InverseJacobianTensor**](classInverseJacobianMatrix.md#typedef-inversejacobiantensor) | [**operator()**](#function-operator) (PositionCoordinate const & coord) const<br>_Compute full inverse Jacobian matrix._  |




























## Detailed Description


A class to calculate the inverse of the Jacobian matrix. If specialised methods are available then these are used by the class. Otherwise the inverse is calculated from the Jacobian matrix.




**Template parameters:**


* `Mapping` The mapping whose inverse we are interested in. 
* `PositionCoordinate` The coordinate system in which the inverse should be calculated. 




    
## Public Types Documentation




### typedef InverseJacobianTensor 

_The type of the tensor representing the inverse Jacobian._ 
```C++
using InverseJacobianMatrix< Mapping, PositionCoordinate >::InverseJacobianTensor =  DTensor<ValidArgIndices, vector_index_set_dual_t<ValidResultIndices> >;
```




<hr>
## Public Functions Documentation




### function InverseJacobianMatrix 

_A constructor for the_ [_**InverseJacobianMatrix**_](classInverseJacobianMatrix.md) _._
```C++
inline explicit KOKKOS_FUNCTION InverseJacobianMatrix::InverseJacobianMatrix (
    Mapping const & mapping
) 
```





**Parameters:**


* `mapping` The mapping whose inverse we are interested in. 




        

<hr>



### function inv\_jacobian\_component 

_Compute the (i,j) coefficient of the inverse Jacobian matrix._ 
```C++
template<class IndexTag1, class IndexTag2>
inline KOKKOS_INLINE_FUNCTION double InverseJacobianMatrix::inv_jacobian_component (
    PositionCoordinate const & coord
) const
```



Be careful because not all mappings are invertible, especially at the centre point.




**Parameters:**


* `coord` The coordinate where we evaluate the inverse Jacobian matrix.



**Returns:**

A double with the value of the (i,j) coefficient of the inverse Jacobian matrix. 





        

<hr>



### function operator() 

_Compute full inverse Jacobian matrix._ 
```C++
inline KOKKOS_INLINE_FUNCTION InverseJacobianTensor InverseJacobianMatrix::operator() (
    PositionCoordinate const & coord
) const
```



For some computations, we need the complete inverse Jacobian matrix or just the coefficients. The coefficients can be given independently with the function inv\_jacobian\_component.




**Parameters:**


* `coord` The coordinate where we evaluate the Jacobian matrix. 



**Returns:**

The inverse Jacobian matrix returned. 





        

<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/mapping/inverse_jacobian_matrix.hpp`

