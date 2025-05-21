

# Class IdentityCoordinateChange

**template &lt;class ArgBasis, class ResultBasis&gt;**



[**ClassList**](annotated.md) **>** [**IdentityCoordinateChange**](classIdentityCoordinateChange.md)



_A class describing an identity transformation._ [More...](#detailed-description)

* `#include <identity_coordinate_change.hpp>`

















## Public Types

| Type | Name |
| ---: | :--- |
| typedef ddcHelper::to\_coord\_t&lt; ArgBasis &gt; | [**CoordArg**](#typedef-coordarg)  <br>_The type of the argument of the function described by this mapping._  |
| typedef ddcHelper::to\_coord\_t&lt; ArgBasis &gt; | [**CoordJacobian**](#typedef-coordjacobian)  <br>_The type of the coordinate that can be used to evaluate the Jacobian of this mapping._  |
| typedef ddcHelper::to\_coord\_t&lt; ResultBasis &gt; | [**CoordResult**](#typedef-coordresult)  <br>_The type of the result of the function described by this mapping._  |




















## Public Functions

| Type | Name |
| ---: | :--- |
|  KOKKOS\_INLINE\_FUNCTION [**IdentityCoordinateChange**](classIdentityCoordinateChange.md)&lt; ResultBasis, ArgBasis &gt; | [**get\_inverse\_mapping**](#function-get_inverse_mapping) () const<br>_Get the inverse mapping._  |
|  KOKKOS\_INLINE\_FUNCTION double | [**inv\_jacobian**](#function-inv_jacobian) ([**CoordArg**](classIdentityCoordinateChange.md#typedef-coordarg) const & coord) const<br>_Compute the Jacobian, the determinant of the Jacobian matrix of the mapping._  |
|  KOKKOS\_INLINE\_FUNCTION double | [**inv\_jacobian\_component**](#function-inv_jacobian_component) ([**CoordArg**](classIdentityCoordinateChange.md#typedef-coordarg) const & coord) const<br>_Compute the (i,j) coefficient of the Jacobian matrix._  |
|  KOKKOS\_FUNCTION [**DTensor**](classTensor.md)&lt; ResultBasis, ArgBasisCov &gt; | [**inv\_jacobian\_matrix**](#function-inv_jacobian_matrix) ([**CoordArg**](classIdentityCoordinateChange.md#typedef-coordarg) const & coord) const<br>_Compute full Jacobian matrix._  |
|  KOKKOS\_INLINE\_FUNCTION double | [**jacobian**](#function-jacobian) ([**CoordArg**](classIdentityCoordinateChange.md#typedef-coordarg) const & coord) const<br>_Compute the Jacobian, the determinant of the Jacobian matrix of the mapping._  |
|  KOKKOS\_INLINE\_FUNCTION double | [**jacobian\_component**](#function-jacobian_component) ([**CoordArg**](classIdentityCoordinateChange.md#typedef-coordarg) const & coord) const<br>_Compute the (i,j) coefficient of the Jacobian matrix._  |
|  KOKKOS\_FUNCTION [**DTensor**](classTensor.md)&lt; ResultBasis, ArgBasisCov &gt; | [**jacobian\_matrix**](#function-jacobian_matrix) ([**CoordArg**](classIdentityCoordinateChange.md#typedef-coordarg) const & coord) const<br>_Compute full Jacobian matrix._  |
|  KOKKOS\_FUNCTION [**CoordResult**](classIdentityCoordinateChange.md#typedef-coordresult) | [**operator()**](#function-operator) (Coord&lt; ArgDims... &gt; const & coord) const<br>_Convert the coordinate in the argument basis to the equivalent coordinate in the result basis._  |




























## Detailed Description


It is not expected that this class appear in a final simulation, but it may be useful when debugging vector calculations on general coordinates.




**Template parameters:**


* `ArgBasis` A VectorIndexSet containing the continuous dimensions on which the argument is described. 
* `ResultBasis` A VectorIndexSet containing the continuous dimensions on which the result is described. 




    
## Public Types Documentation




### typedef CoordArg 

_The type of the argument of the function described by this mapping._ 
```C++
using IdentityCoordinateChange< ArgBasis, ResultBasis >::CoordArg =  ddcHelper::to_coord_t<ArgBasis>;
```




<hr>



### typedef CoordJacobian 

_The type of the coordinate that can be used to evaluate the Jacobian of this mapping._ 
```C++
using IdentityCoordinateChange< ArgBasis, ResultBasis >::CoordJacobian =  ddcHelper::to_coord_t<ArgBasis>;
```




<hr>



### typedef CoordResult 

_The type of the result of the function described by this mapping._ 
```C++
using IdentityCoordinateChange< ArgBasis, ResultBasis >::CoordResult =  ddcHelper::to_coord_t<ResultBasis>;
```




<hr>
## Public Functions Documentation




### function get\_inverse\_mapping 

_Get the inverse mapping._ 
```C++
inline KOKKOS_INLINE_FUNCTION IdentityCoordinateChange < ResultBasis, ArgBasis > IdentityCoordinateChange::get_inverse_mapping () const
```





**Returns:**

The inverse mapping. 





        

<hr>



### function inv\_jacobian 

_Compute the Jacobian, the determinant of the Jacobian matrix of the mapping._ 
```C++
inline KOKKOS_INLINE_FUNCTION double IdentityCoordinateChange::inv_jacobian (
    CoordArg const & coord
) const
```





**Parameters:**


* `coord` The coordinate where we evaluate the Jacobian.



**Returns:**

A double with the value of the determinant of the Jacobian matrix. 





        

<hr>



### function inv\_jacobian\_component 

_Compute the (i,j) coefficient of the Jacobian matrix._ 
```C++
template<class IndexTag1, class IndexTag2>
inline KOKKOS_INLINE_FUNCTION double IdentityCoordinateChange::inv_jacobian_component (
    CoordArg const & coord
) const
```





**Parameters:**


* `coord` The coordinate where we evaluate the Jacobian matrix.



**Returns:**

A double with the value of the (i,j) coefficient of the Jacobian matrix. 





        

<hr>



### function inv\_jacobian\_matrix 

_Compute full Jacobian matrix._ 
```C++
inline KOKKOS_FUNCTION DTensor < ResultBasis, ArgBasisCov > IdentityCoordinateChange::inv_jacobian_matrix (
    CoordArg const & coord
) const
```





**Parameters:**


* `coord` The coordinate where we evaluate the Jacobian matrix. 



**Returns:**

The Jacobian matrix. 





        

<hr>



### function jacobian 

_Compute the Jacobian, the determinant of the Jacobian matrix of the mapping._ 
```C++
inline KOKKOS_INLINE_FUNCTION double IdentityCoordinateChange::jacobian (
    CoordArg const & coord
) const
```





**Parameters:**


* `coord` The coordinate where we evaluate the Jacobian.



**Returns:**

A double with the value of the determinant of the Jacobian matrix. 





        

<hr>



### function jacobian\_component 

_Compute the (i,j) coefficient of the Jacobian matrix._ 
```C++
template<class IndexTag1, class IndexTag2>
inline KOKKOS_INLINE_FUNCTION double IdentityCoordinateChange::jacobian_component (
    CoordArg const & coord
) const
```





**Parameters:**


* `coord` The coordinate where we evaluate the Jacobian matrix.



**Returns:**

A double with the value of the (i,j) coefficient of the Jacobian matrix. 





        

<hr>



### function jacobian\_matrix 

_Compute full Jacobian matrix._ 
```C++
inline KOKKOS_FUNCTION DTensor < ResultBasis, ArgBasisCov > IdentityCoordinateChange::jacobian_matrix (
    CoordArg const & coord
) const
```





**Parameters:**


* `coord` The coordinate where we evaluate the Jacobian matrix. 



**Returns:**

The Jacobian matrix. 





        

<hr>



### function operator() 

_Convert the coordinate in the argument basis to the equivalent coordinate in the result basis._ 
```C++
template<class... ArgDims>
inline KOKKOS_FUNCTION CoordResult IdentityCoordinateChange::operator() (
    Coord< ArgDims... > const & coord
) const
```





**Parameters:**


* `coord` The coordinate to be converted.



**Returns:**

The equivalent coordinate. 





        

<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/coord_transformations/identity_coordinate_change.hpp`

