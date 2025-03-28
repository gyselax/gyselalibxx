

# Class CombinedMapping

**template &lt;class Mapping1, class Mapping2&gt;**



[**ClassList**](annotated.md) **>** [**CombinedMapping**](classCombinedMapping.md)



_A class which describes a mapping which is constructed by combining two mappings. Let us denote Mapping1 as_  _and Mapping2 as_ _then this mapping represents:_ _._[More...](#detailed-description)

* `#include <combined_mapping.hpp>`

















## Public Types

| Type | Name |
| ---: | :--- |
| typedef typename Mapping2::CoordArg | [**CoordArg**](#typedef-coordarg)  <br>_The type of the argument of the function described by this mapping._  |
| typedef typename Mapping2::CoordResult | [**CoordJacobian**](#typedef-coordjacobian)  <br>_The coordinate system on which the Jacobian is described._  |
| typedef typename Mapping1::CoordResult | [**CoordResult**](#typedef-coordresult)  <br>_The type of the result of the function described by this mapping._  |
| typedef [**DTensor**](classTensor.md)&lt; ddc::to\_type\_seq\_t&lt; [**CoordArg**](classCombinedMapping.md#typedef-coordarg) &gt;, vector\_index\_set\_dual\_t&lt; ddc::to\_type\_seq\_t&lt; [**CoordResult**](classCombinedMapping.md#typedef-coordresult) &gt; &gt; &gt; | [**InvJacobianMatrixType**](#typedef-invjacobianmatrixtype)  <br>_The type of the inverse Jacobian matrix._  |
| typedef [**DTensor**](classTensor.md)&lt; ddc::to\_type\_seq\_t&lt; [**CoordResult**](classCombinedMapping.md#typedef-coordresult) &gt;, vector\_index\_set\_dual\_t&lt; ddc::to\_type\_seq\_t&lt; [**CoordArg**](classCombinedMapping.md#typedef-coordarg) &gt; &gt; &gt; | [**JacobianMatrixType**](#typedef-jacobianmatrixtype)  <br>_The type of the Jacobian matrix._  |




















## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**CombinedMapping**](#function-combinedmapping-12) (Map1 mapping\_1, Mapping2 mapping\_2, double epsilon) <br>_Build a_ [_**CombinedMapping**_](classCombinedMapping.md) _from the component mappings. This constructor should be used if the inverse jacobian of the first mapping, or the jacobian of the second mapping cannot be evaluated at the O-point._ |
|   | [**CombinedMapping**](#function-combinedmapping-22) (Map1 mapping\_1, Mapping2 mapping\_2) <br>_Build a_ [_**CombinedMapping**_](classCombinedMapping.md) _from the component mappings. This constructor should be used if both mappings can be safely evaluated at all points in space._ |
|  KOKKOS\_INLINE\_FUNCTION Mapping const & | [**get**](#function-get) () const<br>_Get access to one of the internal mappings._  |
|  KOKKOS\_INLINE\_FUNCTION double | [**inv\_jacobian**](#function-inv_jacobian) ([**CoordJacobian**](classCombinedMapping.md#typedef-coordjacobian) const & coord\_rtheta) const<br>_Compute the determinant of the Jacobian matrix._  |
|  KOKKOS\_INLINE\_FUNCTION double | [**inv\_jacobian\_11**](#function-inv_jacobian_11) ([**CoordJacobian**](classCombinedMapping.md#typedef-coordjacobian) const & coord\_rtheta) const<br>_Compute the (1,1) coefficient of the Jacobian matrix._  |
|  KOKKOS\_INLINE\_FUNCTION double | [**inv\_jacobian\_12**](#function-inv_jacobian_12) ([**CoordJacobian**](classCombinedMapping.md#typedef-coordjacobian) const & coord\_rtheta) const<br>_Compute the (1,2) coefficient of the Jacobian matrix._  |
|  KOKKOS\_INLINE\_FUNCTION double | [**inv\_jacobian\_21**](#function-inv_jacobian_21) ([**CoordJacobian**](classCombinedMapping.md#typedef-coordjacobian) const & coord\_rtheta) const<br>_Compute the (2,1) coefficient of the Jacobian matrix._  |
|  KOKKOS\_INLINE\_FUNCTION double | [**inv\_jacobian\_22**](#function-inv_jacobian_22) ([**CoordJacobian**](classCombinedMapping.md#typedef-coordjacobian) const & coord\_rtheta) const<br>_Compute the (2,2) coefficient of the Jacobian matrix._  |
|  KOKKOS\_FUNCTION [**InvJacobianMatrixType**](classCombinedMapping.md#typedef-invjacobianmatrixtype) | [**inv\_jacobian\_matrix**](#function-inv_jacobian_matrix) ([**CoordJacobian**](classCombinedMapping.md#typedef-coordjacobian) const & coord) const<br>_Compute the full inverse Jacobian matrix._  |
|  KOKKOS\_INLINE\_FUNCTION double | [**jacobian**](#function-jacobian) ([**CoordJacobian**](classCombinedMapping.md#typedef-coordjacobian) const & coord\_rtheta) const<br>_Compute the determinant of the Jacobian matrix._  |
|  KOKKOS\_INLINE\_FUNCTION double | [**jacobian\_component**](#function-jacobian_component) ([**CoordJacobian**](classCombinedMapping.md#typedef-coordjacobian) const & coord\_rtheta) const<br>_Compute the (i,j) component of the Jacobian matrix._  |
|  KOKKOS\_INLINE\_FUNCTION [**JacobianMatrixType**](classCombinedMapping.md#typedef-jacobianmatrixtype) | [**jacobian\_matrix**](#function-jacobian_matrix) ([**CoordJacobian**](classCombinedMapping.md#typedef-coordjacobian) const & coord) const<br>_Compute full Jacobian matrix._  |
|  [**CoordResult**](classCombinedMapping.md#typedef-coordresult) | [**operator()**](#function-operator) ([**CoordArg**](classCombinedMapping.md#typedef-coordarg) coord) <br>_Convert the argument coordinate to the equivalent result coordinate._  |




























## Detailed Description


There are therefore 3 domains in this calculation ,  and  with  and 


The functions in this mapping are defined on the coordinate system associated with the domain . 


    
## Public Types Documentation




### typedef CoordArg 

_The type of the argument of the function described by this mapping._ 
```C++
using CombinedMapping< Mapping1, Mapping2 >::CoordArg =  typename Mapping2::CoordArg;
```




<hr>



### typedef CoordJacobian 

_The coordinate system on which the Jacobian is described._ 
```C++
using CombinedMapping< Mapping1, Mapping2 >::CoordJacobian =  typename Mapping2::CoordResult;
```




<hr>



### typedef CoordResult 

_The type of the result of the function described by this mapping._ 
```C++
using CombinedMapping< Mapping1, Mapping2 >::CoordResult =  typename Mapping1::CoordResult;
```




<hr>



### typedef InvJacobianMatrixType 

_The type of the inverse Jacobian matrix._ 
```C++
using CombinedMapping< Mapping1, Mapping2 >::InvJacobianMatrixType =  DTensor< ddc::to_type_seq_t<CoordArg>, vector_index_set_dual_t<ddc::to_type_seq_t<CoordResult> >>;
```




<hr>



### typedef JacobianMatrixType 

_The type of the Jacobian matrix._ 
```C++
using CombinedMapping< Mapping1, Mapping2 >::JacobianMatrixType =  DTensor< ddc::to_type_seq_t<CoordResult>, vector_index_set_dual_t<ddc::to_type_seq_t<CoordArg> >>;
```




<hr>
## Public Functions Documentation




### function CombinedMapping [1/2]

_Build a_ [_**CombinedMapping**_](classCombinedMapping.md) _from the component mappings. This constructor should be used if the inverse jacobian of the first mapping, or the jacobian of the second mapping cannot be evaluated at the O-point._
```C++
template<class Map1, std::enable_if_t<(has_singular_o_point_inv_jacobian_v< Map1 >)||(has_singular_o_point_inv_jacobian_v< InverseMapping2 >), bool >>
inline CombinedMapping::CombinedMapping (
    Map1 mapping_1,
    Mapping2 mapping_2,
    double epsilon
) 
```





**Parameters:**


* `mapping_1` The first mapping. 
* `mapping_2` The second mapping. 
* `epsilon` The parameter  which determines when a point is close enough to the central O-point for linearisation to be required when calculating the inverse of the Jacobian. The Jacobian is linearised on . 




        

<hr>



### function CombinedMapping [2/2]

_Build a_ [_**CombinedMapping**_](classCombinedMapping.md) _from the component mappings. This constructor should be used if both mappings can be safely evaluated at all points in space._
```C++
template<class Map1, std::enable_if_t< !((has_singular_o_point_inv_jacobian_v< Map1 >)||(has_singular_o_point_inv_jacobian_v< InverseMapping2 >)), bool >>
inline CombinedMapping::CombinedMapping (
    Map1 mapping_1,
    Mapping2 mapping_2
) 
```





**Parameters:**


* `mapping_1` The first mapping. 
* `mapping_2` The second mapping. 




        

<hr>



### function get 

_Get access to one of the internal mappings._ 
```C++
template<class Mapping>
inline KOKKOS_INLINE_FUNCTION Mapping const & CombinedMapping::get () const
```





**Template parameters:**


* `Mapping` The mapping that we want to access. 



**Returns:**

A constant reference to the internal mapping. 





        

<hr>



### function inv\_jacobian 

_Compute the determinant of the Jacobian matrix._ 
```C++
inline KOKKOS_INLINE_FUNCTION double CombinedMapping::inv_jacobian (
    CoordJacobian const & coord_rtheta
) const
```





**See also:** [**inv\_jacobian\_matrix**](classCombinedMapping.md#function-inv_jacobian_matrix) 


**Parameters:**


* `coord_rtheta` The coordinate where we evaluate the Jacobian matrix. 



**Returns:**

The determinant of the Jacobian matrix. 





        

<hr>



### function inv\_jacobian\_11 

_Compute the (1,1) coefficient of the Jacobian matrix._ 
```C++
inline KOKKOS_INLINE_FUNCTION double CombinedMapping::inv_jacobian_11 (
    CoordJacobian const & coord_rtheta
) const
```





**See also:** [**inv\_jacobian\_matrix**](classCombinedMapping.md#function-inv_jacobian_matrix) 


**Parameters:**


* `coord_rtheta` The coordinate where we evaluate the Jacobian matrix. 



**Returns:**

The (1,1) coefficient of the Jacobian matrix. 





        

<hr>



### function inv\_jacobian\_12 

_Compute the (1,2) coefficient of the Jacobian matrix._ 
```C++
inline KOKKOS_INLINE_FUNCTION double CombinedMapping::inv_jacobian_12 (
    CoordJacobian const & coord_rtheta
) const
```





**See also:** [**inv\_jacobian\_matrix**](classCombinedMapping.md#function-inv_jacobian_matrix) 


**Parameters:**


* `coord_rtheta` The coordinate where we evaluate the Jacobian matrix. 



**Returns:**

The (1,2) coefficient of the Jacobian matrix. 





        

<hr>



### function inv\_jacobian\_21 

_Compute the (2,1) coefficient of the Jacobian matrix._ 
```C++
inline KOKKOS_INLINE_FUNCTION double CombinedMapping::inv_jacobian_21 (
    CoordJacobian const & coord_rtheta
) const
```





**See also:** [**inv\_jacobian\_matrix**](classCombinedMapping.md#function-inv_jacobian_matrix) 


**Parameters:**


* `coord_rtheta` The coordinate where we evaluate the Jacobian matrix. 



**Returns:**

The (2,1) coefficient of the Jacobian matrix. 





        

<hr>



### function inv\_jacobian\_22 

_Compute the (2,2) coefficient of the Jacobian matrix._ 
```C++
inline KOKKOS_INLINE_FUNCTION double CombinedMapping::inv_jacobian_22 (
    CoordJacobian const & coord_rtheta
) const
```





**See also:** [**inv\_jacobian\_matrix**](classCombinedMapping.md#function-inv_jacobian_matrix) 


**Parameters:**


* `coord_rtheta` The coordinate where we evaluate the Jacobian matrix. 



**Returns:**

The (2,2) coefficient of the Jacobian matrix. 





        

<hr>



### function inv\_jacobian\_matrix 

_Compute the full inverse Jacobian matrix._ 
```C++
inline KOKKOS_FUNCTION InvJacobianMatrixType CombinedMapping::inv_jacobian_matrix (
    CoordJacobian const & coord
) const
```



If one of the mappings is singular then this function linearises the inverse Jacobian between the O-point and . The inverse Jacobian at the O-point is calculated using the class [**InvJacobianOPoint**](classInvJacobianOPoint.md) in this case. The inverse Jacobian at a point where the matrices are not singular is calculated using non\_singular\_inverse\_jacobian\_matrix




**Parameters:**


* `coord` The coordinate where we evaluate the inverse Jacobian matrix. 



**Returns:**

The calculated inverse Jacobian matrix.




**See also:** non\_singular\_inverse\_jacobian\_matrix 



        

<hr>



### function jacobian 

_Compute the determinant of the Jacobian matrix._ 
```C++
inline KOKKOS_INLINE_FUNCTION double CombinedMapping::jacobian (
    CoordJacobian const & coord_rtheta
) const
```





**See also:** [**jacobian\_matrix**](classCombinedMapping.md#function-jacobian_matrix) 


**Parameters:**


* `coord_rtheta` The coordinate where we evaluate the Jacobian matrix. 



**Returns:**

The determinant of the Jacobian matrix. 





        

<hr>



### function jacobian\_component 

_Compute the (i,j) component of the Jacobian matrix._ 
```C++
template<class IndexTag1, class IndexTag2>
inline KOKKOS_INLINE_FUNCTION double CombinedMapping::jacobian_component (
    CoordJacobian const & coord_rtheta
) const
```





**Parameters:**


* `coord_rtheta` The coordinate where we evaluate the Jacobian matrix. 



**Returns:**

The (i,j) component of the Jacobian matrix. 





        

<hr>



### function jacobian\_matrix 

_Compute full Jacobian matrix._ 
```C++
inline KOKKOS_INLINE_FUNCTION JacobianMatrixType CombinedMapping::jacobian_matrix (
    CoordJacobian const & coord
) const
```



For some computations, we need the complete Jacobian matrix or just the coefficients. This is calculated by combining the Jacobian matrices of the 2 curvilinear mappings.


  


The Jacobians that are used for the calculation must be mappings from  so they can be calculated on the correct coordinate system.




**Parameters:**


* `coord` The coordinate where we evaluate the Jacobian matrix. 



**Returns:**

The calculated Jacobian matrix. 





        

<hr>



### function operator() 

_Convert the argument coordinate to the equivalent result coordinate._ 
```C++
inline CoordResult CombinedMapping::operator() (
    CoordArg coord
) 
```





**Parameters:**


* `coord` The coordinate to be converted.



**Returns:**

The equivalent coordinate. 





        

<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/mapping/combined_mapping.hpp`

