

# Class IdentityTensor

**template &lt;class ElementType, class ValidIndexSetRow, class ValidIndexSetCol&gt;**



[**ClassList**](annotated.md) **>** [**IdentityTensor**](classIdentityTensor.md)



[More...](#detailed-description)

* `#include <static_tensors.hpp>`

















## Public Types

| Type | Name |
| ---: | :--- |
| typedef ElementType | [**element\_type**](#typedef-element_type)  <br>_The type of the elements of the tensor._  |
| typedef ddc::detail::TypeSeq&lt; ValidIndexSetRow, ValidIndexSetCol &gt; | [**index\_set**](#typedef-index_set)  <br>_The TensorIndexSet describing the possible indices._  |
| typedef std::conditional\_t&lt; dim==0, ValidIndexSetRow, ValidIndexSetCol &gt; | [**vector\_index\_set\_t**](#typedef-vector_index_set_t)  <br>_A helper type alias to get all possible indices along a specified dimension._  |




















## Public Functions

| Type | Name |
| ---: | :--- |
|  KOKKOS\_DEFAULTED\_FUNCTION | [**IdentityTensor**](#function-identitytensor) () = default<br>_Construct an uninitialised tensor object._  |


## Public Static Functions

| Type | Name |
| ---: | :--- |
|  constexpr KOKKOS\_FUNCTION ElementType | [**get**](#function-get) () <br>_Get an element of the tensor._  |
|  KOKKOS\_FUNCTION constexpr std::size\_t | [**rank**](#function-rank) () <br>_The rank of the tensor. This is equivalent to the number of indices required to access an element of the tensor._  |
|  KOKKOS\_FUNCTION constexpr std::size\_t | [**size**](#function-size) () <br>_The size of the tensor. This is the number of elements in the tensor._  |


























## Detailed Description


A class containing only static constexpr methods which describes the identity tensor. 

**Template parameters:**


* `ElementType` The type of the elements of the tensor (usually double/complex). 
* `ValidIndexSetRow` The indices of the rows. 
* `ValidIndexSetCol` The indices of the columns. 




    
## Public Types Documentation




### typedef element\_type 

_The type of the elements of the tensor._ 
```C++
using IdentityTensor< ElementType, ValidIndexSetRow, ValidIndexSetCol >::element_type =  ElementType;
```




<hr>



### typedef index\_set 

_The TensorIndexSet describing the possible indices._ 
```C++
using IdentityTensor< ElementType, ValidIndexSetRow, ValidIndexSetCol >::index_set =  ddc::detail::TypeSeq<ValidIndexSetRow, ValidIndexSetCol>;
```




<hr>



### typedef vector\_index\_set\_t 

_A helper type alias to get all possible indices along a specified dimension._ 
```C++
using IdentityTensor< ElementType, ValidIndexSetRow, ValidIndexSetCol >::vector_index_set_t =  std::conditional_t<dim == 0, ValidIndexSetRow, ValidIndexSetCol>;
```





**Template parameters:**


* `Dim` The dimension of interest (0 &lt;= dim &lt; [**rank()**](classIdentityTensor.md#function-rank)). 




        

<hr>
## Public Functions Documentation




### function IdentityTensor 

_Construct an uninitialised tensor object._ 
```C++
KOKKOS_DEFAULTED_FUNCTION IdentityTensor::IdentityTensor () = default
```




<hr>
## Public Static Functions Documentation




### function get 

_Get an element of the tensor._ 
```C++
template<class QueryTensorIndexElement>
static inline constexpr KOKKOS_FUNCTION ElementType IdentityTensor::get () 
```





**Template parameters:**


* `QueryIndexTag` A type describing the relevant index. 



**Returns:**

The relevant element of the tensor. 





        

<hr>



### function rank 

_The rank of the tensor. This is equivalent to the number of indices required to access an element of the tensor._ 
```C++
static inline KOKKOS_FUNCTION constexpr std::size_t IdentityTensor::rank () 
```





**Returns:**

The rank of the tensor. 





        

<hr>



### function size 

_The size of the tensor. This is the number of elements in the tensor._ 
```C++
static inline KOKKOS_FUNCTION constexpr std::size_t IdentityTensor::size () 
```





**Returns:**

The number of elements in the tensor. 





        

<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/data_types/static_tensors.hpp`

