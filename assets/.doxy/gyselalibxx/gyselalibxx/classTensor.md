

# Class Tensor

**template &lt;class ElementType, class... ValidIndexSet&gt;**



[**ClassList**](annotated.md) **>** [**Tensor**](classTensor.md)



_A class representing a_ [_**Tensor**_](classTensor.md) _._[More...](#detailed-description)

* `#include <tensor.hpp>`

















## Public Types

| Type | Name |
| ---: | :--- |
| typedef ElementType | [**element\_type**](#typedef-element_type)  <br>_The type of the elements of the tensor._  |
| typedef ddc::detail::TypeSeq&lt; ValidIndexSet... &gt; | [**index\_set**](#typedef-index_set)  <br>_The TensorIndexSet describing the possible indices._  |
| typedef ddc::type\_seq\_element\_t&lt; dim, [**index\_set**](classTensor.md#typedef-index_set) &gt; | [**vector\_index\_set\_t**](#typedef-vector_index_set_t)  <br>_A helper type alias to get all possible indices along a specified dimension._  |




















## Public Functions

| Type | Name |
| ---: | :--- |
|  KOKKOS\_DEFAULTED\_FUNCTION | [**Tensor**](#function-tensor-15) () = default<br>_Construct an uninitialised tensor object._  |
|  KOKKOS\_FUNCTION | [**Tensor**](#function-tensor-25) (ElementType fill\_value) <br>_Construct a tensor object initialised with a value._  |
|  KOKKOS\_FUNCTION | [**Tensor**](#function-tensor-35) (Params... elements) <br>_Construct a 1D tensor object by providing the elements that should be saved in the vector._  |
|  KOKKOS\_FUNCTION | [**Tensor**](#function-tensor-45) (Coord&lt; Dims... &gt; coord) <br>_Construct a 1D tensor object from a coordinate._  |
|  KOKKOS\_FUNCTION | [**Tensor**](#function-tensor-55) ([**Tensor**](classTensor.md)&lt; OElementType, ValidIndexSet... &gt; const & o\_tensor) <br>_Construct a tensor object by copying an existing tensor._  |
|  KOKKOS\_FUNCTION ElementType & | [**get**](#function-get-12) () <br>_Get a modifiable reference to an element of the tensor._  |
|  KOKKOS\_FUNCTION ElementType const & | [**get**](#function-get-22) () const<br>_Get an element of the tensor._  |
|  KOKKOS\_FUNCTION bool | [**operator!=**](#function-operator) ([**Tensor**](classTensor.md) const & o\_tensor) const<br>_An operator to compare one tensor to another elementwise._  |
|  KOKKOS\_FUNCTION [**Tensor**](classTensor.md) | [**operator\***](#function-operator_1) (OElementType val) const<br>_An operator to multiply all the element of the current tensor by a value._  |
|  KOKKOS\_FUNCTION [**Tensor**](classTensor.md) & | [**operator\*=**](#function-operator_2) (OElementType val) <br>_An operator to multiply all the element of the current tensor by a value._  |
|  KOKKOS\_FUNCTION [**Tensor**](classTensor.md) | [**operator+**](#function-operator_3) ([**Tensor**](classTensor.md) const & val) const<br>_An operator to add two tensors elementwise._  |
|  KOKKOS\_FUNCTION [**Tensor**](classTensor.md) & | [**operator+=**](#function-operator_4) ([**Tensor**](classTensor.md) const & val) <br>_An operator to add two tensors elementwise._  |
|  KOKKOS\_FUNCTION [**Tensor**](classTensor.md) | [**operator-**](#function-operator-) ([**Tensor**](classTensor.md) const & val) const<br>_An operator to subtract one tensor from another elementwise._  |
|  KOKKOS\_FUNCTION [**Tensor**](classTensor.md) | [**operator-**](#function-operator-_1) () const<br>_An operator to get the negation of a tensor elementwise._  |
|  KOKKOS\_FUNCTION [**Tensor**](classTensor.md) & | [**operator-=**](#function-operator-_2) ([**Tensor**](classTensor.md) const & val) <br>_An operator to subtract one tensor from another elementwise._  |
|  KOKKOS\_FUNCTION [**Tensor**](classTensor.md) | [**operator/**](#function-operator_5) (OElementType val) const<br>_An operator to multiply all the element of the current tensor by a value._  |
|  KOKKOS\_FUNCTION [**Tensor**](classTensor.md) & | [**operator/=**](#function-operator_6) (OElementType val) <br>_An operator to divide all the element of the current tensor by a value._  |
|  KOKKOS\_DEFAULTED\_FUNCTION [**Tensor**](classTensor.md) & | [**operator=**](#function-operator_7) ([**Tensor**](classTensor.md) const & other) = default<br>_A copy operator._  |
|  KOKKOS\_FUNCTION [**Tensor**](classTensor.md) & | [**operator=**](#function-operator_8) (Coord&lt; Dims... &gt; coord) <br>_A copy operator._  |
|  KOKKOS\_FUNCTION bool | [**operator==**](#function-operator_9) ([**Tensor**](classTensor.md) const & o\_tensor) const<br>_An operator to compare one tensor to another elementwise._  |


## Public Static Functions

| Type | Name |
| ---: | :--- |
|  KOKKOS\_FUNCTION constexpr std::size\_t | [**rank**](#function-rank) () <br>_The rank of the tensor. This is equivalent to the number of indices required to access an element of the tensor._  |
|  KOKKOS\_FUNCTION constexpr std::size\_t | [**size**](#function-size) () <br>_The size of the tensor. This is the number of elements in the tensor._  |


























## Detailed Description




**Template parameters:**


* `ElementType` The type of the elements of the tensor (usually double/complex). 
* `ValidIndexSet` The indices that can be used along each dimension of the tensor. 




    
## Public Types Documentation




### typedef element\_type 

_The type of the elements of the tensor._ 
```C++
using Tensor< ElementType, ValidIndexSet >::element_type =  ElementType;
```




<hr>



### typedef index\_set 

_The TensorIndexSet describing the possible indices._ 
```C++
using Tensor< ElementType, ValidIndexSet >::index_set =  ddc::detail::TypeSeq<ValidIndexSet...>;
```




<hr>



### typedef vector\_index\_set\_t 

_A helper type alias to get all possible indices along a specified dimension._ 
```C++
using Tensor< ElementType, ValidIndexSet >::vector_index_set_t =  ddc::type_seq_element_t<dim, index_set>;
```





**Template parameters:**


* `Dim` The dimension of interest (0 &lt;= dim &lt; [**rank()**](classTensor.md#function-rank)). 




        

<hr>
## Public Functions Documentation




### function Tensor [1/5]

_Construct an uninitialised tensor object._ 
```C++
KOKKOS_DEFAULTED_FUNCTION Tensor::Tensor () = default
```




<hr>



### function Tensor [2/5]

_Construct a tensor object initialised with a value._ 
```C++
inline explicit KOKKOS_FUNCTION Tensor::Tensor (
    ElementType fill_value
) 
```





**Parameters:**


* `fill_value` The value with which the tensor should be filled. 




        

<hr>



### function Tensor [3/5]

_Construct a 1D tensor object by providing the elements that should be saved in the vector._ 
```C++
template<class... Params, class, class>
inline explicit KOKKOS_FUNCTION Tensor::Tensor (
    Params... elements
) 
```





**Parameters:**


* `elements` The elements of the tensor. 




        

<hr>



### function Tensor [4/5]

_Construct a 1D tensor object from a coordinate._ 
```C++
template<class... Dims>
inline explicit KOKKOS_FUNCTION Tensor::Tensor (
    Coord< Dims... > coord
) 
```





**Parameters:**


* `coord` The coordinate. 




        

<hr>



### function Tensor [5/5]

_Construct a tensor object by copying an existing tensor._ 
```C++
template<class OElementType>
inline explicit KOKKOS_FUNCTION Tensor::Tensor (
    Tensor < OElementType, ValidIndexSet... > const & o_tensor
) 
```





**Parameters:**


* `o_tensor` The tensor to be copied. 




        

<hr>



### function get [1/2]

_Get a modifiable reference to an element of the tensor._ 
```C++
template<class QueryTensorIndexElement>
inline KOKKOS_FUNCTION ElementType & Tensor::get () 
```





**Template parameters:**


* `QueryIndexTag` A type describing the relevant index. 



**Returns:**

The relevant element of the tensor. 





        

<hr>



### function get [2/2]

_Get an element of the tensor._ 
```C++
template<class QueryTensorIndexElement>
inline KOKKOS_FUNCTION ElementType const & Tensor::get () const
```





**Template parameters:**


* `QueryIndexTag` A type describing the relevant index. 



**Returns:**

The relevant element of the tensor. 





        

<hr>



### function operator!= 

_An operator to compare one tensor to another elementwise._ 
```C++
inline KOKKOS_FUNCTION bool Tensor::operator!= (
    Tensor const & o_tensor
) const
```





**Parameters:**


* `o_tensor` The tensor that should be compared with the current tensor. 



**Returns:**

False if the tensors are equal, true otherwise. 





        

<hr>



### function operator\* 

_An operator to multiply all the element of the current tensor by a value._ 
```C++
template<class OElementType>
inline KOKKOS_FUNCTION Tensor Tensor::operator* (
    OElementType val
) const
```





**Parameters:**


* `val` The value by which the elements should be multiplied. 



**Returns:**

A new tensor containing the result of the multiplication. 





        

<hr>



### function operator\*= 

_An operator to multiply all the element of the current tensor by a value._ 
```C++
template<class OElementType>
inline KOKKOS_FUNCTION Tensor & Tensor::operator*= (
    OElementType val
) 
```





**Parameters:**


* `val` The value by which the elements should be multiplied. 



**Returns:**

A reference to the current modified tensor. 





        

<hr>



### function operator+ 

_An operator to add two tensors elementwise._ 
```C++
inline KOKKOS_FUNCTION Tensor Tensor::operator+ (
    Tensor const & val
) const
```





**Parameters:**


* `val` The tensor that should be added to the current tensor. 



**Returns:**

A new tensor containing the result of the addition. 





        

<hr>



### function operator+= 

_An operator to add two tensors elementwise._ 
```C++
inline KOKKOS_FUNCTION Tensor & Tensor::operator+= (
    Tensor const & val
) 
```





**Parameters:**


* `val` The tensor that should be added to the current tensor. 



**Returns:**

A reference to the current modified tensor. 





        

<hr>



### function operator- 

_An operator to subtract one tensor from another elementwise._ 
```C++
inline KOKKOS_FUNCTION Tensor Tensor::operator- (
    Tensor const & val
) const
```





**Parameters:**


* `val` The tensor that should be subtracted from the current tensor. 



**Returns:**

A new tensor containing the result of the subtraction. 





        

<hr>



### function operator- 

_An operator to get the negation of a tensor elementwise._ 
```C++
inline KOKKOS_FUNCTION Tensor Tensor::operator- () const
```





**Returns:**

A new tensor containing the result of the subtraction. 





        

<hr>



### function operator-= 

_An operator to subtract one tensor from another elementwise._ 
```C++
inline KOKKOS_FUNCTION Tensor & Tensor::operator-= (
    Tensor const & val
) 
```





**Parameters:**


* `val` The tensor that should be subtracted from the current tensor. 



**Returns:**

A reference to the current modified tensor. 





        

<hr>



### function operator/ 

_An operator to multiply all the element of the current tensor by a value._ 
```C++
template<class OElementType>
inline KOKKOS_FUNCTION Tensor Tensor::operator/ (
    OElementType val
) const
```





**Parameters:**


* `val` The value by which the elements should be multiplied. 



**Returns:**

A new tensor containing the result of the multiplication. 





        

<hr>



### function operator/= 

_An operator to divide all the element of the current tensor by a value._ 
```C++
template<class OElementType>
inline KOKKOS_FUNCTION Tensor & Tensor::operator/= (
    OElementType val
) 
```





**Parameters:**


* `val` The value by which the elements should be multiplied. 



**Returns:**

A reference to the current modified tensor. 





        

<hr>



### function operator= 

_A copy operator._ 
```C++
KOKKOS_DEFAULTED_FUNCTION Tensor & Tensor::operator= (
    Tensor const & other
) = default
```





**Parameters:**


* `other` The tensor to be copied. 



**Returns:**

A reference to the current tensor. 





        

<hr>



### function operator= 

_A copy operator._ 
```C++
template<class... Dims>
inline KOKKOS_FUNCTION Tensor & Tensor::operator= (
    Coord< Dims... > coord
) 
```





**Parameters:**


* `coord` The coordinate to be copied into the vector. 



**Returns:**

A reference to the current tensor. 





        

<hr>



### function operator== 

_An operator to compare one tensor to another elementwise._ 
```C++
inline KOKKOS_FUNCTION bool Tensor::operator== (
    Tensor const & o_tensor
) const
```





**Parameters:**


* `o_tensor` The tensor that should be compared with the current tensor. 



**Returns:**

True if the tensors are equal, false otherwise. 





        

<hr>
## Public Static Functions Documentation




### function rank 

_The rank of the tensor. This is equivalent to the number of indices required to access an element of the tensor._ 
```C++
static inline KOKKOS_FUNCTION constexpr std::size_t Tensor::rank () 
```





**Returns:**

The rank of the tensor. 





        

<hr>



### function size 

_The size of the tensor. This is the number of elements in the tensor._ 
```C++
static inline KOKKOS_FUNCTION constexpr std::size_t Tensor::size () 
```





**Returns:**

The number of elements in the tensor. 





        

<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/data_types/tensor.hpp`

