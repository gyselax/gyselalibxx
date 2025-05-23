

# Class Tensor

**template &lt;class ElementType, class ValidIndexSetFirstDim, class... ValidIndexSet&gt;**



[**ClassList**](annotated.md) **>** [**Tensor**](classTensor.md)



_A class representing a_ [_**Tensor**_](classTensor.md) _._[More...](#detailed-description)

* `#include <tensor.hpp>`



Inherits the following classes: [TensorCommon](classTensorCommon.md)
















## Public Types inherited from TensorCommon

See [TensorCommon](classTensorCommon.md)

| Type | Name |
| ---: | :--- |
| typedef typename DataStorageType::element\_type | [**element\_type**](classTensorCommon.md#typedef-element_type)  <br>_The type of the elements of the tensor._  |
| typedef ddc::detail::TypeSeq&lt; ValidIndexSet... &gt; | [**index\_set**](classTensorCommon.md#typedef-index_set)  <br>_The TensorIndexSet describing the possible indices._  |
| typedef ddc::type\_seq\_element\_t&lt; dim, [**index\_set**](classTensorCommon.md#typedef-index_set) &gt; | [**vector\_index\_set\_t**](classTensorCommon.md#typedef-vector_index_set_t)  <br>_A helper type alias to get all possible indices along a specified dimension._  |






































## Public Functions

| Type | Name |
| ---: | :--- |
|  KOKKOS\_DEFAULTED\_FUNCTION | [**Tensor**](#function-tensor-17) () = default<br>_Construct an uninitialised tensor object._  |
|  KOKKOS\_FUNCTION | [**Tensor**](#function-tensor-27) (ElementType fill\_value) <br>_Construct a tensor object initialised with a value._  |
|  KOKKOS\_FUNCTION | [**Tensor**](#function-tensor-37) (Params... elements) <br>_Construct a 1D tensor object by providing the elements that should be saved in the vector._  |
|  KOKKOS\_FUNCTION | [**Tensor**](#function-tensor-47) (Coord&lt; Dims... &gt; coord) noexcept<br>_Construct a 1D tensor object from a coordinate._  |
|  KOKKOS\_FUNCTION | [**Tensor**](#function-tensor-57) (const OTensorType & o\_tensor) noexcept<br>_Construct a tensor object by copying an existing compatible tensor. A tensor is compatible if it is defined on the same dimensions._  |
|  KOKKOS\_DEFAULTED\_FUNCTION | [**Tensor**](#function-tensor-67) ([**Tensor**](classTensor.md) const & o\_tensor) = default<br>_Construct a tensor object by copying an existing tensor of exactly the same type. This method can be called implicitly._  |
|  KOKKOS\_DEFAULTED\_FUNCTION | [**Tensor**](#function-tensor-77) ([**Tensor**](classTensor.md) && o\_tensor) = default<br>_Construct a tensor object by moving an existing tensor of exactly the same type. This method can be called implicitly._  |
|  KOKKOS\_DEFAULTED\_FUNCTION [**Tensor**](classTensor.md) & | [**operator=**](#function-operator) ([**Tensor**](classTensor.md) const & other) = default<br>_A copy assign operator._  |
|  KOKKOS\_DEFAULTED\_FUNCTION [**Tensor**](classTensor.md) & | [**operator=**](#function-operator_1) ([**Tensor**](classTensor.md) && other) = default<br>_A move assign operator._  |


## Public Functions inherited from TensorCommon

See [TensorCommon](classTensorCommon.md)

| Type | Name |
| ---: | :--- |
|  KOKKOS\_FUNCTION [**element\_type**](classTensorCommon.md#typedef-element_type) & | [**get**](classTensorCommon.md#function-get-12) () <br>_Get a modifiable reference to an element of the tensor._  |
|  KOKKOS\_FUNCTION [**element\_type**](classTensorCommon.md#typedef-element_type) const & | [**get**](classTensorCommon.md#function-get-22) () const<br>_Get an element of the tensor._  |
|  KOKKOS\_FUNCTION bool | [**operator!=**](classTensorCommon.md#function-operator) ([**TensorCommon**](classTensorCommon.md) const & o\_tensor) const<br>_An operator to compare one tensor to another elementwise._  |
|  KOKKOS\_FUNCTION [**TensorCommon**](classTensorCommon.md) & | [**operator\*=**](classTensorCommon.md#function-operator_1) (Oelement\_type val) <br>_An operator to multiply all the element of the current tensor by a value._  |
|  KOKKOS\_FUNCTION [**TensorCommon**](classTensorCommon.md) & | [**operator+=**](classTensorCommon.md#function-operator_2) ([**TensorCommon**](classTensorCommon.md) const & val) <br>_An operator to add two tensors elementwise._  |
|  KOKKOS\_FUNCTION [**TensorCommon**](classTensorCommon.md) & | [**operator-=**](classTensorCommon.md#function-operator-) ([**TensorCommon**](classTensorCommon.md) const & val) <br>_An operator to subtract one tensor from another elementwise._  |
|  KOKKOS\_FUNCTION [**TensorCommon**](classTensorCommon.md) & | [**operator/=**](classTensorCommon.md#function-operator_3) (Oelement\_type val) <br>_An operator to divide all the element of the current tensor by a value._  |
|  KOKKOS\_FUNCTION [**TensorCommon**](classTensorCommon.md) & | [**operator=**](classTensorCommon.md#function-operator_4) ([**TensorCommon**](classTensorCommon.md) const & other) <br>_A copy operator._  |
|  KOKKOS\_FUNCTION [**TensorCommon**](classTensorCommon.md) & | [**operator=**](classTensorCommon.md#function-operator_5) ([**TensorCommon**](classTensorCommon.md) && other) <br>_A move assign operator._  |
|  KOKKOS\_FUNCTION bool | [**operator==**](classTensorCommon.md#function-operator_6) ([**TensorCommon**](classTensorCommon.md) const & o\_tensor) const<br>_An operator to compare one tensor to another elementwise._  |




## Public Static Functions inherited from TensorCommon

See [TensorCommon](classTensorCommon.md)

| Type | Name |
| ---: | :--- |
|  KOKKOS\_FUNCTION constexpr std::size\_t | [**rank**](classTensorCommon.md#function-rank) () <br>_The rank of the tensor. This is equivalent to the number of indices required to access an element of the tensor._  |
|  KOKKOS\_FUNCTION constexpr std::size\_t | [**size**](classTensorCommon.md#function-size) () <br>_The size of the tensor. This is the number of elements in the tensor._  |












## Protected Attributes inherited from TensorCommon

See [TensorCommon](classTensorCommon.md)

| Type | Name |
| ---: | :--- |
|  DataStorageType | [**m\_data**](classTensorCommon.md#variable-m_data)  <br> |




## Protected Static Attributes inherited from TensorCommon

See [TensorCommon](classTensorCommon.md)

| Type | Name |
| ---: | :--- |
|  constexpr std::size\_t | [**s\_n\_elements**](classTensorCommon.md#variable-s_n_elements)   = `(ddc::type\_seq\_size\_v&lt;ValidIndexSet&gt; \* ...)`<br>_The number of elements in the mdspan._  |




























## Protected Functions inherited from TensorCommon

See [TensorCommon](classTensorCommon.md)

| Type | Name |
| ---: | :--- |
|  KOKKOS\_DEFAULTED\_FUNCTION | [**TensorCommon**](classTensorCommon.md#function-tensorcommon-13) () = default<br>_Construct an uninitialised tensor object._  |
|  KOKKOS\_DEFAULTED\_FUNCTION | [**TensorCommon**](classTensorCommon.md#function-tensorcommon-23) ([**TensorCommon**](classTensorCommon.md) const & o\_tensor) = default<br>_Construct a tensor object by copying an existing tensor of exactly the same type. This method can be called implicitly._  |
|  KOKKOS\_DEFAULTED\_FUNCTION | [**TensorCommon**](classTensorCommon.md#function-tensorcommon-33) ([**TensorCommon**](classTensorCommon.md) && o\_tensor) = default<br>_Move-construct a tensor object by copying an existing tensor of exactly the same type. This method can be called implicitly._  |






## Detailed Description




**Template parameters:**


* `ElementType` The type of the elements of the tensor (usually double/complex). 
* `ValidIndexSet` The indices that can be used along each dimension of the tensor. 




    
## Public Functions Documentation




### function Tensor [1/7]

_Construct an uninitialised tensor object._ 
```C++
KOKKOS_DEFAULTED_FUNCTION Tensor::Tensor () = default
```




<hr>



### function Tensor [2/7]

_Construct a tensor object initialised with a value._ 
```C++
inline explicit KOKKOS_FUNCTION Tensor::Tensor (
    ElementType fill_value
) 
```





**Parameters:**


* `fill_value` The value with which the tensor should be filled. 




        

<hr>



### function Tensor [3/7]

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



### function Tensor [4/7]

_Construct a 1D tensor object from a coordinate._ 
```C++
template<class... Dims>
inline explicit KOKKOS_FUNCTION Tensor::Tensor (
    Coord< Dims... > coord
) noexcept
```





**Parameters:**


* `coord` The coordinate. 




        

<hr>



### function Tensor [5/7]

_Construct a tensor object by copying an existing compatible tensor. A tensor is compatible if it is defined on the same dimensions._ 
```C++
template<class OTensorType, std::enable_if_t< is_tensor_type_v< OTensorType >, bool >>
inline explicit KOKKOS_FUNCTION Tensor::Tensor (
    const OTensorType & o_tensor
) noexcept
```





**Parameters:**


* `o_tensor` The tensor to be copied. 




        

<hr>



### function Tensor [6/7]

_Construct a tensor object by copying an existing tensor of exactly the same type. This method can be called implicitly._ 
```C++
KOKKOS_DEFAULTED_FUNCTION Tensor::Tensor (
    Tensor const & o_tensor
) = default
```





**Parameters:**


* `o_tensor` The tensor to be copied. 




        

<hr>



### function Tensor [7/7]

_Construct a tensor object by moving an existing tensor of exactly the same type. This method can be called implicitly._ 
```C++
KOKKOS_DEFAULTED_FUNCTION Tensor::Tensor (
    Tensor && o_tensor
) = default
```





**Parameters:**


* `o_tensor` The tensor to be copied. 




        

<hr>



### function operator= 

_A copy assign operator._ 
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

_A move assign operator._ 
```C++
KOKKOS_DEFAULTED_FUNCTION Tensor & Tensor::operator= (
    Tensor && other
) = default
```





**Parameters:**


* `other` The tensor to be copied. 



**Returns:**

A r-value reference to the current tensor. 





        

<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/data_types/tensor.hpp`

