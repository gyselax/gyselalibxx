

# File tensor\_common.hpp



[**FileList**](files.md) **>** [**data\_types**](dir_eaa769653453aaefd8cc10e98e9bb3eb.md) **>** [**tensor\_common.hpp**](tensor__common_8hpp.md)

[Go to the source code of this file](tensor__common_8hpp_source.md)



* `#include <array>`
* `#include <cassert>`
* `#include <ddc/ddc.hpp>`
* `#include "ddc_aliases.hpp"`
* `#include "tensor_index_tools.hpp"`
* `#include "vector_index_tools.hpp"`













## Namespaces

| Type | Name |
| ---: | :--- |
| namespace | [**ddcHelper**](namespaceddcHelper.md) <br> |


## Classes

| Type | Name |
| ---: | :--- |
| class | [**TensorCommon**](classTensorCommon.md) &lt;class DataStorageType, ValidIndexSet&gt;<br>_A superclass for_ [_**Tensor**_](classTensor.md) _calculations._[_**Tensor**_](classTensor.md) _classes containing data will inherit from this class. The class_[_**Tensor**_](classTensor.md) _will represent most Tensors but other subclasses may be necessary (e.g. to access a Vector in a_[_**VectorField**_](classVectorField.md) _)._ |






## Public Attributes

| Type | Name |
| ---: | :--- |
|  constexpr bool | [**is\_tensor\_type\_v**](#variable-is_tensor_type_v)   = `detail::enable\_tensor\_type&lt;std::remove\_const\_t&lt;std::remove\_reference\_t&lt;Type&gt;&gt;&gt;`<br> |
















## Public Functions

| Type | Name |
| ---: | :--- |
|  KOKKOS\_FUNCTION TensorType | [**operator\***](#function-operator) (Oelement\_type val, TensorType const & tensor) <br>_An operator to multiply all a tensor by a value._  |
|  KOKKOS\_FUNCTION TensorType | [**operator\***](#function-operator_1) (TensorType const & tensor, Oelement\_type val) <br>_An operator to multiply all the element of the current tensor by a value._  |
|  KOKKOS\_INLINE\_FUNCTION Coord&lt; Dims... &gt; | [**operator+**](#function-operator_2) (Coord&lt; Dims... &gt; const & coord, [**TensorCommon**](classTensorCommon.md)&lt; storage\_type, ddc::detail::TypeSeq&lt; Dims... &gt; &gt; const & tensor) <br> |
|  KOKKOS\_FUNCTION TensorType | [**operator+**](#function-operator_3) (TensorType const & tensor, Oelement\_type val) <br>_An operator to add two tensors elementwise._  |
|  KOKKOS\_INLINE\_FUNCTION Coord&lt; Dims... &gt; & | [**operator+=**](#function-operator_4) (Coord&lt; Dims... &gt; & coord, [**TensorCommon**](classTensorCommon.md)&lt; storage\_type, ddc::detail::TypeSeq&lt; Dims... &gt; &gt; const & tensor) <br> |
|  KOKKOS\_INLINE\_FUNCTION Coord&lt; Dims... &gt; | [**operator-**](#function-operator-) (Coord&lt; Dims... &gt; const & coord, [**TensorCommon**](classTensorCommon.md)&lt; storage\_type, ddc::detail::TypeSeq&lt; Dims... &gt; &gt; const & tensor) <br> |
|  KOKKOS\_FUNCTION TensorType | [**operator-**](#function-operator-_1) (TensorType const & tensor, TensorType const & val) <br>_An operator to subtract one tensor from another elementwise._  |
|  KOKKOS\_FUNCTION TensorType | [**operator-**](#function-operator-_2) (TensorType const & tensor) <br>_An operator to get the negation of a tensor elementwise._  |
|  KOKKOS\_INLINE\_FUNCTION Coord&lt; Dims... &gt; & | [**operator-=**](#function-operator-_3) (Coord&lt; Dims... &gt; & coord, [**TensorCommon**](classTensorCommon.md)&lt; storage\_type, ddc::detail::TypeSeq&lt; Dims... &gt; &gt; const & tensor) <br> |
|  KOKKOS\_FUNCTION TensorType | [**operator/**](#function-operator_5) (TensorType const & tensor, Oelement\_type val) <br>_An operator to divide all the element of the current tensor by a value._  |




























## Public Attributes Documentation




### variable is\_tensor\_type\_v 

```C++
constexpr bool is_tensor_type_v;
```




<hr>
## Public Functions Documentation




### function operator\* 

_An operator to multiply all a tensor by a value._ 
```C++
template<class TensorType, class Oelement_type, std::enable_if_t< is_tensor_type_v< TensorType >, bool >>
KOKKOS_FUNCTION TensorType operator* (
    Oelement_type val,
    TensorType const & tensor
) 
```





**Parameters:**


* `val` The value by which the elements should be multiplied. 
* `tensor` The tensor being multiplied. 



**Returns:**

A new tensor containing the result of the multiplication. 





        

<hr>



### function operator\* 

_An operator to multiply all the element of the current tensor by a value._ 
```C++
template<class TensorType, class Oelement_type, std::enable_if_t< is_tensor_type_v< TensorType >, bool >>
KOKKOS_FUNCTION TensorType operator* (
    TensorType const & tensor,
    Oelement_type val
) 
```





**Parameters:**


* `tensor` The tensor being multiplied. 
* `val` The value by which the elements should be multiplied. 



**Returns:**

A new tensor containing the result of the multiplication. 





        

<hr>



### function operator+ 

```C++
template<class storage_type, class... Dims>
KOKKOS_INLINE_FUNCTION Coord< Dims... > operator+ (
    Coord< Dims... > const & coord,
    TensorCommon < storage_type, ddc::detail::TypeSeq< Dims... > > const & tensor
) 
```



An operator to add the elements of a tensor to a coordinate. This can be useful in some calculations, e.g when calculating the foot of a characteristic. 

**Parameters:**


* `coord` The coordinate to which the tensor is added. 
* `tensor` The tensor to be added to the coordinate. 



**Returns:**

The new coordinate. 





        

<hr>



### function operator+ 

_An operator to add two tensors elementwise._ 
```C++
template<class TensorType, class Oelement_type, std::enable_if_t< is_tensor_type_v< TensorType >, bool >>
KOKKOS_FUNCTION TensorType operator+ (
    TensorType const & tensor,
    Oelement_type val
) 
```





**Parameters:**


* `tensor` The first tensor in the addition. 
* `val` The second tensor in the addition. 



**Returns:**

A new tensor containing the result of the addition. 





        

<hr>



### function operator+= 

```C++
template<class storage_type, class... Dims>
KOKKOS_INLINE_FUNCTION Coord< Dims... > & operator+= (
    Coord< Dims... > & coord,
    TensorCommon < storage_type, ddc::detail::TypeSeq< Dims... > > const & tensor
) 
```



An operator to add the elements of a tensor to a coordinate. This can be useful in some calculations, e.g when calculating the foot of a characteristic. 

**Parameters:**


* `coord` The coordinate to which the tensor is added. 
* `tensor` The tensor to be added to the coordinate. 



**Returns:**

The new coordinate. 





        

<hr>



### function operator- 

```C++
template<class storage_type, class... Dims>
KOKKOS_INLINE_FUNCTION Coord< Dims... > operator- (
    Coord< Dims... > const & coord,
    TensorCommon < storage_type, ddc::detail::TypeSeq< Dims... > > const & tensor
) 
```



An operator to add the elements of a tensor to a coordinate. This can be useful in some calculations, e.g when calculating the foot of a characteristic. 

**Parameters:**


* `coord` The coordinate from which the tensor is subtracted. 
* `tensor` The tensor to be subtracted from the coordinate. 



**Returns:**

The new coordinate. 





        

<hr>



### function operator- 

_An operator to subtract one tensor from another elementwise._ 
```C++
template<class TensorType, std::enable_if_t< is_tensor_type_v< TensorType >, bool >>
KOKKOS_FUNCTION TensorType operator- (
    TensorType const & tensor,
    TensorType const & val
) 
```





**Parameters:**


* `tensor` The tensor which is subtracted from. 
* `val` The tensor that should be subtracted. 



**Returns:**

A new tensor containing the result of the subtraction. 





        

<hr>



### function operator- 

_An operator to get the negation of a tensor elementwise._ 
```C++
template<class TensorType, std::enable_if_t< is_tensor_type_v< TensorType >, bool >>
KOKKOS_FUNCTION TensorType operator- (
    TensorType const & tensor
) 
```





**Parameters:**


* `tensor` The tensor to be negated. 



**Returns:**

A new tensor containing the result of the negation. 





        

<hr>



### function operator-= 

```C++
template<class storage_type, class... Dims>
KOKKOS_INLINE_FUNCTION Coord< Dims... > & operator-= (
    Coord< Dims... > & coord,
    TensorCommon < storage_type, ddc::detail::TypeSeq< Dims... > > const & tensor
) 
```



An operator to add the elements of a tensor to a coordinate. This can be useful in some calculations, e.g when calculating the foot of a characteristic. 

**Parameters:**


* `coord` The coordinate from which the tensor is subtracted. 
* `tensor` The tensor to be subtracted from the coordinate. 



**Returns:**

The new coordinate. 





        

<hr>



### function operator/ 

_An operator to divide all the element of the current tensor by a value._ 
```C++
template<class TensorType, class Oelement_type, std::enable_if_t< is_tensor_type_v< TensorType >, bool >>
KOKKOS_FUNCTION TensorType operator/ (
    TensorType const & tensor,
    Oelement_type val
) 
```





**Parameters:**


* `tensor` The tensor being divided. 
* `val` The value by which the elements should be multiplied. 



**Returns:**

A new tensor containing the result of the multiplication. 





        

<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/data_types/tensor_common.hpp`

