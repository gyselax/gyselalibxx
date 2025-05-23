

# Class TensorCommon

**template &lt;class DataStorageType, class... ValidIndexSet&gt;**



[**ClassList**](annotated.md) **>** [**TensorCommon**](classTensorCommon.md)



_A superclass for_ [_**Tensor**_](classTensor.md) _calculations._[_**Tensor**_](classTensor.md) _classes containing data will inherit from this class. The class_[_**Tensor**_](classTensor.md) _will represent most Tensors but other subclasses may be necessary (e.g. to access a Vector in a_[_**VectorField**_](classVectorField.md) _)._[More...](#detailed-description)

* `#include <tensor_common.hpp>`





Inherited by the following classes: [Tensor](classTensor.md)












## Public Types

| Type | Name |
| ---: | :--- |
| typedef typename DataStorageType::element\_type | [**element\_type**](#typedef-element_type)  <br>_The type of the elements of the tensor._  |
| typedef ddc::detail::TypeSeq&lt; ValidIndexSet... &gt; | [**index\_set**](#typedef-index_set)  <br>_The TensorIndexSet describing the possible indices._  |
| typedef ddc::type\_seq\_element\_t&lt; dim, [**index\_set**](classTensorCommon.md#typedef-index_set) &gt; | [**vector\_index\_set\_t**](#typedef-vector_index_set_t)  <br>_A helper type alias to get all possible indices along a specified dimension._  |




















## Public Functions

| Type | Name |
| ---: | :--- |
|  KOKKOS\_FUNCTION [**element\_type**](classTensorCommon.md#typedef-element_type) & | [**get**](#function-get-12) () <br>_Get a modifiable reference to an element of the tensor._  |
|  KOKKOS\_FUNCTION [**element\_type**](classTensorCommon.md#typedef-element_type) const & | [**get**](#function-get-22) () const<br>_Get an element of the tensor._  |
|  KOKKOS\_FUNCTION bool | [**operator!=**](#function-operator) ([**TensorCommon**](classTensorCommon.md) const & o\_tensor) const<br>_An operator to compare one tensor to another elementwise._  |
|  KOKKOS\_FUNCTION [**TensorCommon**](classTensorCommon.md) & | [**operator\*=**](#function-operator_1) (Oelement\_type val) <br>_An operator to multiply all the element of the current tensor by a value._  |
|  KOKKOS\_FUNCTION [**TensorCommon**](classTensorCommon.md) & | [**operator+=**](#function-operator_2) ([**TensorCommon**](classTensorCommon.md) const & val) <br>_An operator to add two tensors elementwise._  |
|  KOKKOS\_FUNCTION [**TensorCommon**](classTensorCommon.md) & | [**operator-=**](#function-operator-) ([**TensorCommon**](classTensorCommon.md) const & val) <br>_An operator to subtract one tensor from another elementwise._  |
|  KOKKOS\_FUNCTION [**TensorCommon**](classTensorCommon.md) & | [**operator/=**](#function-operator_3) (Oelement\_type val) <br>_An operator to divide all the element of the current tensor by a value._  |
|  KOKKOS\_FUNCTION [**TensorCommon**](classTensorCommon.md) & | [**operator=**](#function-operator_4) ([**TensorCommon**](classTensorCommon.md) const & other) <br>_A copy operator._  |
|  KOKKOS\_FUNCTION [**TensorCommon**](classTensorCommon.md) & | [**operator=**](#function-operator_5) ([**TensorCommon**](classTensorCommon.md) && other) <br>_A move assign operator._  |
|  KOKKOS\_FUNCTION bool | [**operator==**](#function-operator_6) ([**TensorCommon**](classTensorCommon.md) const & o\_tensor) const<br>_An operator to compare one tensor to another elementwise._  |


## Public Static Functions

| Type | Name |
| ---: | :--- |
|  KOKKOS\_FUNCTION constexpr std::size\_t | [**rank**](#function-rank) () <br>_The rank of the tensor. This is equivalent to the number of indices required to access an element of the tensor._  |
|  KOKKOS\_FUNCTION constexpr std::size\_t | [**size**](#function-size) () <br>_The size of the tensor. This is the number of elements in the tensor._  |






## Protected Attributes

| Type | Name |
| ---: | :--- |
|  DataStorageType | [**m\_data**](#variable-m_data)  <br> |


## Protected Static Attributes

| Type | Name |
| ---: | :--- |
|  constexpr std::size\_t | [**s\_n\_elements**](#variable-s_n_elements)   = `(ddc::type\_seq\_size\_v&lt;ValidIndexSet&gt; \* ...)`<br>_The number of elements in the mdspan._  |














## Protected Functions

| Type | Name |
| ---: | :--- |
|  KOKKOS\_DEFAULTED\_FUNCTION | [**TensorCommon**](#function-tensorcommon-13) () = default<br>_Construct an uninitialised tensor object._  |
|  KOKKOS\_DEFAULTED\_FUNCTION | [**TensorCommon**](#function-tensorcommon-23) ([**TensorCommon**](classTensorCommon.md) const & o\_tensor) = default<br>_Construct a tensor object by copying an existing tensor of exactly the same type. This method can be called implicitly._  |
|  KOKKOS\_DEFAULTED\_FUNCTION | [**TensorCommon**](#function-tensorcommon-33) ([**TensorCommon**](classTensorCommon.md) && o\_tensor) = default<br>_Move-construct a tensor object by copying an existing tensor of exactly the same type. This method can be called implicitly._  |




## Detailed Description




**Template parameters:**


* `element_type` The type of the elements of the tensor (usually double/complex). 
* `LayoutType` The way in which the underlying mdspan will be laid out in memory, usually Kokkos::layout\_right or Kokkos::layout\_stride. 
* `ValidIndexSet` The indices that can be used along each dimension of the tensor. 




    
## Public Types Documentation




### typedef element\_type 

_The type of the elements of the tensor._ 
```C++
using TensorCommon< DataStorageType, ValidIndexSet >::element_type =  typename DataStorageType::element_type;
```




<hr>



### typedef index\_set 

_The TensorIndexSet describing the possible indices._ 
```C++
using TensorCommon< DataStorageType, ValidIndexSet >::index_set =  ddc::detail::TypeSeq<ValidIndexSet...>;
```




<hr>



### typedef vector\_index\_set\_t 

_A helper type alias to get all possible indices along a specified dimension._ 
```C++
using TensorCommon< DataStorageType, ValidIndexSet >::vector_index_set_t =  ddc::type_seq_element_t<dim, index_set>;
```





**Template parameters:**


* `Dim` The dimension of interest (0 &lt;= dim &lt; [**rank()**](classTensorCommon.md#function-rank)). 




        

<hr>
## Public Functions Documentation




### function get [1/2]

_Get a modifiable reference to an element of the tensor._ 
```C++
template<class QueryTensorIndexElement>
inline KOKKOS_FUNCTION element_type & TensorCommon::get () 
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
inline KOKKOS_FUNCTION element_type const & TensorCommon::get () const
```





**Template parameters:**


* `QueryIndexTag` A type describing the relevant index. 



**Returns:**

The relevant element of the tensor. 





        

<hr>



### function operator!= 

_An operator to compare one tensor to another elementwise._ 
```C++
inline KOKKOS_FUNCTION bool TensorCommon::operator!= (
    TensorCommon const & o_tensor
) const
```





**Parameters:**


* `o_tensor` The tensor that should be compared with the current tensor. 



**Returns:**

False if the tensors are equal, true otherwise. 





        

<hr>



### function operator\*= 

_An operator to multiply all the element of the current tensor by a value._ 
```C++
template<class Oelement_type>
inline KOKKOS_FUNCTION TensorCommon & TensorCommon::operator*= (
    Oelement_type val
) 
```





**Parameters:**


* `val` The value by which the elements should be multiplied. 



**Returns:**

A reference to the current modified tensor. 





        

<hr>



### function operator+= 

_An operator to add two tensors elementwise._ 
```C++
inline KOKKOS_FUNCTION TensorCommon & TensorCommon::operator+= (
    TensorCommon const & val
) 
```





**Parameters:**


* `val` The tensor that should be added to the current tensor. 



**Returns:**

A reference to the current modified tensor. 





        

<hr>



### function operator-= 

_An operator to subtract one tensor from another elementwise._ 
```C++
inline KOKKOS_FUNCTION TensorCommon & TensorCommon::operator-= (
    TensorCommon const & val
) 
```





**Parameters:**


* `val` The tensor that should be subtracted from the current tensor. 



**Returns:**

A reference to the current modified tensor. 





        

<hr>



### function operator/= 

_An operator to divide all the element of the current tensor by a value._ 
```C++
template<class Oelement_type>
inline KOKKOS_FUNCTION TensorCommon & TensorCommon::operator/= (
    Oelement_type val
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
inline KOKKOS_FUNCTION TensorCommon & TensorCommon::operator= (
    TensorCommon const & other
) 
```





**Parameters:**


* `other` The tensor to be copied. 



**Returns:**

A reference to the current tensor. 





        

<hr>



### function operator= 

_A move assign operator._ 
```C++
inline KOKKOS_FUNCTION TensorCommon & TensorCommon::operator= (
    TensorCommon && other
) 
```





**Parameters:**


* `other` The tensor to be copied. 



**Returns:**

A reference to the current tensor. 





        

<hr>



### function operator== 

_An operator to compare one tensor to another elementwise._ 
```C++
inline KOKKOS_FUNCTION bool TensorCommon::operator== (
    TensorCommon const & o_tensor
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
static inline KOKKOS_FUNCTION constexpr std::size_t TensorCommon::rank () 
```





**Returns:**

The rank of the tensor. 





        

<hr>



### function size 

_The size of the tensor. This is the number of elements in the tensor._ 
```C++
static inline KOKKOS_FUNCTION constexpr std::size_t TensorCommon::size () 
```





**Returns:**

The number of elements in the tensor. 





        

<hr>
## Protected Attributes Documentation




### variable m\_data 

```C++
DataStorageType TensorCommon< DataStorageType, ValidIndexSet >::m_data;
```



The object providing access to the underlying mdspan via an operator(). The data in this object may be references or owned data depending on the subclass. 


        

<hr>
## Protected Static Attributes Documentation




### variable s\_n\_elements 

_The number of elements in the mdspan._ 
```C++
constexpr std::size_t TensorCommon< DataStorageType, ValidIndexSet >::s_n_elements;
```




<hr>
## Protected Functions Documentation




### function TensorCommon [1/3]

_Construct an uninitialised tensor object._ 
```C++
KOKKOS_DEFAULTED_FUNCTION TensorCommon::TensorCommon () = default
```




<hr>



### function TensorCommon [2/3]

_Construct a tensor object by copying an existing tensor of exactly the same type. This method can be called implicitly._ 
```C++
KOKKOS_DEFAULTED_FUNCTION TensorCommon::TensorCommon (
    TensorCommon const & o_tensor
) = default
```





**Parameters:**


* `o_tensor` The tensor to be copied. 




        

<hr>



### function TensorCommon [3/3]

_Move-construct a tensor object by copying an existing tensor of exactly the same type. This method can be called implicitly._ 
```C++
KOKKOS_DEFAULTED_FUNCTION TensorCommon::TensorCommon (
    TensorCommon && o_tensor
) = default
```





**Parameters:**


* `o_tensor` The tensor to be copied. 




        

<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/data_types/tensor_common.hpp`

