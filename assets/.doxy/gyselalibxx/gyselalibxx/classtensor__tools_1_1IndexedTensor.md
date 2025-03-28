

# Class tensor\_tools::IndexedTensor

**template &lt;class TensorType, class TypeSeqVectorIndexIdMap&gt;**



[**ClassList**](annotated.md) **>** [**tensor\_tools**](namespacetensor__tools.md) **>** [**IndexedTensor**](classtensor__tools_1_1IndexedTensor.md)



_A class to capture the description of a tensor indexed at a specific component. This class should not be explicitly declared in user code. It is the output of a call to the index&lt;...&gt; function and is an input to the tensor\_mul function._ [More...](#detailed-description)

* `#include <indexed_tensor.hpp>`

















## Public Types

| Type | Name |
| ---: | :--- |
| typedef TypeSeqVectorIndexIdMap | [**index\_pattern**](#typedef-index_pattern)  <br>_The TensorIndexIdMap describing how the tensor is indexed._  |
| typedef TensorType | [**tensor\_type**](#typedef-tensor_type)  <br>_The type of the tensor being indexed._  |




















## Public Functions

| Type | Name |
| ---: | :--- |
|  KOKKOS\_FUNCTION | [**IndexedTensor**](#function-indexedtensor-13) (TensorType & tensor) <br>_A constructor of an indexed tensor._  |
|   | [**IndexedTensor**](#function-indexedtensor-23) ([**IndexedTensor**](classtensor__tools_1_1IndexedTensor.md) const &) = delete<br> |
|   | [**IndexedTensor**](#function-indexedtensor-33) ([**IndexedTensor**](classtensor__tools_1_1IndexedTensor.md) &&) = delete<br> |
|  KOKKOS\_FUNCTION TensorType & | [**operator()**](#function-operator) () <br>_An operator to access the underlying tensor._  |
|  KOKKOS\_FUNCTION TensorType const & | [**operator()**](#function-operator_1) () const<br>_An operator to access the underlying tensor._  |
|  [**IndexedTensor**](classtensor__tools_1_1IndexedTensor.md) & | [**operator=**](#function-operator_2) ([**IndexedTensor**](classtensor__tools_1_1IndexedTensor.md) const &) = delete<br> |
|  [**IndexedTensor**](classtensor__tools_1_1IndexedTensor.md) & | [**operator=**](#function-operator_3) ([**IndexedTensor**](classtensor__tools_1_1IndexedTensor.md) &&) = delete<br> |
|  KOKKOS\_DEFAULTED\_FUNCTION | [**~IndexedTensor**](#function-indexedtensor) () noexcept<br> |




























## Detailed Description




**Template parameters:**


* `TensorType` The tensor being indexed. 
* `TypeSeqVectorIndexIdMap` A TypeSeq containing VectorIndexIdMaps describing the indices used to access the tensor. 




    
## Public Types Documentation




### typedef index\_pattern 

_The TensorIndexIdMap describing how the tensor is indexed._ 
```C++
using tensor_tools::IndexedTensor< TensorType, TypeSeqVectorIndexIdMap >::index_pattern =  TypeSeqVectorIndexIdMap;
```




<hr>



### typedef tensor\_type 

_The type of the tensor being indexed._ 
```C++
using tensor_tools::IndexedTensor< TensorType, TypeSeqVectorIndexIdMap >::tensor_type =  TensorType;
```




<hr>
## Public Functions Documentation




### function IndexedTensor [1/3]

_A constructor of an indexed tensor._ 
```C++
inline explicit KOKKOS_FUNCTION tensor_tools::IndexedTensor::IndexedTensor (
    TensorType & tensor
) 
```





**Parameters:**


* `tensor` The tensor to be indexed. 




        

<hr>



### function IndexedTensor [2/3]

```C++
tensor_tools::IndexedTensor::IndexedTensor (
    IndexedTensor const &
) = delete
```




<hr>



### function IndexedTensor [3/3]

```C++
tensor_tools::IndexedTensor::IndexedTensor (
    IndexedTensor &&
) = delete
```




<hr>



### function operator() 

_An operator to access the underlying tensor._ 
```C++
inline KOKKOS_FUNCTION TensorType & tensor_tools::IndexedTensor::operator() () 
```





**Returns:**

The underlying tensor. 





        

<hr>



### function operator() 

_An operator to access the underlying tensor._ 
```C++
inline KOKKOS_FUNCTION TensorType const & tensor_tools::IndexedTensor::operator() () const
```





**Returns:**

The underlying tensor. 





        

<hr>



### function operator= 

```C++
IndexedTensor & tensor_tools::IndexedTensor::operator= (
    IndexedTensor const &
) = delete
```




<hr>



### function operator= 

```C++
IndexedTensor & tensor_tools::IndexedTensor::operator= (
    IndexedTensor &&
) = delete
```




<hr>



### function ~IndexedTensor 

```C++
KOKKOS_DEFAULTED_FUNCTION tensor_tools::IndexedTensor::~IndexedTensor () noexcept
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/data_types/indexed_tensor.hpp`

