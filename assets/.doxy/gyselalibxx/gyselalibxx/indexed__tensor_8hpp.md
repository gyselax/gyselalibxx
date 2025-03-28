

# File indexed\_tensor.hpp



[**FileList**](files.md) **>** [**data\_types**](dir_eaa769653453aaefd8cc10e98e9bb3eb.md) **>** [**indexed\_tensor.hpp**](indexed__tensor_8hpp.md)

[Go to the source code of this file](indexed__tensor_8hpp_source.md)

[More...](#detailed-description)

* `#include "tensor.hpp"`
* `#include "tensor_index_tools.hpp"`
* `#include "vector_index_tools.hpp"`













## Namespaces

| Type | Name |
| ---: | :--- |
| namespace | [**tensor\_tools**](namespacetensor__tools.md) <br> |


## Classes

| Type | Name |
| ---: | :--- |
| class | [**IndexedTensor**](classtensor__tools_1_1IndexedTensor.md) &lt;class TensorType, class TypeSeqVectorIndexIdMap&gt;<br>_A class to capture the description of a tensor indexed at a specific component. This class should not be explicitly declared in user code. It is the output of a call to the index&lt;...&gt; function and is an input to the tensor\_mul function._  |






















## Public Functions

| Type | Name |
| ---: | :--- |
|  KOKKOS\_FUNCTION auto | [**index**](#function-index) (TensorType const & tensor) <br>_A tool to build an IndexedTensor object from a set of indices. This object can be used to carry out a tensor multiplication using the tensor\_mul function. Example:_  |
|  KOKKOS\_FUNCTION auto | [**tensor\_mul**](#function-tensor_mul) (IndexedTensorType... tensor\_to\_mul) <br>_A function to multiply tensors together. IndexedTensors are used to describe the index pattern of the multiplication. Example:_  |




























## Detailed Description


File providing functions to carry out tensor calculus operations via index patterns. 


    
## Public Functions Documentation




### function index 

_A tool to build an IndexedTensor object from a set of indices. This object can be used to carry out a tensor multiplication using the tensor\_mul function. Example:_ 
```C++
template<char... ids, class TensorType>
KOKKOS_FUNCTION auto index (
    TensorType const & tensor
) 
```




```C++
index<'i', 'j', 'k'>(my_3d_tensor)
```
 

**Parameters:**


* `tensor` The tensor that is being indexed. 



**Template parameters:**


* `ids` The characters describing the index pattern used to index the tensor. 



**Returns:**

An IndexedTensor object. 





        

<hr>



### function tensor\_mul 

_A function to multiply tensors together. IndexedTensors are used to describe the index pattern of the multiplication. Example:_ 
```C++
template<class... IndexedTensorType>
KOKKOS_FUNCTION auto tensor_mul (
    IndexedTensorType... tensor_to_mul
) 
```




```C++
// Create a tensor C such that C_{ik} = A_{ij} B^{j}_{k}
Tensor C = tensor_mul(index<'i', 'j'>(A), index<'j', 'k'>(B));
```
 

**Parameters:**


* `tensor_to_mul` The indexed tensor objects which should be multiplied together. 



**Returns:**

A tensor that is the result of multiplying the objects according to the index pattern 





        

<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/data_types/indexed_tensor.hpp`

