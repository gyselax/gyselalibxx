

# Class tensor\_tools::TensorIndexElement&lt; ddc::TypeSeq&lt; ValidatingVectorIndexSets... &gt;, Dims... &gt;

**template &lt;class... ValidatingVectorIndexSets, class... Dims&gt;**



[**ClassList**](annotated.md) **>** [**tensor\_tools**](namespacetensor__tools.md) **>** [**TensorIndexElement&lt; ddc::TypeSeq&lt; ValidatingVectorIndexSets... &gt;, Dims... &gt;**](classtensor__tools_1_1TensorIndexElement_3_01ddc_1_1TypeSeq_3_01ValidatingVectorIndexSets_8_8_8_01_4_00_01Dims_8_8_8_01_4.md)



_A class describing an index of a tensor. See the general_ [_**TensorIndexElement**_](classtensor__tools_1_1TensorIndexElement.md) _template for more details._[More...](#detailed-description)

* `#include <tensor_index_tools.hpp>`

















## Public Types

| Type | Name |
| ---: | :--- |
| typedef ddc::TypeSeq&lt; Dims... &gt; | [**IdxTypeSeq**](#typedef-idxtypeseq)  <br>_The type sequence containing all the indices applied to the tensor._  |
| typedef ddc::type\_seq\_element\_t&lt; IDim, [**IdxTypeSeq**](classtensor__tools_1_1TensorIndexElement_3_01ddc_1_1TypeSeq_3_01ValidatingVectorIndexSets_8_8_8_01_4_00_01Dims_8_8_8_01_4.md#typedef-idxtypeseq) &gt; | [**index\_on\_dim\_t**](#typedef-index_on_dim_t)  <br>_Get the types that describes the index in dimension IDim._  |
| typedef ValidatingTensorIndexSet | [**tensor\_index\_set**](#typedef-tensor_index_set)  <br>_Get the TensorIndexSet that this element indexes._  |






















## Public Static Functions

| Type | Name |
| ---: | :--- |
|  KOKKOS\_FUNCTION constexpr std::size\_t | [**index**](#function-index) () <br>_Get an integer which describes the position of this element in a 1D array which stores a tensor described by the TensorIndexSet._  |
|  KOKKOS\_FUNCTION constexpr std::size\_t | [**rank**](#function-rank) () <br>_The rank of the tensor. This is equivalent to the number of indices required to access an element of the tensor._  |


























## Detailed Description




**Template parameters:**


* `ValidatingVectorIndexSets` The VectorIndexSets in which each tag comprising this element (along each dimension) can be found. 
* `Dims` The dimensions at which the tensor is accessed. 




    
## Public Types Documentation




### typedef IdxTypeSeq 

_The type sequence containing all the indices applied to the tensor._ 
```C++
using tensor_tools::TensorIndexElement< ddc::TypeSeq< ValidatingVectorIndexSets... >, Dims... >::IdxTypeSeq =  ddc::TypeSeq<Dims...>;
```




<hr>



### typedef index\_on\_dim\_t 

_Get the types that describes the index in dimension IDim._ 
```C++
using tensor_tools::TensorIndexElement< ddc::TypeSeq< ValidatingVectorIndexSets... >, Dims... >::index_on_dim_t =  ddc::type_seq_element_t<IDim, IdxTypeSeq>;
```





**Template parameters:**


* `IDim` The index of interest. 




        

<hr>



### typedef tensor\_index\_set 

_Get the TensorIndexSet that this element indexes._ 
```C++
using tensor_tools::TensorIndexElement< ddc::TypeSeq< ValidatingVectorIndexSets... >, Dims... >::tensor_index_set =  ValidatingTensorIndexSet;
```




<hr>
## Public Static Functions Documentation




### function index 

_Get an integer which describes the position of this element in a 1D array which stores a tensor described by the TensorIndexSet._ 
```C++
static inline KOKKOS_FUNCTION constexpr std::size_t tensor_tools::TensorIndexElement< ddc::TypeSeq< ValidatingVectorIndexSets... >, Dims... >::index () 
```





**Returns:**

The index of this element in a 1D array. 





        

<hr>



### function rank 

_The rank of the tensor. This is equivalent to the number of indices required to access an element of the tensor._ 
```C++
static inline KOKKOS_FUNCTION constexpr std::size_t tensor_tools::TensorIndexElement< ddc::TypeSeq< ValidatingVectorIndexSets... >, Dims... >::rank () 
```





**Returns:**

The rank of the tensor. 





        

<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/data_types/tensor_index_tools.hpp`

