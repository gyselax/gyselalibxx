

# File tensor.hpp



[**FileList**](files.md) **>** [**data\_types**](dir_eaa769653453aaefd8cc10e98e9bb3eb.md) **>** [**tensor.hpp**](tensor_8hpp.md)

[Go to the source code of this file](tensor_8hpp_source.md)



* `#include <array>`
* `#include <cassert>`
* `#include <ddc/ddc.hpp>`
* `#include "ddc_aliases.hpp"`
* `#include "tensor_common.hpp"`
* `#include "tensor_index_tools.hpp"`
* `#include "vector_index_tools.hpp"`













## Namespaces

| Type | Name |
| ---: | :--- |
| namespace | [**Kokkos**](namespaceKokkos.md) <br> |


## Classes

| Type | Name |
| ---: | :--- |
| struct | [**reduction\_identity&lt; Tensor&lt; ElementType, ValidIndexSet... &gt; &gt;**](structKokkos_1_1reduction__identity_3_01Tensor_3_01ElementType_00_01ValidIndexSet_8_8_8_01_4_01_4.md) &lt;class ElementType, ValidIndexSet&gt;<br>_A specialisation of Kokkos::reduction\_identity to allow calling parallel\_transform on tensors._  |
| class | [**Tensor**](classTensor.md) &lt;class ElementType, class ValidIndexSetFirstDim, ValidIndexSet&gt;<br>_A class representing a_ [_**Tensor**_](classTensor.md) _._ |


## Public Types

| Type | Name |
| ---: | :--- |
| typedef [**Tensor**](classTensor.md)&lt; double, ValidIndexSet... &gt; | [**DTensor**](#typedef-dtensor)  <br>_A helper type alias to get a tensor containing doubles._  |
| typedef [**Vector**](classTensor.md)&lt; double, Dims... &gt; | [**DVector**](#typedef-dvector)  <br>_A helper type alias to get a 1D tensor (a vector) of doubles._  |
| typedef [**Tensor**](classTensor.md)&lt; ElementType, VectorIndexSet&lt; Dims... &gt; &gt; | [**Vector**](#typedef-vector)  <br>_A helper type alias to get a 1D tensor (a vector)._  |
| typedef typename detail::ToTensor&lt; ElementType, TypeSeqValidIndexSet &gt;::type | [**to\_tensor\_t**](#typedef-to_tensor_t)  <br> |




















## Public Functions

| Type | Name |
| ---: | :--- |
|  std::ostream & | [**operator&lt;&lt;**](#function-operator) (std::ostream & o, [**Vector**](classTensor.md)&lt; ElementType, Dim1, TailDims... &gt; vec) <br> |




























## Public Types Documentation




### typedef DTensor 

_A helper type alias to get a tensor containing doubles._ 
```C++
using DTensor =  Tensor<double, ValidIndexSet...>;
```





**Template parameters:**


* `ValidIndexSet` The indices that can be used along each dimension of the tensor. 




        

<hr>



### typedef DVector 

_A helper type alias to get a 1D tensor (a vector) of doubles._ 
```C++
using DVector =  Vector<double, Dims...>;
```





**Template parameters:**


* `Dims` The dimensions that can be used to index the vector. 




        

<hr>



### typedef Vector 

_A helper type alias to get a 1D tensor (a vector)._ 
```C++
using Vector =  Tensor<ElementType, VectorIndexSet<Dims...> >;
```





**Template parameters:**


* `ElementType` The type of the elements of the tensor (usually double/complex). 
* `Dims` The dimensions that can be used to index the vector. 




        

<hr>



### typedef to\_tensor\_t 

```C++
using to_tensor_t =  typename detail::ToTensor<ElementType, TypeSeqValidIndexSet>::type;
```




<hr>
## Public Functions Documentation




### function operator&lt;&lt; 

```C++
template<class ElementType, class Dim1, class... TailDims>
std::ostream & operator<< (
    std::ostream & o,
    Vector < ElementType, Dim1, TailDims... > vec
) 
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/data_types/tensor.hpp`

