

# File vector\_index\_tools.hpp



[**FileList**](files.md) **>** [**data\_types**](dir_eaa769653453aaefd8cc10e98e9bb3eb.md) **>** [**vector\_index\_tools.hpp**](vector__index__tools_8hpp.md)

[Go to the source code of this file](vector__index__tools_8hpp_source.md)



* `#include <ddc/ddc.hpp>`













## Namespaces

| Type | Name |
| ---: | :--- |
| namespace | [**tensor\_tools**](namespacetensor__tools.md) <br> |


## Classes

| Type | Name |
| ---: | :--- |
| struct | [**GetContravariantDims&lt; VectorIndexSet&lt; Dims... &gt; &gt;**](structtensor__tools_1_1GetContravariantDims_3_01VectorIndexSet_3_01Dims_8_8_8_01_4_01_4.md) &lt;Dims&gt;<br>_A class to get a VectorIndexSet containing only contravariant dimensions._  |
| struct | [**GetCovariantDims&lt; VectorIndexSet&lt; Dims... &gt; &gt;**](structtensor__tools_1_1GetCovariantDims_3_01VectorIndexSet_3_01Dims_8_8_8_01_4_01_4.md) &lt;Dims&gt;<br>_A class to get a VectorIndexSet containing only covariant dimensions._  |
| struct | [**VectorIndexIdMap**](structtensor__tools_1_1VectorIndexIdMap.md) &lt;Id, class AssociatedVectorIndexSet&gt;<br>_A class representing a vector index identifier._  |
| struct | [**is\_contravariant\_vector\_index\_set&lt; VectorIndexSet&lt; Dims... &gt; &gt;**](structtensor__tools_1_1is__contravariant__vector__index__set_3_01VectorIndexSet_3_01Dims_8_8_8_01_4_01_4.md) &lt;Dims&gt;<br>_A helper structure to check if all the dimensions in a VectorIndexSet can represent contravariant indices._  |
| struct | [**is\_covariant\_vector\_index\_set&lt; VectorIndexSet&lt; Dims... &gt; &gt;**](structtensor__tools_1_1is__covariant__vector__index__set_3_01VectorIndexSet_3_01Dims_8_8_8_01_4_01_4.md) &lt;Dims&gt;<br>_A helper structure to check if all the dimensions in a VectorIndexSet can represent covariant indices._  |
| struct | [**is\_vector\_index\_set**](structtensor__tools_1_1is__vector__index__set.md) &lt;class Type&gt;<br>_A helper structure to recognise a VectorIndexSet type._  |
| struct | [**is\_vector\_index\_set&lt; VectorIndexSet&lt; Dims... &gt; &gt;**](structtensor__tools_1_1is__vector__index__set_3_01VectorIndexSet_3_01Dims_8_8_8_01_4_01_4.md) &lt;Dims&gt;<br> |
| struct | [**vector\_index\_set\_dual&lt; VectorIndexSet&lt; Dims... &gt; &gt;**](structtensor__tools_1_1vector__index__set__dual_3_01VectorIndexSet_3_01Dims_8_8_8_01_4_01_4.md) &lt;Dims&gt;<br>_The implementation of_ [_**vector\_index\_set\_dual**_](structtensor__tools_1_1vector__index__set__dual.md) _for a VectorIndexSet._ |


## Public Types

| Type | Name |
| ---: | :--- |
| typedef ddc::detail::TypeSeq&lt; Dims... &gt; | [**VectorIndexSet**](#typedef-vectorindexset)  <br>_A type alias to describe a set of dimensions that can be used to index a vector (e.g. x to get E\_x)._  |
| typedef typename [**tensor\_tools::GetContravariantDims**](structtensor__tools_1_1GetContravariantDims.md)&lt; AnyVectorIndexSet &gt;::type | [**get\_contravariant\_dims\_t**](#typedef-get_contravariant_dims_t)  <br> |
| typedef typename [**tensor\_tools::GetCovariantDims**](structtensor__tools_1_1GetCovariantDims.md)&lt; AnyVectorIndexSet &gt;::type | [**get\_covariant\_dims\_t**](#typedef-get_covariant_dims_t)  <br> |
| typedef typename [**tensor\_tools::vector\_index\_set\_dual**](structtensor__tools_1_1vector__index__set__dual.md)&lt; VectorIndexSet &gt;::type | [**vector\_index\_set\_dual\_t**](#typedef-vector_index_set_dual_t)  <br>_A type alias to find a VectorIndexSet describing the covariant indices from a VectorIndexSet describing contravariant indices or to find a VectorIndexSet describing the contravariant indices from a VectorIndexSet describing covariant indices._  |






## Public Static Attributes

| Type | Name |
| ---: | :--- |
|  constexpr bool | [**is\_contravariant\_vector\_index\_set\_v**](#variable-is_contravariant_vector_index_set_v)   = `[**tensor\_tools::is\_contravariant\_vector\_index\_set**](structtensor__tools_1_1is__contravariant__vector__index__set.md)&lt;VectorIndexSet&gt;::value`<br>_A compile-time boolean to check if all the dimensions in a VectorIndexSet can represent contravariant indices._  |
|  constexpr bool | [**is\_covariant\_vector\_index\_set\_v**](#variable-is_covariant_vector_index_set_v)   = `[**tensor\_tools::is\_covariant\_vector\_index\_set**](structtensor__tools_1_1is__covariant__vector__index__set.md)&lt;VectorIndexSet&gt;::value`<br>_A compile-time boolean to check if all the dimensions in a VectorIndexSet can represent covariant indices._  |
|  constexpr bool | [**is\_vector\_index\_set\_v**](#variable-is_vector_index_set_v)   = `[**tensor\_tools::is\_vector\_index\_set**](structtensor__tools_1_1is__vector__index__set.md)&lt;Type&gt;::value`<br> |










































## Public Types Documentation




### typedef VectorIndexSet 

_A type alias to describe a set of dimensions that can be used to index a vector (e.g. x to get E\_x)._ 
```C++
using VectorIndexSet =  ddc::detail::TypeSeq<Dims...>;
```




<hr>



### typedef get\_contravariant\_dims\_t 

```C++
using get_contravariant_dims_t =  typename tensor_tools::GetContravariantDims<AnyVectorIndexSet>::type;
```




<hr>



### typedef get\_covariant\_dims\_t 

```C++
using get_covariant_dims_t =  typename tensor_tools::GetCovariantDims<AnyVectorIndexSet>::type;
```




<hr>



### typedef vector\_index\_set\_dual\_t 

_A type alias to find a VectorIndexSet describing the covariant indices from a VectorIndexSet describing contravariant indices or to find a VectorIndexSet describing the contravariant indices from a VectorIndexSet describing covariant indices._ 
```C++
using vector_index_set_dual_t =  typename tensor_tools::vector_index_set_dual<VectorIndexSet>::type;
```




<hr>
## Public Static Attributes Documentation




### variable is\_contravariant\_vector\_index\_set\_v 

_A compile-time boolean to check if all the dimensions in a VectorIndexSet can represent contravariant indices._ 
```C++
constexpr bool is_contravariant_vector_index_set_v;
```




<hr>



### variable is\_covariant\_vector\_index\_set\_v 

_A compile-time boolean to check if all the dimensions in a VectorIndexSet can represent covariant indices._ 
```C++
constexpr bool is_covariant_vector_index_set_v;
```




<hr>



### variable is\_vector\_index\_set\_v 

```C++
constexpr bool is_vector_index_set_v;
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/data_types/vector_index_tools.hpp`

