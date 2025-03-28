

# Namespace tensor\_tools



[**Namespace List**](namespaces.md) **>** [**tensor\_tools**](namespacetensor__tools.md)



[More...](#detailed-description)
















## Classes

| Type | Name |
| ---: | :--- |
| struct | [**GetContravariantDims**](structtensor__tools_1_1GetContravariantDims.md) &lt;class AnyVectorIndexSet&gt;<br>_A class to get a VectorIndexSet containing only contravariant dimensions._  |
| struct | [**GetContravariantDims&lt; VectorIndexSet&lt; Dims... &gt; &gt;**](structtensor__tools_1_1GetContravariantDims_3_01VectorIndexSet_3_01Dims_8_8_8_01_4_01_4.md) &lt;Dims&gt;<br>_A class to get a VectorIndexSet containing only contravariant dimensions._  |
| struct | [**GetCovariantDims**](structtensor__tools_1_1GetCovariantDims.md) &lt;class AnyVectorIndexSet&gt;<br>_A class to get a VectorIndexSet containing only covariant dimensions._  |
| struct | [**GetCovariantDims&lt; VectorIndexSet&lt; Dims... &gt; &gt;**](structtensor__tools_1_1GetCovariantDims_3_01VectorIndexSet_3_01Dims_8_8_8_01_4_01_4.md) &lt;Dims&gt;<br>_A class to get a VectorIndexSet containing only covariant dimensions._  |
| class | [**IndexedTensor**](classtensor__tools_1_1IndexedTensor.md) &lt;class TensorType, class TypeSeqVectorIndexIdMap&gt;<br>_A class to capture the description of a tensor indexed at a specific component. This class should not be explicitly declared in user code. It is the output of a call to the index&lt;...&gt; function and is an input to the tensor\_mul function._  |
| class | [**TensorIndexElement**](classtensor__tools_1_1TensorIndexElement.md) &lt;class ValidatingTensorIndexSet, Dims&gt;<br>_A class describing an index of a tensor. For example for a 2x2 metric tensor on an (x,y) plane the element_  _would have the index TensorIndexElement&lt;TensorIndexSetXY, X, X&gt;._ |
| struct | [**VectorIndexIdMap**](structtensor__tools_1_1VectorIndexIdMap.md) &lt;Id, class AssociatedVectorIndexSet&gt;<br>_A class representing a vector index identifier._  |
| struct | [**is\_contravariant\_vector\_index\_set**](structtensor__tools_1_1is__contravariant__vector__index__set.md) &lt;class VectorIndexSet&gt;<br> |
| struct | [**is\_contravariant\_vector\_index\_set&lt; VectorIndexSet&lt; Dims... &gt; &gt;**](structtensor__tools_1_1is__contravariant__vector__index__set_3_01VectorIndexSet_3_01Dims_8_8_8_01_4_01_4.md) &lt;Dims&gt;<br>_A helper structure to check if all the dimensions in a VectorIndexSet can represent contravariant indices._  |
| struct | [**is\_covariant\_vector\_index\_set**](structtensor__tools_1_1is__covariant__vector__index__set.md) &lt;class VectorIndexSet&gt;<br> |
| struct | [**is\_covariant\_vector\_index\_set&lt; VectorIndexSet&lt; Dims... &gt; &gt;**](structtensor__tools_1_1is__covariant__vector__index__set_3_01VectorIndexSet_3_01Dims_8_8_8_01_4_01_4.md) &lt;Dims&gt;<br>_A helper structure to check if all the dimensions in a VectorIndexSet can represent covariant indices._  |
| struct | [**is\_vector\_index\_set**](structtensor__tools_1_1is__vector__index__set.md) &lt;class Type&gt;<br>_A helper structure to recognise a VectorIndexSet type._  |
| struct | [**is\_vector\_index\_set&lt; VectorIndexSet&lt; Dims... &gt; &gt;**](structtensor__tools_1_1is__vector__index__set_3_01VectorIndexSet_3_01Dims_8_8_8_01_4_01_4.md) &lt;Dims&gt;<br> |
| struct | [**vector\_index\_set\_dual**](structtensor__tools_1_1vector__index__set__dual.md) &lt;class VectorIndexSet&gt;<br>_A helper structure to find a VectorIndexSet describing the covariant indices from a VectorIndexSet describing contravariant indices or to find a VectorIndexSet describing the contravariant indices from a VectorIndexSet describing covariant indices._  |
| struct | [**vector\_index\_set\_dual&lt; VectorIndexSet&lt; Dims... &gt; &gt;**](structtensor__tools_1_1vector__index__set__dual_3_01VectorIndexSet_3_01Dims_8_8_8_01_4_01_4.md) &lt;Dims&gt;<br>_The implementation of_ [_**vector\_index\_set\_dual**_](structtensor__tools_1_1vector__index__set__dual.md) _for a VectorIndexSet._ |


## Public Types

| Type | Name |
| ---: | :--- |
| typedef typename details::ExtractSubTensorElement&lt; TypeSeqVectorIndexIdMapGlobal, TypeSeqVectorIndexIdMapLocal, GlobalTensorIndexElement &gt;::type | [**extract\_sub\_tensor\_element\_t**](#typedef-extract_sub_tensor_element_t)  <br>_Extract the relevant elements of a_ [_**TensorIndexElement**_](classtensor__tools_1_1TensorIndexElement.md) _to create a sub-TensorIndexElement using a global and a local TypeSeq of VectorIndexIdMaps to identify the relevant elements. For example: for GlobalTensorIndexElement = TensorIndexElement&lt;X,Y&gt; with TypeSeqVectorIndexIdMapGlobal = TypeSeq&lt;_[_**VectorIndexIdMap**_](structtensor__tools_1_1VectorIndexIdMap.md) _&lt;'i', VectorIndexSet&lt;X, Y&gt;&gt;,_[_**VectorIndexIdMap**_](structtensor__tools_1_1VectorIndexIdMap.md) _&lt;'j', VectorIndexSet&lt;X, Y&gt;&gt;&gt; and TypeSeqVectorIndexIdMapLocal = TypeSeq&lt;_[_**VectorIndexIdMap**_](structtensor__tools_1_1VectorIndexIdMap.md) _&lt;'j', VectorIndexSet&lt;X, Y&gt;&gt;&gt; we obtain TensorIndexElement&lt;Y&gt;_ |
| typedef [**to\_tensor\_index\_element\_t**](namespacetensor__tools.md#typedef-to_tensor_index_element_t)&lt; [**get\_type\_seq\_vector\_index\_set\_t**](namespacetensor__tools.md#typedef-get_type_seq_vector_index_set_t)&lt; TypeSeqVectorIndexIdMap &gt;, typename details::GetNthTensorIndexElementFromMap&lt; TypeSeqVectorIndexIdMap, [**get\_nth\_tensor\_index\_element\_t**](namespacetensor__tools.md#typedef-get_nth_tensor_index_element_t)&lt; Elem, [**get\_type\_seq\_vector\_index\_set\_t**](namespacetensor__tools.md#typedef-get_type_seq_vector_index_set_t)&lt; [**unique\_indices\_t**](namespacetensor__tools.md#typedef-unique_indices_t)&lt; TypeSeqVectorIndexIdMap &gt; &gt; &gt;, ddc::type\_seq\_size\_v&lt; TypeSeqVectorIndexIdMap &gt; &gt;::type &gt; | [**get\_nth\_tensor\_index\_element\_from\_map\_t**](#typedef-get_nth_tensor_index_element_from_map_t)  <br>_Get the n-th valid index for a tensor which is accessed according to the pattern described by a TypeSeq of VectorIndexIdMaps. E.g. for a 2D tensor with components A\_{xx}, A\_{xy}, A\_{yx}, A\_{yy}, indexed with._  |
| typedef typename details::GetNthTensorIndexElement&lt; IndexPosition, ddc::type\_seq\_size\_v&lt; TypeSeqVectorIndexSet &gt;, TypeSeqVectorIndexSet &gt;::type | [**get\_nth\_tensor\_index\_element\_t**](#typedef-get_nth_tensor_index_element_t)  <br>_Get the_ [_**TensorIndexElement**_](classtensor__tools_1_1TensorIndexElement.md) _which indexes a_[_**Tensor**_](classTensor.md) _at the n-th position of its internal array. E.g. for a 2x2_[_**Tensor**_](classTensor.md) _, get\_nth\_tensor\_index\_element\_t&lt;1, TypeSeqVectorIndexSet&gt; returns the_[_**TensorIndexElement**_](classtensor__tools_1_1TensorIndexElement.md) _which indexes element 1 of the array, so the element {0,1} of the tensor._ |
| typedef typename details::ExtractTypeSeqIndexSet&lt; TypeSeqVectorIndexIdMap &gt;::type | [**get\_type\_seq\_vector\_index\_set\_t**](#typedef-get_type_seq_vector_index_set_t)  <br>_Get a TypeSeq of valid VectorIndexSets from a TypeSeq of VectorIndexIdMaps._  |
| typedef typename details::GetIndexIds&lt; TypeSeqVectorIndexIdMap &gt;::type | [**index\_identifiers\_t**](#typedef-index_identifiers_t)  <br>_Create a TypeSeq of integral\_constants containing all the character ids found in the TypeSeqVectorIndexIdMap._  |
| typedef typename details::GetNonRepeatedIndices&lt; TypeSeqVectorIndexIdMap, ddc::type\_seq\_size\_v&lt; TypeSeqVectorIndexIdMap &gt; - 1, [**index\_identifiers\_t**](namespacetensor__tools.md#typedef-index_identifiers_t)&lt; TypeSeqVectorIndexIdMap &gt; &gt;::type | [**non\_repeated\_indices\_t**](#typedef-non_repeated_indices_t)  <br>_Extract the VectorIndexIdMaps whose character id only appears once in a TypeSeq of VectorIndexIdMaps._  |
| typedef typename details::GetRelevantVectorIndexSets&lt; ID, TypeSeqVectorIndexIdMap &gt;::type | [**relevant\_vector\_index\_sets\_t**](#typedef-relevant_vector_index_sets_t)  <br>_Extract all VectorIndexSets which are described by the same character identifier._  |
| typedef typename details::ToTensorIndexElement&lt; ValidatingTensorIndexSet, TypeSeqTensorIndexTag &gt;::type | [**to\_tensor\_index\_element\_t**](#typedef-to_tensor_index_element_t)  <br>_Get a_ [_**TensorIndexElement**_](classtensor__tools_1_1TensorIndexElement.md) _from a TypeSeq of valid VectorIndexSets and a TypeSeq of indices._ |
| typedef typename details::GetUniqueIndices&lt; TypeSeqVectorIndexIdMap &gt;::type | [**unique\_indices\_t**](#typedef-unique_indices_t)  <br>_Create a TypeSeq of_ [_**VectorIndexIdMap**_](structtensor__tools_1_1VectorIndexIdMap.md) _in which each character id only appears once from a TypeSeq of VectorIndexIdMaps with repeat character ids._ |




## Public Attributes

| Type | Name |
| ---: | :--- |
|  constexpr std::size\_t | [**char\_occurrences\_v**](#variable-char_occurrences_v)   = `details::CountChar&lt;search\_char, CharTypeSeq&gt;::value`<br>_Count the number of instances of a character in a TypeSeq of literals._  |
|  constexpr bool | [**enable\_indexed\_tensor**](#variable-enable_indexed_tensor)   = `false`<br>_A boolean, true if the type is an_ [_**IndexedTensor**_](classtensor__tools_1_1IndexedTensor.md) _, false otherwise._ |
|  constexpr bool | [**enable\_indexed\_tensor&lt; IndexedTensor&lt; TensorType, TypeSeqVectorIndexIdMap &gt; &gt;**](#variable-enable_indexed_tensor-indexedtensor-tensortype-typeseqvectorindexidmap)   = `true`<br>_A boolean, true if the type is an_ [_**IndexedTensor**_](classtensor__tools_1_1IndexedTensor.md) _._ |
|  constexpr bool | [**is\_indexed\_tensor\_v**](#variable-is_indexed_tensor_v)   = `[**enable\_indexed\_tensor**](namespacetensor__tools.md#variable-enable_indexed_tensor)&lt;std::remove\_const\_t&lt;std::remove\_reference\_t&lt;Type&gt;&gt;&gt;`<br>_A tool to check if a type is an_ [_**IndexedTensor**_](classtensor__tools_1_1IndexedTensor.md) _._ |
|  constexpr bool | [**is\_tensor\_index\_element\_v**](#variable-is_tensor_index_element_v)   = `details::enable\_tensor\_index\_element&lt;std::remove\_const\_t&lt;std::remove\_reference\_t&lt;Type&gt;&gt;&gt;`<br>_A helper to check if a type is a_ [_**TensorIndexElement**_](classtensor__tools_1_1TensorIndexElement.md) _._ |












































## Detailed Description


A namespace to group all the tools that are useful to carry out non-trivial operations on tensors. 


    
## Public Types Documentation




### typedef extract\_sub\_tensor\_element\_t 

_Extract the relevant elements of a_ [_**TensorIndexElement**_](classtensor__tools_1_1TensorIndexElement.md) _to create a sub-TensorIndexElement using a global and a local TypeSeq of VectorIndexIdMaps to identify the relevant elements. For example: for GlobalTensorIndexElement = TensorIndexElement&lt;X,Y&gt; with TypeSeqVectorIndexIdMapGlobal = TypeSeq&lt;_[_**VectorIndexIdMap**_](structtensor__tools_1_1VectorIndexIdMap.md) _&lt;'i', VectorIndexSet&lt;X, Y&gt;&gt;,_[_**VectorIndexIdMap**_](structtensor__tools_1_1VectorIndexIdMap.md) _&lt;'j', VectorIndexSet&lt;X, Y&gt;&gt;&gt; and TypeSeqVectorIndexIdMapLocal = TypeSeq&lt;_[_**VectorIndexIdMap**_](structtensor__tools_1_1VectorIndexIdMap.md) _&lt;'j', VectorIndexSet&lt;X, Y&gt;&gt;&gt; we obtain TensorIndexElement&lt;Y&gt;_
```C++
using tensor_tools::extract_sub_tensor_element_t = typedef typename details::ExtractSubTensorElement< TypeSeqVectorIndexIdMapGlobal, TypeSeqVectorIndexIdMapLocal, GlobalTensorIndexElement>::type;
```





**Template parameters:**


* `TypeSeqVectorIndexIdMapGlobal` The global TypeSeq of VectorIndexIdMaps describing how the [**TensorIndexElement**](classtensor__tools_1_1TensorIndexElement.md) was indexed. 
* `TypeSeqVectorIndexIdMapLocal` The local TypeSeq of VectorIndexIdMaps describing how to identify the relevant indices for the output. 
* `GlobalTensorIndexElement` The starting [**TensorIndexElement**](classtensor__tools_1_1TensorIndexElement.md). 




        

<hr>



### typedef get\_nth\_tensor\_index\_element\_from\_map\_t 

_Get the n-th valid index for a tensor which is accessed according to the pattern described by a TypeSeq of VectorIndexIdMaps. E.g. for a 2D tensor with components A\_{xx}, A\_{xy}, A\_{yx}, A\_{yy}, indexed with._ 
```C++
using tensor_tools::get_nth_tensor_index_element_from_map_t = typedef to_tensor_index_element_t< get_type_seq_vector_index_set_t<TypeSeqVectorIndexIdMap>, typename details::GetNthTensorIndexElementFromMap< TypeSeqVectorIndexIdMap, get_nth_tensor_index_element_t< Elem, get_type_seq_vector_index_set_t<unique_indices_t<TypeSeqVectorIndexIdMap> >>, ddc::type_seq_size_v<TypeSeqVectorIndexIdMap> >::type>;
```




```C++
TypeSeq<VectorIndexIdMap<'i', VectorIndexSet<X, Y>>, VectorIndexIdMap<'i', VectorIndexSet<X, Y>>>
```




* the 1st element is [**TensorIndexElement**](classtensor__tools_1_1TensorIndexElement.md)&lt;TypeSeq&lt;VectorIndexSet&lt;X, Y&gt;, VectorIndexSet&lt;X, Y&gt;&gt;, [**X**](structX.md), [**X**](structX.md)&gt;
* the 2nd element is [**TensorIndexElement**](classtensor__tools_1_1TensorIndexElement.md)&lt;TypeSeq&lt;VectorIndexSet&lt;X, Y&gt;, VectorIndexSet&lt;X, Y&gt;&gt;, [**Y**](structY.md), [**Y**](structY.md)&gt;




A\_{xy} and A\_{yx} are not valid components as they do not respect the index pattern.




**Template parameters:**


* `Elem` The element of interest. 
* `TypeSeqVectorIndexIdMap` A TypeSeq containing a [**VectorIndexIdMap**](structtensor__tools_1_1VectorIndexIdMap.md) for each dimension of the tensor. 




        

<hr>



### typedef get\_nth\_tensor\_index\_element\_t 

_Get the_ [_**TensorIndexElement**_](classtensor__tools_1_1TensorIndexElement.md) _which indexes a_[_**Tensor**_](classTensor.md) _at the n-th position of its internal array. E.g. for a 2x2_[_**Tensor**_](classTensor.md) _, get\_nth\_tensor\_index\_element\_t&lt;1, TypeSeqVectorIndexSet&gt; returns the_[_**TensorIndexElement**_](classtensor__tools_1_1TensorIndexElement.md) _which indexes element 1 of the array, so the element {0,1} of the tensor._
```C++
using tensor_tools::get_nth_tensor_index_element_t = typedef typename details::GetNthTensorIndexElement< IndexPosition, ddc::type_seq_size_v<TypeSeqVectorIndexSet>, TypeSeqVectorIndexSet>::type;
```





**Template parameters:**


* `IndexPosition` The index of the underlying array for which we want to collect the [**TensorIndexElement**](classtensor__tools_1_1TensorIndexElement.md). 
* `TypeSeqVectorIndexSet` A TypeSeq containing the VectorIndexSets describing the valid indices along each dimension of the tensor. 




        

<hr>



### typedef get\_type\_seq\_vector\_index\_set\_t 

_Get a TypeSeq of valid VectorIndexSets from a TypeSeq of VectorIndexIdMaps._ 
```C++
using tensor_tools::get_type_seq_vector_index_set_t = typedef typename details::ExtractTypeSeqIndexSet<TypeSeqVectorIndexIdMap>::type;
```





**Template parameters:**


* `TypeSeqVectorIndexIdMap` A TypeSeq containing a [**VectorIndexIdMap**](structtensor__tools_1_1VectorIndexIdMap.md) for each dimension of the tensor. 




        

<hr>



### typedef index\_identifiers\_t 

_Create a TypeSeq of integral\_constants containing all the character ids found in the TypeSeqVectorIndexIdMap._ 
```C++
using tensor_tools::index_identifiers_t = typedef typename details::GetIndexIds<TypeSeqVectorIndexIdMap>::type;
```





**Template parameters:**


* `TypeSeqVectorIndexIdMap` A TypeSeq containing a [**VectorIndexIdMap**](structtensor__tools_1_1VectorIndexIdMap.md) for each dimension of the tensor. 




        

<hr>



### typedef non\_repeated\_indices\_t 

_Extract the VectorIndexIdMaps whose character id only appears once in a TypeSeq of VectorIndexIdMaps._ 
```C++
using tensor_tools::non_repeated_indices_t = typedef typename details::GetNonRepeatedIndices< TypeSeqVectorIndexIdMap, ddc::type_seq_size_v<TypeSeqVectorIndexIdMap> - 1, index_identifiers_t<TypeSeqVectorIndexIdMap> >::type;
```





**Template parameters:**


* `TypeSeqVectorIndexIdMap` A TypeSeq containing a [**VectorIndexIdMap**](structtensor__tools_1_1VectorIndexIdMap.md) for each dimension of the tensor. 




        

<hr>



### typedef relevant\_vector\_index\_sets\_t 

_Extract all VectorIndexSets which are described by the same character identifier._ 
```C++
using tensor_tools::relevant_vector_index_sets_t = typedef typename details::GetRelevantVectorIndexSets<ID, TypeSeqVectorIndexIdMap>::type;
```





**Template parameters:**


* `ID` The character identifying the VectorIndexSets we are searching for. 
* `TypeSeqVectorIndexIdMap` A TypeSeq containing a [**VectorIndexIdMap**](structtensor__tools_1_1VectorIndexIdMap.md) for each dimension of the tensor. 




        

<hr>



### typedef to\_tensor\_index\_element\_t 

_Get a_ [_**TensorIndexElement**_](classtensor__tools_1_1TensorIndexElement.md) _from a TypeSeq of valid VectorIndexSets and a TypeSeq of indices._
```C++
using tensor_tools::to_tensor_index_element_t = typedef typename details:: ToTensorIndexElement<ValidatingTensorIndexSet, TypeSeqTensorIndexTag>::type;
```





**Template parameters:**


* `ValidatingTensorIndexSet` A TypeSeq containing the VectorIndexSets describing the tags that can be used as indices in each dimension. 
* `TypeSeqTensorIndexTag` A TypeSeq containing the tags used to index the tensor. 




        

<hr>



### typedef unique\_indices\_t 

_Create a TypeSeq of_ [_**VectorIndexIdMap**_](structtensor__tools_1_1VectorIndexIdMap.md) _in which each character id only appears once from a TypeSeq of VectorIndexIdMaps with repeat character ids._
```C++
using tensor_tools::unique_indices_t = typedef typename details::GetUniqueIndices<TypeSeqVectorIndexIdMap>::type;
```





**Template parameters:**


* `TypeSeqVectorIndexIdMap` A TypeSeq containing a [**VectorIndexIdMap**](structtensor__tools_1_1VectorIndexIdMap.md) for each dimension of the tensor. 




        

<hr>
## Public Attributes Documentation




### variable char\_occurrences\_v 

_Count the number of instances of a character in a TypeSeq of literals._ 
```C++
constexpr std::size_t tensor_tools::char_occurrences_v;
```





**Template parameters:**


* `search_char` The character that you are searching for. 
* `TupleType` The type of the TypeSeq of integral\_constants of characters. 




        

<hr>



### variable enable\_indexed\_tensor 

_A boolean, true if the type is an_ [_**IndexedTensor**_](classtensor__tools_1_1IndexedTensor.md) _, false otherwise._
```C++
constexpr bool tensor_tools::enable_indexed_tensor;
```




<hr>



### variable enable\_indexed\_tensor&lt; IndexedTensor&lt; TensorType, TypeSeqVectorIndexIdMap &gt; &gt; 

_A boolean, true if the type is an_ [_**IndexedTensor**_](classtensor__tools_1_1IndexedTensor.md) _._
```C++
constexpr bool tensor_tools::enable_indexed_tensor< IndexedTensor< TensorType, TypeSeqVectorIndexIdMap > >;
```




<hr>



### variable is\_indexed\_tensor\_v 

_A tool to check if a type is an_ [_**IndexedTensor**_](classtensor__tools_1_1IndexedTensor.md) _._
```C++
constexpr bool tensor_tools::is_indexed_tensor_v;
```




<hr>



### variable is\_tensor\_index\_element\_v 

_A helper to check if a type is a_ [_**TensorIndexElement**_](classtensor__tools_1_1TensorIndexElement.md) _._
```C++
constexpr bool tensor_tools::is_tensor_index_element_v;
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/data_types/indexed_tensor.hpp`

