

# File type\_seq\_tools.hpp



[**FileList**](files.md) **>** [**src**](dir_68267d1309a1af8e8297ef4c3efbcdba.md) **>** [**utils**](dir_313caf1132e152dd9b58bea13a4052ca.md) **>** [**type\_seq\_tools.hpp**](type__seq__tools_8hpp.md)

[Go to the source code of this file](type__seq__tools_8hpp_source.md)




















## Public Types

| Type | Name |
| ---: | :--- |
| typedef typename detail::FindAllGrids&lt; TypeSeqDim, TypeSeqGrid &gt;::type | [**find\_all\_grids\_t**](#typedef-find_all_grids_t)  <br>_A tool to find the a sequence of grids that are defined along the specified dimensions (e.g. get_ [_**GridX**_](structGridX.md) _from_[_**X**_](structX.md) _)_ |
| typedef typename detail::FindGrid&lt; Dim, TypeSeqGrid &gt;::type | [**find\_grid\_t**](#typedef-find_grid_t)  <br>_A tool to find the grid that is defined along the specified dimension (e.g. get_ [_**GridX**_](structGridX.md) _from_[_**X**_](structX.md) _)_ |
| typedef typename detail::FindIdxType&lt; CoordType, IdxRangeType &gt;::type | [**find\_idx\_t**](#typedef-find_idx_t)  <br>_Find the type of an index which allows access to a Coordinate of the specified type._  |
| typedef typename detail::TypeSeqCat&lt; TypeSeqs... &gt;::type | [**type\_seq\_cat\_t**](#typedef-type_seq_cat_t)  <br>_Concatenate type sequences into a new type sequence. This is similar to type\_seq\_merge\_t but it does not remove duplicate elements._  |
| typedef typename detail::TypeSeqDuplicate&lt; Element, std::make\_index\_sequence&lt; n\_elements &gt; &gt;::type | [**type\_seq\_duplicate\_t**](#typedef-type_seq_duplicate_t)  <br>_Create a type sequence containing the element Element, repeated n times._  |
| typedef typename detail::TypeSeqRange&lt; TypeSeqIn, Start, End, Start &gt;::type | [**type\_seq\_range\_t**](#typedef-type_seq_range_t)  <br>_A tool to get a subset of a TypeSeq by slicing [Start:End]._  |
| typedef typename detail::GetUnique&lt; StartTypeSeq &gt;::type | [**type\_seq\_unique\_t**](#typedef-type_seq_unique_t)  <br>_Get a TypeSeq containing all unique types from the original TypeSeq._  |




## Public Attributes

| Type | Name |
| ---: | :--- |
|  constexpr bool | [**type\_seq\_has\_unique\_elements\_v**](#variable-type_seq_has_unique_elements_v)   = `std::is\_same\_v&lt;TypeSeqType, type\_seq\_unique\_t&lt;TypeSeqType&gt;&gt;`<br>_Determine if a type sequence only contains unique elements._  |
|  constexpr int | [**type\_seq\_permutation\_parity\_v**](#variable-type_seq_permutation_parity_v)   = `detail::GetPermutationParity&lt;TypeSeqType, OrderedTypeSeq&gt;::value`<br>_Determine if the permutation parity of a type sequence._  |












































## Public Types Documentation




### typedef find\_all\_grids\_t 

_A tool to find the a sequence of grids that are defined along the specified dimensions (e.g. get_ [_**GridX**_](structGridX.md) _from_[_**X**_](structX.md) _)_
```C++
using find_all_grids_t =  typename detail::FindAllGrids<TypeSeqDim, TypeSeqGrid>::type;
```




<hr>



### typedef find\_grid\_t 

_A tool to find the grid that is defined along the specified dimension (e.g. get_ [_**GridX**_](structGridX.md) _from_[_**X**_](structX.md) _)_
```C++
using find_grid_t =  typename detail::FindGrid<Dim, TypeSeqGrid>::type;
```




<hr>



### typedef find\_idx\_t 

_Find the type of an index which allows access to a Coordinate of the specified type._ 
```C++
using find_idx_t =  typename detail::FindIdxType<CoordType, IdxRangeType>::type;
```





**Template parameters:**


* `CoordType` The type of the coordinate 
* `IdxRangeType` The type of the index range that the index will come from. 




        

<hr>



### typedef type\_seq\_cat\_t 

_Concatenate type sequences into a new type sequence. This is similar to type\_seq\_merge\_t but it does not remove duplicate elements._ 
```C++
using type_seq_cat_t =  typename detail::TypeSeqCat<TypeSeqs...>::type;
```





**Template parameters:**


* `TypeSeqs` The type sequences to be concatenated. 




        

<hr>



### typedef type\_seq\_duplicate\_t 

_Create a type sequence containing the element Element, repeated n times._ 
```C++
using type_seq_duplicate_t =  typename detail::TypeSeqDuplicate<Element, std::make_index_sequence<n_elements> >::type;
```





**Template parameters:**


* `Element` The element to be placed in the type sequence. 
* `n_element` The number of times the element should appear in the type sequence. 




        

<hr>



### typedef type\_seq\_range\_t 

_A tool to get a subset of a TypeSeq by slicing [Start:End]._ 
```C++
using type_seq_range_t =  typename detail::TypeSeqRange<TypeSeqIn, Start, End, Start>::type;
```




<hr>



### typedef type\_seq\_unique\_t 

_Get a TypeSeq containing all unique types from the original TypeSeq._ 
```C++
using type_seq_unique_t =  typename detail::GetUnique<StartTypeSeq>::type;
```





**Template parameters:**


* `StartTypeSeq` The original TypeSeq which may contain duplicate types. 




        

<hr>
## Public Attributes Documentation




### variable type\_seq\_has\_unique\_elements\_v 

_Determine if a type sequence only contains unique elements._ 
```C++
constexpr bool type_seq_has_unique_elements_v;
```





**Template parameters:**


* `TypeSeqType` The type sequence being examined. 




        

<hr>



### variable type\_seq\_permutation\_parity\_v 

_Determine if the permutation parity of a type sequence._ 
```C++
constexpr int type_seq_permutation_parity_v;
```





**Template parameters:**


* `TypeSeqType` The type sequence whose permutation parity is calculated. 
* `OrderedTypeSeq` The final order of the indices.whose permutation parity is calculated. 
* `OrderedTypeSeq` The final order of the indices. 




        

<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/utils/type_seq_tools.hpp`

