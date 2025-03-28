

# File type\_seq\_tools.hpp



[**FileList**](files.md) **>** [**src**](dir_68267d1309a1af8e8297ef4c3efbcdba.md) **>** [**utils**](dir_313caf1132e152dd9b58bea13a4052ca.md) **>** [**type\_seq\_tools.hpp**](type__seq__tools_8hpp.md)

[Go to the source code of this file](type__seq__tools_8hpp_source.md)




















## Public Types

| Type | Name |
| ---: | :--- |
| typedef typename detail::FindGrid&lt; Dim, TypeSeqGrid &gt;::type | [**find\_grid\_t**](#typedef-find_grid_t)  <br>_A tool to find the grid that is defined along the specified dimension (e.g. get_ [_**GridX**_](structGridX.md) _from_[_**X**_](structX.md) _)_ |
| typedef typename detail::TypeSeqCat&lt; TypeSeqs... &gt;::type | [**type\_seq\_cat\_t**](#typedef-type_seq_cat_t)  <br>_Concatenate type sequences into a new type sequence. This is similar to type\_seq\_merge\_t but it does not remove duplicate elements._  |
| typedef typename detail::TypeSeqRange&lt; TypeSeqIn, Start, End, Start &gt;::type | [**type\_seq\_range\_t**](#typedef-type_seq_range_t)  <br>_A tool to get a subset of a TypeSeq by slicing [Start:End]._  |
| typedef typename detail::GetUnique&lt; StartTypeSeq &gt;::type | [**type\_seq\_unique\_t**](#typedef-type_seq_unique_t)  <br>_Get a TypeSeq containing all unique types from the original TypeSeq._  |




## Public Attributes

| Type | Name |
| ---: | :--- |
|  constexpr bool | [**type\_seq\_has\_unique\_elements\_v**](#variable-type_seq_has_unique_elements_v)   = `std::is\_same\_v&lt;TypeSeqType, type\_seq\_unique\_t&lt;TypeSeqType&gt;&gt;`<br>_Determine if a type sequence only contains unique elements._  |












































## Public Types Documentation




### typedef find\_grid\_t 

_A tool to find the grid that is defined along the specified dimension (e.g. get_ [_**GridX**_](structGridX.md) _from_[_**X**_](structX.md) _)_
```C++
using find_grid_t =  typename detail::FindGrid<Dim, TypeSeqGrid>::type;
```




<hr>



### typedef type\_seq\_cat\_t 

_Concatenate type sequences into a new type sequence. This is similar to type\_seq\_merge\_t but it does not remove duplicate elements._ 
```C++
using type_seq_cat_t =  typename detail::TypeSeqCat<TypeSeqs...>::type;
```





**Template parameters:**


* `TypeSeqs` The type sequences to be concatenated. 




        

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

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/utils/type_seq_tools.hpp`

