

# File mpitools.hpp



[**FileList**](files.md) **>** [**mpi\_parallelisation**](dir_a35b8fd75f8fad0c2619b083ab571e51.md) **>** [**mpitools.hpp**](mpitools_8hpp.md)

[Go to the source code of this file](mpitools_8hpp_source.md)



* `#include <mpi.h>`
* `#include "ddc_aliases.hpp"`















## Classes

| Type | Name |
| ---: | :--- |
| struct | [**MPIDim**](structMPIDim.md) &lt;class DistributedDim&gt;<br>_An internal tag used to dsecribe an artificial dimension describing the MPI rank where the scattered information will be sent to or where the gathered information will be collected from._  |


## Public Types

| Type | Name |
| ---: | :--- |
| typedef typename detail::InsertIntoTypeSeq&lt; TagToInsert, ddc::type\_seq\_rank\_v&lt; ddc::type\_seq\_element\_t&lt; 0, SubSeq &gt;, TypeSeq &gt;, TypeSeq &gt;::type | [**insert\_into\_seq\_before\_t**](#typedef-insert_into_seq_before_t)  <br>_A tool to insert a tag into an existing TypeSeq immediately preceding an existing subset of tags._  |
| typedef typename detail::InsertIntoTypeSeq&lt; TagToInsert, POS, TypeSeq &gt;::type | [**insert\_into\_type\_seq\_t**](#typedef-insert_into_type_seq_t)  <br>_A tool to insert a tag into an existing TypeSeq at a specified position._  |
| typedef typename detail::InsertMPITags&lt; MPISeq, TypeSeq &gt;::type | [**insert\_mpi\_tags\_into\_seq\_t**](#typedef-insert_mpi_tags_into_seq_t)  <br>_Insert MPI distribution tags into an existing TypeSeq._  |




## Public Attributes

| Type | Name |
| ---: | :--- |
|  constexpr bool | [**dimensions\_are\_adjacent\_v**](#variable-dimensions_are_adjacent_v)   = `detail::DimensionsAreAdjacent&lt;Query, ContainerTypeSeq&gt;::value`<br>_A helper constant to determine if a set of types found within a TypeSeq are adjacent to one another._  |


## Public Static Attributes

| Type | Name |
| ---: | :--- |
|  const MPI\_Datatype | [**MPI\_type\_descriptor\_t**](#variable-mpi_type_descriptor_t)   = `detail::MPITypeDescriptor&lt;ElementType&gt;::get\_type()`<br>_A helper to get the MPI type descriptor from an element type._  |










































## Public Types Documentation




### typedef insert\_into\_seq\_before\_t 

_A tool to insert a tag into an existing TypeSeq immediately preceding an existing subset of tags._ 
```C++
using insert_into_seq_before_t =  typename detail::InsertIntoTypeSeq< TagToInsert, ddc::type_seq_rank_v<ddc::type_seq_element_t<0, SubSeq>, TypeSeq>, TypeSeq>::type;
```





**Template parameters:**


* `TagToInsert` The tag to be inserted. 
* `SubSeq` The subset of tags found in the TypeSeq before which the tag should be inserted. 
* `TypeSeq` The TypeSeq into which the tag should be inserted. 




        

<hr>



### typedef insert\_into\_type\_seq\_t 

_A tool to insert a tag into an existing TypeSeq at a specified position._ 
```C++
using insert_into_type_seq_t =  typename detail::InsertIntoTypeSeq<TagToInsert, POS, TypeSeq>::type;
```





**Template parameters:**


* `TagToInsert` The tag to be inserted. 
* `POS` The position at which the tag should be inserted. 
* `TypeSeq` The TypeSeq into which the tag should be inserted. 




        

<hr>



### typedef insert\_mpi\_tags\_into\_seq\_t 

_Insert MPI distribution tags into an existing TypeSeq._ 
```C++
using insert_mpi_tags_into_seq_t =  typename detail::InsertMPITags<MPISeq, TypeSeq>::type;
```



The MPI tags are each associated with an index range. The MPI tags are inserted into the TypeSeq immediately preceding the tag with which they are associated. This allows an index range to be split along the axes on which it will be scattered.


E.g. if MPI&lt;Phi&gt; is inserted into &lt;[**R**](structR.md), [**Theta**](structTheta.md), Phi, VPar&gt; we would obtain: &lt;[**R**](structR.md), [**Theta**](structTheta.md), MPI&lt;Phi&gt;, Phi, VPar&gt;




**Template parameters:**


* `MPISeq` A type sequence containing the MPI tags which should be inserted 
* `TypeSeq` The type sequence which should be inserted. 




        

<hr>
## Public Attributes Documentation




### variable dimensions\_are\_adjacent\_v 

_A helper constant to determine if a set of types found within a TypeSeq are adjacent to one another._ 
```C++
constexpr bool dimensions_are_adjacent_v;
```





**Template parameters:**


* `Query` The TypeSeq describing the possibly adjacent subset. 
* `ContainerTypeSeq` The TypeSeq which contains the subset. 




        

<hr>
## Public Static Attributes Documentation




### variable MPI\_type\_descriptor\_t 

_A helper to get the MPI type descriptor from an element type._ 
```C++
const MPI_Datatype MPI_type_descriptor_t;
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/mpi_parallelisation/mpitools.hpp`

