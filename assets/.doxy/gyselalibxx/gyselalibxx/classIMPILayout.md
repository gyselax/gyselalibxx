

# Class IMPILayout

**template &lt;class IdxRangeData, class... DistributedDim&gt;**



[**ClassList**](annotated.md) **>** [**IMPILayout**](classIMPILayout.md)



_A super class describing a way in which data may be laid out across MPI processes._ [More...](#detailed-description)

* `#include <impilayout.hpp>`

















## Public Types

| Type | Name |
| ---: | :--- |
| typedef IdxRangeData | [**discrete\_domain\_type**](#typedef-discrete_domain_type)  <br>_The index range of the data._  |
| typedef IdxRange&lt; DistributedDim... &gt; | [**distributed\_sub\_idx\_range**](#typedef-distributed_sub_idx_range)  <br>_The index range of the distributed section of the data._  |
| typedef ddc::detail::TypeSeq&lt; DistributedDim... &gt; | [**distributed\_type\_seq**](#typedef-distributed_type_seq)  <br>_A type sequence describing the dimensions which are distributed across MPI processes._  |






## Public Static Attributes

| Type | Name |
| ---: | :--- |
|  constexpr bool | [**distributed\_idx\_ranges\_are\_first**](#variable-distributed_idx_ranges_are_first)   = `/* multi line expression */`<br>_A flag to indicate whether the distributed dimensions are the dimensions which are the furthest from being contiguous in memory._  |
|  constexpr int | [**n\_distributed\_dimensions**](#variable-n_distributed_dimensions)   = `sizeof...(DistributedDim)`<br>_The number of dimensions that are distributed across MPI processes._  |










































## Detailed Description




**Template parameters:**


* `IdxRangeData` The index range on which the data is defined. 
* `DistributedDim` The tags of the discrete dimensions which are distributed across MPI processes. 




    
## Public Types Documentation




### typedef discrete\_domain\_type 

_The index range of the data._ 
```C++
using IMPILayout< IdxRangeData, DistributedDim >::discrete_domain_type =  IdxRangeData;
```




<hr>



### typedef distributed\_sub\_idx\_range 

_The index range of the distributed section of the data._ 
```C++
using IMPILayout< IdxRangeData, DistributedDim >::distributed_sub_idx_range =  IdxRange<DistributedDim...>;
```




<hr>



### typedef distributed\_type\_seq 

_A type sequence describing the dimensions which are distributed across MPI processes._ 
```C++
using IMPILayout< IdxRangeData, DistributedDim >::distributed_type_seq =  ddc::detail::TypeSeq<DistributedDim...>;
```




<hr>
## Public Static Attributes Documentation




### variable distributed\_idx\_ranges\_are\_first 

_A flag to indicate whether the distributed dimensions are the dimensions which are the furthest from being contiguous in memory._ 
```C++
constexpr bool IMPILayout< IdxRangeData, DistributedDim >::distributed_idx_ranges_are_first;
```




<hr>



### variable n\_distributed\_dimensions 

_The number of dimensions that are distributed across MPI processes._ 
```C++
constexpr int IMPILayout< IdxRangeData, DistributedDim >::n_distributed_dimensions;
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/mpi_parallelisation/impilayout.hpp`

