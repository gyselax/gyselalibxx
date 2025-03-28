

# Class MPILayout

**template &lt;class IdxRangeData, class... DistributedDim&gt;**



[**ClassList**](annotated.md) **>** [**MPILayout**](classMPILayout.md)



_A class describing a way in which data may be laid out across MPI processes._ [More...](#detailed-description)

* `#include <mpilayout.hpp>`



Inherits the following classes: [IMPILayout](classIMPILayout.md)














## Public Types

| Type | Name |
| ---: | :--- |
| typedef typename [**base\_type::distributed\_sub\_idx\_range**](classIMPILayout.md#typedef-distributed_sub_idx_range) | [**distributed\_sub\_idx\_range**](#typedef-distributed_sub_idx_range)  <br>_The index range of the distributed section of the data._  |
| typedef typename [**base\_type::distributed\_type\_seq**](classIMPILayout.md#typedef-distributed_type_seq) | [**distributed\_type\_seq**](#typedef-distributed_type_seq)  <br>_A type sequence describing the dimensions which are distributed across MPI processes._  |
| typedef IdxRangeData | [**idx\_range\_type**](#typedef-idx_range_type)  <br>_The index range of the data._  |


## Public Types inherited from IMPILayout

See [IMPILayout](classIMPILayout.md)

| Type | Name |
| ---: | :--- |
| typedef IdxRangeData | [**discrete\_domain\_type**](classIMPILayout.md#typedef-discrete_domain_type)  <br>_The index range of the data._  |
| typedef IdxRange&lt; DistributedDim... &gt; | [**distributed\_sub\_idx\_range**](classIMPILayout.md#typedef-distributed_sub_idx_range)  <br>_The index range of the distributed section of the data._  |
| typedef ddc::detail::TypeSeq&lt; DistributedDim... &gt; | [**distributed\_type\_seq**](classIMPILayout.md#typedef-distributed_type_seq)  <br>_A type sequence describing the dimensions which are distributed across MPI processes._  |












## Public Static Attributes inherited from IMPILayout

See [IMPILayout](classIMPILayout.md)

| Type | Name |
| ---: | :--- |
|  constexpr bool | [**distributed\_idx\_ranges\_are\_first**](classIMPILayout.md#variable-distributed_idx_ranges_are_first)   = `/* multi line expression */`<br>_A flag to indicate whether the distributed dimensions are the dimensions which are the furthest from being contiguous in memory._  |
|  constexpr int | [**n\_distributed\_dimensions**](classIMPILayout.md#variable-n_distributed_dimensions)   = `sizeof...(DistributedDim)`<br>_The number of dimensions that are distributed across MPI processes._  |






























## Public Static Functions

| Type | Name |
| ---: | :--- |
|  [**idx\_range\_type**](classMPILayout.md#typedef-idx_range_type) | [**distribute\_idx\_range**](#function-distribute_idx_range) ([**idx\_range\_type**](classMPILayout.md#typedef-idx_range_type) global\_idx\_range, int comm\_size, int rank) <br>_Get the distributed index range which follows the chosen layout._  |
















































## Protected Static Functions

| Type | Name |
| ---: | :--- |
|  IdxRange&lt; HeadTag &gt; | [**internal\_distribute\_idx\_range**](#function-internal_distribute_idx_range-12) (IdxRange&lt; HeadTag &gt; global\_idx\_range, int comm\_size, int rank) <br>_Distribute a 1D index range over the MPI processes._  |
|  IdxRange&lt; HeadTag, Tags... &gt; | [**internal\_distribute\_idx\_range**](#function-internal_distribute_idx_range-22) (IdxRange&lt; HeadTag, Tags... &gt; idx\_range, int comm\_size, int rank) <br>_Distribute the index range over the MPI processes._  |




## Detailed Description


This class describes the simplest way of laying out data. In this layout the data is distributed along dimensions in order until the data is distributed. It is possible that the data may not be distributed along some of the requested dimensions if such distribution is not necessary. For example if we distribute dimensions [**R**](structR.md) and [**Theta**](structTheta.md) of a ([**R**](structR.md), [**Theta**](structTheta.md), Phi, [**Vpar**](structVpar.md), [**Mu**](structMu.md)) grid of size 256 x 1024 x 256 x 64 x 8 over 64 processes, the theta dimension will not be distributed.


The data is distributed such that it is maximally distributed over each dimension (in order) such that the size of each local index range along that dimension is the same for all processes. For example if we distribute dimensions [**X**](structX.md) and [**Y**](structY.md) of a ([**X**](structX.md), [**Y**](structY.md), Z) grid of size (10, 15, 4) over 6 processes, the [**X**](structX.md) dimension will be distributed over 2 processes and the [**Y**](structY.md) dimension will be distributed over 3 processes.




**Template parameters:**


* `IdxRangeData` The IdxRange on which the data is defined. 
* `DistributedDim` The tags of the discrete dimensions which are distributed across MPI processes. 




    
## Public Types Documentation




### typedef distributed\_sub\_idx\_range 

_The index range of the distributed section of the data._ 
```C++
using MPILayout< IdxRangeData, DistributedDim >::distributed_sub_idx_range =  typename base_type::distributed_sub_idx_range;
```




<hr>



### typedef distributed\_type\_seq 

_A type sequence describing the dimensions which are distributed across MPI processes._ 
```C++
using MPILayout< IdxRangeData, DistributedDim >::distributed_type_seq =  typename base_type::distributed_type_seq;
```




<hr>



### typedef idx\_range\_type 

_The index range of the data._ 
```C++
using MPILayout< IdxRangeData, DistributedDim >::idx_range_type =  IdxRangeData;
```




<hr>
## Public Static Functions Documentation




### function distribute\_idx\_range 

_Get the distributed index range which follows the chosen layout._ 
```C++
static inline idx_range_type MPILayout::distribute_idx_range (
    idx_range_type global_idx_range,
    int comm_size,
    int rank
) 
```





**Parameters:**


* `global_idx_range` The global (non-distributed) index range. 
* `comm_size` The number of MPI processes over which the data is distributed. 
* `rank` The rank of the current MPI process.



**Returns:**

The distributed index range. 





        

<hr>
## Protected Static Functions Documentation




### function internal\_distribute\_idx\_range [1/2]

_Distribute a 1D index range over the MPI processes._ 
```C++
template<class HeadTag>
static inline IdxRange< HeadTag > MPILayout::internal_distribute_idx_range (
    IdxRange< HeadTag > global_idx_range,
    int comm_size,
    int rank
) 
```





**Parameters:**


* `global_idx_range` The index range to be distributed. 
* `comm_size` The number of processes over which the data should be disctributed. 
* `rank` The rank of the process within the current group of processes



**Returns:**

The distributed index range. 





        

<hr>



### function internal\_distribute\_idx\_range [2/2]

_Distribute the index range over the MPI processes._ 
```C++
template<class HeadTag, class... Tags, std::enable_if_t<(sizeof...(Tags) > 0), bool >>
static inline IdxRange< HeadTag, Tags... > MPILayout::internal_distribute_idx_range (
    IdxRange< HeadTag, Tags... > idx_range,
    int comm_size,
    int rank
) 
```



This function is called recursively. At each pass it distributes the index range over the first dimension in the discrete index range. The remaining dimensions and processes are then handled in the recursive call.




**Parameters:**


* `idx_range` The index range to be distributed. 
* `comm_size` The number of processes over which the data should be disctributed. 
* `rank` The rank of the process within the current group of processes



**Returns:**

The distributed index range. 





        

<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/mpi_parallelisation/mpilayout.hpp`

