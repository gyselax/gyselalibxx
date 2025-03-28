

# Class MPITransposeAllToAll

**template &lt;class Layout1, class Layout2&gt;**



[**ClassList**](annotated.md) **>** [**MPITransposeAllToAll**](classMPITransposeAllToAll.md)



_A class describing an operator for converting from/to different MPI layouts using AlltoAll._ [More...](#detailed-description)

* `#include <mpitransposealltoall.hpp>`



Inherits the following classes: [IMPITranspose](classIMPITranspose.md)














## Public Types

| Type | Name |
| ---: | :--- |
| typedef typename Layout1::distributed\_sub\_idx\_range | [**distributed\_idx\_range\_type1**](#typedef-distributed_idx_range_type1)  <br>_The type of the index range of the first MPI layout._  |
| typedef typename Layout2::distributed\_sub\_idx\_range | [**distributed\_idx\_range\_type2**](#typedef-distributed_idx_range_type2)  <br>_The type of the index range of the second MPI layout._  |
| typedef typename Layout1::discrete\_domain\_type | [**idx\_range\_type1**](#typedef-idx_range_type1)  <br>_The type of the index range of the first MPI layout._  |
| typedef typename Layout2::discrete\_domain\_type | [**idx\_range\_type2**](#typedef-idx_range_type2)  <br>_The type of the index range of the second MPI layout._  |


## Public Types inherited from IMPITranspose

See [IMPITranspose](classIMPITranspose.md)

| Type | Name |
| ---: | :--- |
| typedef typename Layout1::discrete\_domain\_type | [**idx\_range\_type1**](classIMPITranspose.md#typedef-idx_range_type1)  <br>_The type of the index range of the first MPI layout._  |
| typedef typename Layout2::discrete\_domain\_type | [**idx\_range\_type2**](classIMPITranspose.md#typedef-idx_range_type2)  <br>_The type of the index range of the second MPI layout._  |






































## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**MPITransposeAllToAll**](#function-mpitransposealltoall) (IdxRange global\_idx\_range, MPI\_Comm comm) <br>_A constructor for the transpose operator._  |
|  auto | [**get\_local\_idx\_range**](#function-get_local_idx_range) () const<br>_Getter for the local index range._  |
|  void | [**operator()**](#function-operator) (ExecSpace const & execution\_space, Field&lt; ElementType, IdxRangeOut, MemSpace &gt; recv\_field, ConstField&lt; ElementType, InIdxRange, MemSpace &gt; send\_field) const<br>_An operator which transposes from one layout to another._  |
|  void | [**transpose\_to**](#function-transpose_to) (ExecSpace const & execution\_space, Field&lt; ElementType, typename OutLayout::discrete\_domain\_type, MemSpace &gt; recv\_field, ConstField&lt; ElementType, InIdxRange, MemSpace &gt; send\_field) const<br>_An operator which transposes from one layout to another._  |


## Public Functions inherited from IMPITranspose

See [IMPITranspose](classIMPITranspose.md)

| Type | Name |
| ---: | :--- |
|   | [**IMPITranspose**](classIMPITranspose.md#function-impitranspose) (MPI\_Comm comm) <br>_A constructor for the class._  |
















## Protected Attributes inherited from IMPITranspose

See [IMPITranspose](classIMPITranspose.md)

| Type | Name |
| ---: | :--- |
|  MPI\_Comm | [**m\_comm**](classIMPITranspose.md#variable-m_comm)  <br>_The MPI communicator._  |






































## Detailed Description


This class implements a basic AlltoAll operator and currently only works with a basic MPIBlockLayout.




**Template parameters:**


* `Layout1` One of the MPI layouts. 
* `Layout2` The other MPI layouts. 




    
## Public Types Documentation




### typedef distributed\_idx\_range\_type1 

_The type of the index range of the first MPI layout._ 
```C++
using MPITransposeAllToAll< Layout1, Layout2 >::distributed_idx_range_type1 =  typename Layout1::distributed_sub_idx_range;
```




<hr>



### typedef distributed\_idx\_range\_type2 

_The type of the index range of the second MPI layout._ 
```C++
using MPITransposeAllToAll< Layout1, Layout2 >::distributed_idx_range_type2 =  typename Layout2::distributed_sub_idx_range;
```




<hr>



### typedef idx\_range\_type1 

_The type of the index range of the first MPI layout._ 
```C++
using MPITransposeAllToAll< Layout1, Layout2 >::idx_range_type1 =  typename Layout1::discrete_domain_type;
```




<hr>



### typedef idx\_range\_type2 

_The type of the index range of the second MPI layout._ 
```C++
using MPITransposeAllToAll< Layout1, Layout2 >::idx_range_type2 =  typename Layout2::discrete_domain_type;
```




<hr>
## Public Functions Documentation




### function MPITransposeAllToAll 

_A constructor for the transpose operator._ 
```C++
template<class IdxRange>
inline MPITransposeAllToAll::MPITransposeAllToAll (
    IdxRange global_idx_range,
    MPI_Comm comm
) 
```





**Parameters:**


* `global_idx_range` The global index range of the data across processes provided in either layout. 
* `comm` The MPI communicator. 




        

<hr>



### function get\_local\_idx\_range 

_Getter for the local index range._ 
```C++
template<class Layout>
inline auto MPITransposeAllToAll::get_local_idx_range () const
```





**Template parameters:**


* `Layout` The layout whose index range should be retrieved.



**Returns:**

The local index range in the specified MPI layout. 





        

<hr>



### function operator() 

_An operator which transposes from one layout to another._ 
```C++
template<class ElementType, class InIdxRange, class IdxRangeOut, class MemSpace, class ExecSpace>
inline void MPITransposeAllToAll::operator() (
    ExecSpace const & execution_space,
    Field< ElementType, IdxRangeOut, MemSpace > recv_field,
    ConstField< ElementType, InIdxRange, MemSpace > send_field
) const
```



If the dimensions are ordered differently in the index ranges of the two layouts then this function can be used. If the index ranges have the same type then the transpose\_to function must be used to disambiguate.




**Parameters:**


* `execution_space` The execution space (Host/Device) where the code will run. 
* `recv_field` The chunk which will describe the data in the new layout. This data is gathered from the MPI processes. 
* `send_field` The chunk describing the data in the current layout. This data will be scattered to other MPI processes. 




        

<hr>



### function transpose\_to 

_An operator which transposes from one layout to another._ 
```C++
template<class OutLayout, class ElementType, class MemSpace, class ExecSpace, class InIdxRange>
inline void MPITransposeAllToAll::transpose_to (
    ExecSpace const & execution_space,
    Field< ElementType, typename OutLayout::discrete_domain_type, MemSpace > recv_field,
    ConstField< ElementType, InIdxRange, MemSpace > send_field
) const
```



If the dimensions are ordered differently in the index ranges of the two layouts then this function can be used. If the index ranges have the same type then the transpose\_to function must be used to disambiguate.




**Template parameters:**


* `OutLayout` The layout that the data should be transposed to.



**Parameters:**


* `execution_space` The execution space (Host/Device) where the code will run. 
* `recv_field` The chunk which will describe the data in the new layout. This data is gathered from the MPI processes. 
* `send_field` The chunk describing the data in the current layout. This data will be scattered to other MPI processes. 




        

<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/mpi_parallelisation/mpitransposealltoall.hpp`

