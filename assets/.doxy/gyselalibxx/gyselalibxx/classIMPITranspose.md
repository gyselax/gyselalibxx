

# Class IMPITranspose

**template &lt;class Layout1, class Layout2&gt;**



[**ClassList**](annotated.md) **>** [**IMPITranspose**](classIMPITranspose.md)



_A superclass describing an operator for converting from/to different MPI layouts._ [More...](#detailed-description)

* `#include <impitranspose.hpp>`





Inherited by the following classes: [MPITransposeAllToAll](classMPITransposeAllToAll.md),  [MPITransposeAllToAll](classMPITransposeAllToAll.md)












## Public Types

| Type | Name |
| ---: | :--- |
| typedef typename Layout1::discrete\_domain\_type | [**idx\_range\_type1**](#typedef-idx_range_type1)  <br>_The type of the index range of the first MPI layout._  |
| typedef typename Layout2::discrete\_domain\_type | [**idx\_range\_type2**](#typedef-idx_range_type2)  <br>_The type of the index range of the second MPI layout._  |




















## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**IMPITranspose**](#function-impitranspose) (MPI\_Comm comm) <br>_A constructor for the class._  |








## Protected Attributes

| Type | Name |
| ---: | :--- |
|  MPI\_Comm | [**m\_comm**](#variable-m_comm)  <br>_The MPI communicator._  |




















## Detailed Description




**Template parameters:**


* `Layout1` One of the MPI layouts. 
* `Layout2` The other MPI layouts. 




    
## Public Types Documentation




### typedef idx\_range\_type1 

_The type of the index range of the first MPI layout._ 
```C++
using IMPITranspose< Layout1, Layout2 >::idx_range_type1 =  typename Layout1::discrete_domain_type;
```




<hr>



### typedef idx\_range\_type2 

_The type of the index range of the second MPI layout._ 
```C++
using IMPITranspose< Layout1, Layout2 >::idx_range_type2 =  typename Layout2::discrete_domain_type;
```




<hr>
## Public Functions Documentation




### function IMPITranspose 

_A constructor for the class._ 
```C++
inline explicit IMPITranspose::IMPITranspose (
    MPI_Comm comm
) 
```





**Parameters:**


* `comm` The MPI communicator 




        

<hr>
## Protected Attributes Documentation




### variable m\_comm 

_The MPI communicator._ 
```C++
MPI_Comm IMPITranspose< Layout1, Layout2 >::m_comm;
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/mpi_parallelisation/impitranspose.hpp`

