

# Struct MPIDim

**template &lt;class DistributedDim&gt;**



[**ClassList**](annotated.md) **>** [**MPIDim**](structMPIDim.md)



_An internal tag used to dsecribe an artificial dimension describing the MPI rank where the scattered information will be sent to or where the gathered information will be collected from._ 

* `#include <mpitools.hpp>`

















## Public Types

| Type | Name |
| ---: | :--- |
| typedef DistributedDim | [**distributed\_dim**](#typedef-distributed_dim)  <br>_The dimension which is distributed over the MPI ranks._  |
















































## Public Types Documentation




### typedef distributed\_dim 

_The dimension which is distributed over the MPI ranks._ 
```C++
using MPIDim< DistributedDim >::distributed_dim =  DistributedDim;
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/mpi_parallelisation/mpitools.hpp`

