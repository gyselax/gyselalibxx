

# File matrix\_batch\_tridiag.hpp



[**FileList**](files.md) **>** [**matrix\_tools**](dir_8cedd1260cc2f2819c8df2fc66ad98b5.md) **>** [**matrix\_batch\_tridiag.hpp**](matrix__batch__tridiag_8hpp.md)

[Go to the source code of this file](matrix__batch__tridiag_8hpp_source.md)



* `#include <Kokkos_Core.hpp>`
* `#include "matrix_batch.hpp"`















## Classes

| Type | Name |
| ---: | :--- |
| class | [**MatrixBatchTridiag**](classMatrixBatchTridiag.md) &lt;class ExecSpace&gt;<br>_A structure for solving a set of independent tridiagonal systems using a direct method. The parallelism operates on the whole collection by dispatching to threads. Each problem is treated sequentially, by the tridiagonal matrix algorithm (TDMA). This solver is stable for tridiagonal matrices which satisfy one of the following conditions:_  |



















































------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/matrix_tools/matrix_batch_tridiag.hpp`

