

# File matrix\_batch\_ell.hpp



[**FileList**](files.md) **>** [**matrix\_tools**](dir_8cedd1260cc2f2819c8df2fc66ad98b5.md) **>** [**matrix\_batch\_ell.hpp**](matrix__batch__ell_8hpp.md)

[Go to the source code of this file](matrix__batch__ell_8hpp_source.md)



* `#include <ginkgo/extensions/kokkos.hpp>`
* `#include <ginkgo/ginkgo.hpp>`
* `#include <Kokkos_Core.hpp>`
* `#include "matrix_batch.hpp"`
* `#include "matrix_utils.hpp"`















## Classes

| Type | Name |
| ---: | :--- |
| class | [**MatrixBatchEll**](classMatrixBatchEll.md) &lt;class ExecSpace&gt;<br>[_**Matrix**_](classMatrix.md) _class which is able to manage and solve a batch of sparse linear systems. Executes on either CPU or GPU. It takes advantage of the sparse structure, and the only batched solver available in Ginkgo : Stabilised Bicg. The sparsity pattern is assumed to be the same for all matrices. ie the non-zero components are located at the same places for all matrices. This class uses the ELL storage format which needs two 1D arrays, one stores values the other column indices. The class returns these arrays (as Kokkos views) with the get\_batch\_idx\_and\_vals function, it is then possible to fill them outside the class. Tolerance and maximal number of iterations, which are parameters for the iterative solver, are set in the constructor. It is possible to get convergence information by activating the logger at constructor call._ |



















































------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/matrix_tools/matrix_batch_ell.hpp`

