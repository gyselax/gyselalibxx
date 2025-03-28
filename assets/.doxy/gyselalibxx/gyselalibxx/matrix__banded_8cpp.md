

# File matrix\_banded.cpp



[**FileList**](files.md) **>** [**matrix\_tools**](dir_8cedd1260cc2f2819c8df2fc66ad98b5.md) **>** [**matrix\_banded.cpp**](matrix__banded_8cpp.md)

[Go to the source code of this file](matrix__banded_8cpp_source.md)



* `#include <algorithm>`
* `#include <cassert>`
* `#include <cmath>`
* `#include "matrix_banded.hpp"`





































## Public Functions

| Type | Name |
| ---: | :--- |
|  int | [**dgbtrf\_**](#function-dgbtrf_) (int const \* m, int const \* n, int const \* kl, int const \* ku, double \* a\_b, int const \* lda\_b, int \* ipiv, int \* info) <br> |
|  int | [**dgbtrs\_**](#function-dgbtrs_) (char const \* trans, int const \* n, int const \* kl, int const \* ku, int const \* nrhs, double \* a\_b, int const \* lda\_b, int \* ipiv, double \* b, int const \* ldb, int \* info) <br> |




























## Public Functions Documentation




### function dgbtrf\_ 

```C++
int dgbtrf_ (
    int const * m,
    int const * n,
    int const * kl,
    int const * ku,
    double * a_b,
    int const * lda_b,
    int * ipiv,
    int * info
) 
```




<hr>



### function dgbtrs\_ 

```C++
int dgbtrs_ (
    char const * trans,
    int const * n,
    int const * kl,
    int const * ku,
    int const * nrhs,
    double * a_b,
    int const * lda_b,
    int * ipiv,
    double * b,
    int const * ldb,
    int * info
) 
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/matrix_tools/matrix_banded.cpp`

