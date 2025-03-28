

# File matrix\_pds\_tridiag.cpp

[**File List**](files.md) **>** [**matrix\_tools**](dir_8cedd1260cc2f2819c8df2fc66ad98b5.md) **>** [**matrix\_pds\_tridiag.cpp**](matrix__pds__tridiag_8cpp.md)

[Go to the documentation of this file](matrix__pds__tridiag_8cpp.md)


```C++
#include <cassert>
#include <cmath>

#include <string.h>

#include "matrix_pds_tridiag.hpp"

extern "C" int dpttrf_(int const* n, double* d, double* e, int* info);
extern "C" int dpttrs_(
        int const* n,
        int const* nrhs,
        double* d,
        double* e,
        double* b,
        int const* ldb,
        int* info);

Matrix_PDS_Tridiag::Matrix_PDS_Tridiag(int const n)
    : Matrix(n)
    , d(std::make_unique<double[]>(n))
    , l(std::make_unique<double[]>(n - 1))
{
    memset(d.get(), 0, n * sizeof(double));
    memset(l.get(), 0, (n - 1) * sizeof(double));
}

double Matrix_PDS_Tridiag::get_element(int i, int j) const
{
    if (i == j) {
        return d[i];
    }
    if (i > j) {
        std::swap(i, j);
    }
    if (i + 1 == j) {
        return l[i];
    }
    return 0.0;
}

void Matrix_PDS_Tridiag::set_element(int i, int j, double const a_ij)
{
    if (i == j) {
        d[i] = a_ij;
        return;
    }
    if (i > j) {
        std::swap(i, j);
    }
    if (i + 1 != j) {
        assert(std::fabs(a_ij) < 1e-20);
    } else {
        l[i] = a_ij;
    }
}

int Matrix_PDS_Tridiag::factorise_method()
{
    int info;
    dpttrf_(&n, d.get(), l.get(), &info);
    return info;
}

int Matrix_PDS_Tridiag::solve_inplace_method(double* b, char, int const n_equations) const
{
    int info;
    dpttrs_(&n, &n_equations, d.get(), l.get(), b, &n, &info);
    return info;
}
```


