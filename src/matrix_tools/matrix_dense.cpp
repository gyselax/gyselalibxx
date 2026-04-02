#include <cassert>

#if __has_include(<mkl_lapacke.h>)
#include <mkl_lapacke.h>
#else
#include <lapacke.h>
#endif

#include "matrix_dense.hpp"

Matrix_Dense::Matrix_Dense(int const n) : Matrix(n)
{
    assert(n > 0);
    ipiv = std::make_unique<int[]>(n);
    a = std::make_unique<double[]>(n * n);
    for (int i = 0; i < n * n; ++i) {
        a[i] = 0;
    }
}

void Matrix_Dense::set_element(int const i, int const j, double const aij)
{
    a[j * n + i] = aij;
}

double Matrix_Dense::get_element(int const i, int const j) const
{
    assert(i < n);
    assert(j < n);
    return a[j * n + i];
}

int Matrix_Dense::factorise_method()
{
    return LAPACKE_dgetrf(LAPACK_COL_MAJOR, n, n, a.get(), n, ipiv.get());
}

int Matrix_Dense::solve_inplace_method(double* b, char const transpose, int const n_equations) const
{
    return LAPACKE_dgetrs(
            LAPACK_COL_MAJOR,
            transpose,
            n,
            n_equations,
            a.get(),
            n,
            ipiv.get(),
            b,
            n);
}
