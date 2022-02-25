#include <cassert>

#include "sll/matrix_dense.hpp"

extern "C" int dgetrf_(int const* m, int const* n, double* a, int const* lda, int* ipiv, int* info);
extern "C" int dgetrs_(
        char const* trans,
        int const* n,
        int const* nrhs,
        double* a,
        int const* lda,
        int* ipiv,
        double* b,
        int const* ldb,
        int* info);

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

int Matrix_Dense::factorize_method()
{
    int info;
    dgetrf_(&n, &n, a.get(), &n, ipiv.get(), &info);
    return info;
}

int Matrix_Dense::solve_inplace_method(double* b, char const transpose, int const n_equations) const
{
    int info;
    dgetrs_(&transpose, &n, &n_equations, a.get(), &n, ipiv.get(), b, &n, &info);
    return info;
}
