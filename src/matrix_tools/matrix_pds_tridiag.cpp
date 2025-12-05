#include <cassert>
#include <cmath>

#include <string.h>

#if __has_include(<mkl_lapacke.h>)
#include <mkl_lapacke.h>
#else
#include <lapacke.h>
#endif

#include "matrix_pds_tridiag.hpp"

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
    LAPACKE_dpttrf(&n, d.get(), l.get(), &info);
    return info;
}

int Matrix_PDS_Tridiag::solve_inplace_method(double* b, char, int const n_equations) const
{
    int info;
    LAPACKE_dpttrs(&n, &n_equations, d.get(), l.get(), b, &n, &info);
    return info;
}
