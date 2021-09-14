#include <cassert>
#include <cmath>

#include <string.h>

#include "sll/matrix_pds_tridiag.h"

extern "C" int dpttrf_(const int* n, double* d, double* e, int* info);
extern "C" int dpttrs_(
        const int* n,
        const int* nrhs,
        double* d,
        double* e,
        double* b,
        const int* ldb,
        int* info);

Matrix_PDS_Tridiag::Matrix_PDS_Tridiag(int n)
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
        int tmp(i);
        i = j;
        j = tmp;
    }
    if (i + 1 == j) {
        return l[i];
    }
    return 0.0;
}

void Matrix_PDS_Tridiag::set_element(int i, int j, double a_ij)
{
    if (i == j) {
        d[i] = a_ij;
        return;
    }
    if (i > j) {
        int tmp(i);
        i = j;
        j = tmp;
    }
    if (i + 1 != j) {
        assert(fabs(a_ij) < 1e-20);
    } else {
        l[i] = a_ij;
    }
}

int Matrix_PDS_Tridiag::factorize_method()
{
    int info;
    dpttrf_(&n, d.get(), l.get(), &info);
    return info;
}

int Matrix_PDS_Tridiag::solve_inplace_method(const char transpose, double* b, int nrows, int ncols)
        const
{
    int info;
    dpttrs_(&n, &ncols, d.get(), l.get(), b, &nrows, &info);
    return info;
}
