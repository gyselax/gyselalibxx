#include "matrix_pds_tridiag.h"
#include <string.h>
#include <cassert>

extern "C" int dpttrf_(const int *n, double *d, double* e, int *info);
extern "C" int dpttrs_(const int *n, const int* nrhs,
        double *d, double *e, double *b, const int *ldb, int *info);

Matrix_PDS_Tridiag::Matrix_PDS_Tridiag(int n)
    : Matrix(n), d(new double[n]), l(new double[n-1])
{
    memset(d, 0, n*sizeof(double));
    memset(l, 0, (n-1)*sizeof(double));
}

Matrix_PDS_Tridiag::~Matrix_PDS_Tridiag()
{
    delete[] d;
    delete[] l;
}

double Matrix_PDS_Tridiag::get_element(int i, int j) const
{
    if (i==j)
    {
        return d[i];
    }
    if (i>j)
    {
        int tmp(i);
        i=j;
        j=tmp;
    }
    if (i+1==j)
    {
        return l[i];
    }
    return 0.0;
}

void Matrix_PDS_Tridiag::set_element(int i, int j, double a_ij)
{
    if (i==j)
    {
        d[i] = a_ij;
        return;
    }
    if (i>j)
    {
        int tmp(i);
        i=j;
        j=tmp;
    }
    if (i+1!=j)
    {
        assert(fabs(a_ij) < 1e-20);
    }
    else
    {
        l[i] = a_ij;
    }
}

int Matrix_PDS_Tridiag::factorize_method()
{
    int info;
    dpttrf_( &n, d, l, &info);
    return info;
}

int Matrix_PDS_Tridiag::solve_inplace_method(const char transpose, double* b, int nrows, int ncols) const
{
    int info;
    dpttrs_( &n, &ncols, d, l, b, &nrows, &info);
    return info;
}
