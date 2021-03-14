#include <cassert>

#include "matrix_banded.h"

extern "C" int dgbtrf_(
        const int* m,
        const int* n,
        const int* kl,
        const int* ku,
        double* a_b,
        const int* lda_b,
        int* ipiv,
        int* info);
extern "C" int dgbtrs_(
        const char* trans,
        const int* n,
        const int* kl,
        const int* ku,
        const int* nrhs,
        double* a_b,
        const int* lda_b,
        int* ipiv,
        double* b,
        const int* ldb,
        int* info);

Matrix_Banded::Matrix_Banded(int n, int kl, int ku) : Matrix(n), kl(kl), ku(ku), c(2 * kl + ku + 1)
{
    assert(n > 0);
    assert(kl >= 0);
    assert(ku >= 0);
    assert(kl <= n);
    assert(ku <= n);
    ipiv = std::make_unique<int[]>(n);
    q = std::make_unique<double[]>(c * n);

    /*
     * Given the linear system A*x=b, we assume that A is a square (n by n)
     * with ku super-diagonals and kl sub-diagonals.
     * All non-zero elements of A are stored in the rectangular matrix q, using
     * the format required by DGBTRF (LAPACK): diagonals of A are rows of q.
     * q has 2*kl rows for the subdiagonals, 1 row for the diagonal, and ku rows
     * for the superdiagonals. (The kl additional rows are needed for pivoting.)
     * The term A(i,j) of the full matrix is stored in q(i-j+2*kl+1,j).
     */
    for (int i(0); i < c * n; ++i) {
        q[i] = 0.0;
    }
}

double Matrix_Banded::get_element(int i, int j) const
{
    if (i >= max(0, j - ku) and i < min(n, j + kl + 1)) {
        return q[i * c + kl + ku + j - i];
    } else {
        return 0.0;
    }
}

void Matrix_Banded::set_element(int i, int j, double a_ij)
{
    if (i >= max(0, j - ku) and i < min(n, j + kl + 1)) {
        q[i * c + kl + ku + j - i] = a_ij;
    } else {
        assert(fabs(a_ij) < 1e-20);
    }
}

int Matrix_Banded::factorize_method()
{
    int info;
    dgbtrf_(&n, &n, &kl, &ku, q.get(), &c, ipiv.get(), &info);
    return info;
}

int Matrix_Banded::solve_inplace_method(const char transpose, double* b, int nrows, int ncols) const
{
    int info;

    dgbtrs_(&transpose, &n, &kl, &ku, &ncols, q.get(), &c, ipiv.get(), b, &nrows, &info);
    return info;
}
