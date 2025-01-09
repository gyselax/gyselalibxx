#include <algorithm>
#include <cassert>
#include <cmath>

#include "matrix_banded.hpp"

extern "C" int dgbtrf_(
        int const* m,
        int const* n,
        int const* kl,
        int const* ku,
        double* a_b,
        int const* lda_b,
        int* ipiv,
        int* info);
extern "C" int dgbtrs_(
        char const* trans,
        int const* n,
        int const* kl,
        int const* ku,
        int const* nrhs,
        double* a_b,
        int const* lda_b,
        int* ipiv,
        double* b,
        int const* ldb,
        int* info);

Matrix_Banded::Matrix_Banded(int const n, int const kl, int const ku)
    : Matrix(n)
    , kl(kl)
    , ku(ku)
    , c(2 * kl + ku + 1)
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
    for (int i = 0; i < c * n; ++i) {
        q[i] = 0.0;
    }
}

double Matrix_Banded::get_element(int const i, int const j) const
{
    if (i >= std::max(0, j - ku) && i < std::min(n, j + kl + 1)) {
        return q[j * c + kl + ku + i - j];
    } else {
        return 0.0;
    }
}

void Matrix_Banded::set_element(int const i, int const j, double const a_ij)
{
    if (i >= std::max(0, j - ku) && i < std::min(n, j + kl + 1)) {
        q[j * c + kl + ku + i - j] = a_ij;
    } else {
        assert(std::fabs(a_ij) < 1e-20);
    }
}

int Matrix_Banded::factorize_method()
{
    int info;
    dgbtrf_(&n, &n, &kl, &ku, q.get(), &c, ipiv.get(), &info);
    return info;
}

int Matrix_Banded::solve_inplace_method(double* b, char const transpose, int const n_equations)
        const
{
    int info;
    dgbtrs_(&transpose, &n, &kl, &ku, &n_equations, q.get(), &c, ipiv.get(), b, &n, &info);
    return info;
}
