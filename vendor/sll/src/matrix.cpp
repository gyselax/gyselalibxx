#include <algorithm>
#include <cassert>
#include <iomanip>
#include <iostream>
#include <memory>
#include <utility>

#include "sll/matrix.hpp"
#include "sll/matrix_banded.hpp"
#include "sll/matrix_center_block.hpp"
#include "sll/matrix_corner_block.hpp"
#include "sll/matrix_dense.hpp"
#include "sll/matrix_pds_tridiag.hpp"
#include "sll/matrix_periodic_banded.hpp"
#include "sll/view.hpp"

using std::max;
using std::min;

void Matrix::solve_inplace(DSpan1D& b) const
{
    assert(b.extent(0) == n);
    int info = solve_inplace_method('N', b.data(), b.extent(0), 1);

    if (info < 0) {
        std::cerr << -info << "-th argument had an illegal value" << std::endl;
        // TODO: Add LOG_FATAL_ERROR
    }
}

void Matrix::solve_transpose_inplace(DSpan1D& b) const
{
    assert(b.extent(0) == n);
    int info = solve_inplace_method('T', b.data(), b.extent(0), 1);

    if (info < 0) {
        std::cerr << -info << "-th argument had an illegal value" << std::endl;
        // TODO: Add LOG_FATAL_ERROR
    }
}

void Matrix::solve_inplace_matrix(DSpan2D& bx) const
{
    assert(bx.extent(1) == n);
    int info = solve_inplace_method('N', bx.data(), bx.extent(1), bx.extent(0));

    if (info < 0) {
        std::cerr << -info << "-th argument had an illegal value" << std::endl;
        // TODO: Add LOG_FATAL_ERROR
    }
}

void Matrix::factorize()
{
    int info = factorize_method();

    if (info < 0) {
        std::cerr << -info << "-th argument had an illegal value" << std::endl;
        // TODO: Add LOG_FATAL_ERROR
    } else if (info > 0) {
        std::cerr << "U(" << info << "," << info << ") is exactly zero.";
        std::cerr << " The factorization has been completed, but the factor";
        std::cerr << " U is exactly singular, and division by zero will occur "
                     "if "
                     "it is used to";
        std::cerr << " solve a system of equations.";
        // TODO: Add LOG_FATAL_ERROR
    }
}

std::unique_ptr<Matrix> Matrix::make_new_banded(int n, int kl, int ku, bool pds)
{
    if (kl == ku && kl == 1 && pds) {
        return std::make_unique<Matrix_PDS_Tridiag>(n);
    } else if (2 * kl + 1 + ku >= n) {
        return std::make_unique<Matrix_Dense>(n);
    } else {
        return std::make_unique<Matrix_Banded>(n, kl, ku);
    }
}

std::unique_ptr<Matrix> Matrix::make_new_periodic_banded(int n, int kl, int ku, bool pds)
{
    int border_size(max(kl, ku));
    int banded_size(n - border_size);
    std::unique_ptr<Matrix> block_mat;
    if (pds && kl == ku && kl == 1) {
        block_mat = std::make_unique<Matrix_PDS_Tridiag>(banded_size);
    } else if (
            border_size * n + border_size * (border_size + 1) + (2 * kl + 1 + ku) * banded_size
            >= n * n) {
        return std::make_unique<Matrix_Dense>(n);
    } else {
        block_mat = std::make_unique<Matrix_Banded>(banded_size, kl, ku);
    }
    return std::make_unique<Matrix_Periodic_Banded>(n, kl, ku, std::move(block_mat));
}

std::unique_ptr<Matrix> Matrix::make_new_block_with_banded_region(
        int n,
        int kl,
        int ku,
        bool pds,
        int block1_size,
        int block2_size)
{
    int banded_size(n - block1_size - block2_size);
    std::unique_ptr<Matrix> block_mat;
    if (pds && kl == ku && kl == 1) {
        block_mat = std::make_unique<Matrix_PDS_Tridiag>(banded_size);
    } else if (2 * kl + 1 + ku >= banded_size) {
        return std::make_unique<Matrix_Dense>(n);
    } else {
        block_mat = std::make_unique<Matrix_Banded>(banded_size, kl, ku);
    }
    if (block2_size == 0) {
        return std::make_unique<Matrix_Corner_Block>(n, block1_size, std::move(block_mat));
    } else {
        return std::make_unique<
                Matrix_Center_Block>(n, block1_size, block2_size, std::move(block_mat));
    }
}

std::ostream& operator<<(std::ostream& o, const Matrix& m)
{
    const int n(m.get_size());
    for (int i(0); i < n; ++i) {
        for (int j(0); j < n; ++j) {
            o << std::fixed << std::setprecision(3) << std::setw(10) << m.get_element(i, j);
        }
        o << std::endl;
    }
    return o;
}
