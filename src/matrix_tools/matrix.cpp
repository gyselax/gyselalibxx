#include <algorithm>
#include <cassert>
#include <iomanip>
#include <iostream>
#include <memory>
#include <utility>

#include "matrix.hpp"
#include "matrix_banded.hpp"
#include "matrix_centre_block.hpp"
#include "matrix_corner_block.hpp"
#include "matrix_dense.hpp"
#include "matrix_pds_tridiag.hpp"
#include "matrix_periodic_banded.hpp"
#include "view.hpp"

using std::max;
using std::min;

DSpan1D Matrix::solve_inplace(DSpan1D const b) const
{
    assert(int(b.extent(0)) == n);
    int const info = solve_inplace_method(b.data_handle(), 'N', 1);

    if (info < 0) {
        std::cerr << -info << "-th argument had an illegal value" << std::endl;
        // TODO: Add LOG_FATAL_ERROR
    }
    return b;
}

DSpan1D Matrix::solve_transpose_inplace(DSpan1D const b) const
{
    assert(int(b.extent(0)) == n);
    int const info = solve_inplace_method(b.data_handle(), 'T', 1);

    if (info < 0) {
        std::cerr << -info << "-th argument had an illegal value" << std::endl;
        // TODO: Add LOG_FATAL_ERROR
    }
    return b;
}

DSpan2D Matrix::solve_multiple_inplace(DSpan2D const bx) const
{
    assert(int(bx.extent(1)) == n);
    int const info = solve_inplace_method(bx.data_handle(), 'N', bx.extent(0));

    if (info < 0) {
        std::cerr << -info << "-th argument had an illegal value" << std::endl;
        // TODO: Add LOG_FATAL_ERROR
    }
    return bx;
}

void Matrix::factorise()
{
    int const info = factorise_method();

    if (info < 0) {
        std::cerr << -info << "-th argument had an illegal value" << std::endl;
        // TODO: Add LOG_FATAL_ERROR
    } else if (info > 0) {
        std::cerr << "U(" << info << "," << info << ") is exactly zero.";
        std::cerr << " The factorisation has been completed, but the factor";
        std::cerr << " U is exactly singular, and division by zero will occur "
                     "if "
                     "it is used to";
        std::cerr << " solve a system of equations.";
        // TODO: Add LOG_FATAL_ERROR
    }
}

std::unique_ptr<Matrix> Matrix::make_new_banded(
        int const n,
        int const kl,
        int const ku,
        bool const pds)
{
    if (kl == ku && kl == 1 && pds) {
        return std::make_unique<Matrix_PDS_Tridiag>(n);
    } else if (2 * kl + 1 + ku >= n) {
        return std::make_unique<Matrix_Dense>(n);
    } else {
        return std::make_unique<Matrix_Banded>(n, kl, ku);
    }
}

std::unique_ptr<Matrix> Matrix::make_new_periodic_banded(
        int const n,
        int const kl,
        int const ku,
        bool const pds)
{
    int const border_size = max(kl, ku);
    int const banded_size = n - border_size;
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
        int const n,
        int const kl,
        int const ku,
        bool const pds,
        int const block1_size,
        int const block2_size)
{
    int const banded_size = n - block1_size - block2_size;
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
                Matrix_Centre_Block>(n, block1_size, block2_size, std::move(block_mat));
    }
}

std::ostream& operator<<(std::ostream& os, Matrix const& m)
{
    int const n = m.get_size();
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            os << std::fixed << std::setprecision(3) << std::setw(10) << m.get_element(i, j);
        }
        os << std::endl;
    }
    return os;
}
