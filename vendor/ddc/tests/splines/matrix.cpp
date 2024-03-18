// Copyright (C) The DDC development team, see COPYRIGHT.md file
//
// SPDX-License-Identifier: MIT

#include <algorithm>
#include <cmath>
#include <memory>

#include <ddc/ddc.hpp>
#include <ddc/kernels/splines.hpp>

#include <gtest/gtest.h>

#include <Kokkos_DualView.hpp>

#include "ddc/kernels/splines/view.hpp"

#include "test_utils.hpp"

namespace {

void fill_identity(ddc::DSpan2D mat)
{
    assert(mat.extent(0) == mat.extent(1));
    for (std::size_t i(0); i < mat.extent(0); ++i) {
        for (std::size_t j(0); j < mat.extent(1); ++j) {
            mat(i, j) = int(i == j);
        }
    }
}

/*
void copy_matrix(ddc::DSpan2D copy, std::unique_ptr<ddc::detail::Matrix>& mat)
{
    assert(mat->get_size() == int(copy.extent(0)));
    assert(mat->get_size() == int(copy.extent(1)));

    for (std::size_t i(0); i < copy.extent(0); ++i) {
        for (std::size_t j(0); j < copy.extent(1); ++j) {
            copy(i, j) = mat->get_element(i, j);
        }
    }
}
*/

void check_inverse(ddc::DSpan2D matrix, ddc::DSpan2D inv)
{
    double TOL = 1e-10;
    std::size_t N = matrix.extent(0);

    for (std::size_t i(0); i < N; ++i) {
        for (std::size_t j(0); j < N; ++j) {
            double id_val = 0.0;
            for (std::size_t k(0); k < N; ++k) {
                id_val += matrix(i, k) * inv(j, k);
            }
            EXPECT_NEAR(id_val, static_cast<double>(i == j), TOL);
        }
    }
}

void check_inverse_transpose(ddc::DSpan2D matrix, ddc::DSpan2D inv)
{
    double TOL = 1e-10;
    std::size_t N = matrix.extent(0);

    for (std::size_t i(0); i < N; ++i) {
        for (std::size_t j(0); j < N; ++j) {
            double id_val = 0.0;
            for (std::size_t k(0); k < N; ++k) {
                id_val += matrix(i, k) * inv(k, j);
            }
            EXPECT_NEAR(id_val, static_cast<double>(i == j), TOL);
        }
    }
}
} // namespace

class MatrixSizesFixture : public testing::TestWithParam<std::tuple<std::size_t, std::size_t>>
{
};

TEST_P(MatrixSizesFixture, Sparse)
{
    auto const [N, k] = GetParam();


    std::unique_ptr<ddc::detail::Matrix> matrix
            = ddc::detail::MatrixMaker::make_new_sparse<Kokkos::DefaultExecutionSpace>(N);

    std::vector<double> val_ptr(N * N);
    ddc::DSpan2D val(val_ptr.data(), N, N);
    for (std::size_t i(0); i < N; ++i) {
        for (std::size_t j(0); j < N; ++j) {
            if (i == j) {
                matrix->set_element(i, j, 3. / 4);
                val(i, j) = 3. / 4;
            } else if (std::abs((std::ptrdiff_t)(j - i)) <= (std::ptrdiff_t)k) {
                matrix->set_element(i, j, -(1. / 4) / k);
                val(i, j) = -(1. / 4) / k;
            } else {
                val(i, j) = 0.;
            }
        }
    }
    // copy_matrix(val, matrix); // copy_matrix is not available for sparse matrix because of a limitation of Ginkgo API (get_element is not implemented). The workaround is to fill val directly in the loop

    matrix->factorize();

    Kokkos::DualView<double*> inv_ptr("inv_ptr", N * N);
    ddc::DSpan2D inv(inv_ptr.h_view.data(), N, N);
    fill_identity(inv);
    inv_ptr.modify_host();
    inv_ptr.sync_device();
    matrix->solve_multiple_inplace(ddc::DSpan2D(inv_ptr.d_view.data(), N, N));
    inv_ptr.modify_device();
    inv_ptr.sync_host();

    Kokkos::DualView<double*> inv_tr_ptr("inv_tr_ptr", N * N);
    ddc::DSpan2D inv_tr(inv_tr_ptr.h_view.data(), N, N);
    fill_identity(inv_tr);
    inv_tr_ptr.modify_host();
    inv_tr_ptr.sync_device();
    matrix->solve_multiple_transpose_inplace(ddc::DSpan2D(inv_tr_ptr.d_view.data(), N, N));
    inv_tr_ptr.modify_device();
    inv_tr_ptr.sync_host();

    check_inverse(val, inv);
    check_inverse_transpose(val, inv_tr);
}

INSTANTIATE_TEST_SUITE_P(
        MyGroup,
        MatrixSizesFixture,
        testing::Combine(testing::Values<std::size_t>(10, 20), testing::Range<std::size_t>(1, 7)));
