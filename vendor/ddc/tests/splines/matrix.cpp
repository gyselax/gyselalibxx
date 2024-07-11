// Copyright (C) The DDC development team, see COPYRIGHT.md file
//
// SPDX-License-Identifier: MIT

#include <algorithm>
#include <cmath>
#include <memory>
#include <sstream>

#include <ddc/ddc.hpp>
#include <ddc/kernels/splines.hpp>

#include <gtest/gtest.h>

#include <Kokkos_DualView.hpp>

#include "test_utils.hpp"

namespace DDC_HIP_5_7_ANONYMOUS_NAMESPACE_WORKAROUND(MATRIX_CPP)
{
    void fill_identity(
            ddc::detail::SplinesLinearProblem<Kokkos::DefaultHostExecutionSpace>::MultiRHS mat)
    {
        assert(mat.extent(0) == mat.extent(1));
        for (std::size_t i(0); i < mat.extent(0); ++i) {
            for (std::size_t j(0); j < mat.extent(1); ++j) {
                mat(i, j) = int(i == j);
            }
        }
    }

    void copy_matrix(
            ddc::detail::SplinesLinearProblem<Kokkos::DefaultHostExecutionSpace>::MultiRHS copy,
            ddc::detail::SplinesLinearProblem<Kokkos::DefaultExecutionSpace> const& mat)
    {
        assert(mat.size() == copy.extent(0));
        assert(mat.size() == copy.extent(1));

        for (std::size_t i(0); i < copy.extent(0); ++i) {
            for (std::size_t j(0); j < copy.extent(1); ++j) {
                copy(i, j) = mat.get_element(i, j);
            }
        }
    }

    void check_inverse(
            ddc::detail::SplinesLinearProblem<Kokkos::DefaultHostExecutionSpace>::MultiRHS matrix,
            ddc::detail::SplinesLinearProblem<Kokkos::DefaultHostExecutionSpace>::MultiRHS inv)
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

    void check_inverse_transpose(
            ddc::detail::SplinesLinearProblem<Kokkos::DefaultHostExecutionSpace>::MultiRHS matrix,
            ddc::detail::SplinesLinearProblem<Kokkos::DefaultHostExecutionSpace>::MultiRHS inv)
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

    void solve_and_validate(
            ddc::detail::SplinesLinearProblem<Kokkos::DefaultExecutionSpace> & matrix)
    {
        const std::size_t N = matrix.size();

        std::vector<double> val_ptr(N * N);
        ddc::detail::SplinesLinearProblem<Kokkos::DefaultHostExecutionSpace>::MultiRHS
                val(val_ptr.data(), N, N);

        copy_matrix(val, matrix);

        matrix.setup_solver();

        Kokkos::DualView<double*> inv_ptr("inv_ptr", N * N);
        ddc::detail::SplinesLinearProblem<Kokkos::DefaultHostExecutionSpace>::MultiRHS
                inv(inv_ptr.h_view.data(), N, N);
        fill_identity(inv);
        inv_ptr.modify_host();
        inv_ptr.sync_device();
        matrix.solve(ddc::detail::SplinesLinearProblem<
                     Kokkos::DefaultExecutionSpace>::MultiRHS(inv_ptr.d_view.data(), N, N));
        inv_ptr.modify_device();
        inv_ptr.sync_host();

        Kokkos::DualView<double*> inv_tr_ptr("inv_tr_ptr", N * N);
        ddc::detail::SplinesLinearProblem<Kokkos::DefaultHostExecutionSpace>::MultiRHS
                inv_tr(inv_tr_ptr.h_view.data(), N, N);
        fill_identity(inv_tr);
        inv_tr_ptr.modify_host();
        inv_tr_ptr.sync_device();
        matrix
                .solve(ddc::detail::SplinesLinearProblem<Kokkos::DefaultExecutionSpace>::
                               MultiRHS(inv_tr_ptr.d_view.data(), N, N),
                       true);
        inv_tr_ptr.modify_device();
        inv_tr_ptr.sync_host();

        check_inverse(val, inv);
        check_inverse_transpose(val, inv_tr);
    }

} // namespace )

TEST(MatrixSparse, Formatting)
{
    ddc::detail::SplinesLinearProblemSparse<Kokkos::DefaultExecutionSpace> matrix(2);
    matrix.set_element(0, 0, 1);
    matrix.set_element(0, 1, 2);
    matrix.set_element(1, 0, 3);
    matrix.set_element(1, 1, 4);
    std::stringstream ss;
    ss << matrix;
    EXPECT_EQ(ss.str(), "     1.000     2.000\n     3.000     4.000\n");
}


TEST(Matrix, Dense)
{
    std::size_t const N = 10;
    std::size_t const k = 10;
    std::unique_ptr<ddc::detail::SplinesLinearProblem<Kokkos::DefaultExecutionSpace>> matrix
            = std::make_unique<
                    ddc::detail::SplinesLinearProblemDense<Kokkos::DefaultExecutionSpace>>(N);

    // Build a non-symmetric full-rank matrix (without zero)
    for (std::size_t i(0); i < N; ++i) {
        matrix->set_element(i, i, 3. / 4 * ((N + 1) * i + 1));
        for (std::size_t j(std::max(0, int(i) - int(k))); j < i; ++j) {
            matrix->set_element(i, j, -(1. / 4) / k * (N * i + j + 1));
        }
        for (std::size_t j(i + 1); j < std::min(N, i + k + 1); ++j) {
            matrix->set_element(i, j, -(1. / 4) / k * (N * i + j + 1));
        }
    }

    solve_and_validate(*matrix);
}

TEST(Matrix, Band)
{
    std::size_t const N = 10;
    std::size_t const k = 3;
    std::unique_ptr<ddc::detail::SplinesLinearProblem<Kokkos::DefaultExecutionSpace>> matrix
            = std::make_unique<
                    ddc::detail::SplinesLinearProblemBand<Kokkos::DefaultExecutionSpace>>(N, k, k);

    // Build a non-symmetric full-rank band matrix
    for (std::size_t i(0); i < N; ++i) {
        matrix->set_element(i, i, 3. / 4 * ((N + 1) * i + 1));
        for (std::size_t j(std::max(0, int(i) - int(k))); j < i; ++j) {
            matrix->set_element(i, j, -(1. / 4) / k * (N * i + j + 1));
        }
        for (std::size_t j(i + 1); j < std::min(N, i + k + 1); ++j) {
            matrix->set_element(i, j, -(1. / 4) / k * (N * i + j + 1));
        }
    }

    solve_and_validate(*matrix);
}

TEST(Matrix, PDSBand)
{
    std::size_t const N = 10;
    std::size_t const k = 3;
    std::unique_ptr<ddc::detail::SplinesLinearProblem<Kokkos::DefaultExecutionSpace>> matrix
            = std::make_unique<
                    ddc::detail::SplinesLinearProblemPDSBand<Kokkos::DefaultExecutionSpace>>(N, k);

    // Build a positive-definite symmetric full-rank band matrix
    for (std::size_t i(0); i < N; ++i) {
        matrix->set_element(i, i, 2.0 * k + 1);
        for (std::size_t j(std::max(0, int(i) - int(k))); j < i; ++j) {
            matrix->set_element(i, j, -1.0);
        }
        for (std::size_t j(i + 1); j < std::min(N, i + k + 1); ++j) {
            matrix->set_element(i, j, -1.0);
        }
    }

    solve_and_validate(*matrix);
}

TEST(Matrix, PDSTridiag)
{
    std::size_t const N = 10;
    std::size_t const k = 1;
    std::unique_ptr<ddc::detail::SplinesLinearProblem<Kokkos::DefaultExecutionSpace>> matrix
            = std::make_unique<
                    ddc::detail::SplinesLinearProblemPDSTridiag<Kokkos::DefaultExecutionSpace>>(N);

    // Build a positive-definite symmetric full-rank tridiagonal matrix
    for (std::size_t i(0); i < N; ++i) {
        matrix->set_element(i, i, 2.0 * k + 1);
        for (std::size_t j(std::max(0, int(i) - int(k))); j < i; ++j) {
            matrix->set_element(i, j, -1.0);
        }
        for (std::size_t j(i + 1); j < std::min(N, i + k + 1); ++j) {
            matrix->set_element(i, j, -1.0);
        }
    }

    solve_and_validate(*matrix);
}

TEST(Matrix, 2x2Blocks)
{
    std::size_t const N = 10;
    std::size_t const k = 10;
    std::unique_ptr<ddc::detail::SplinesLinearProblem<Kokkos::DefaultExecutionSpace>> top_left_block
            = std::make_unique<
                    ddc::detail::SplinesLinearProblemDense<Kokkos::DefaultExecutionSpace>>(7);
    std::unique_ptr<ddc::detail::SplinesLinearProblem<Kokkos::DefaultExecutionSpace>> matrix
            = std::make_unique<ddc::detail::SplinesLinearProblem2x2Blocks<
                    Kokkos::DefaultExecutionSpace>>(N, std::move(top_left_block));

    // Build a non-symmetric full-rank matrix (without zero)
    for (std::size_t i(0); i < N; ++i) {
        matrix->set_element(i, i, 3. / 4 * ((N + 1) * i + 1));
        for (std::size_t j(std::max(0, int(i) - int(k))); j < i; ++j) {
            matrix->set_element(i, j, -(1. / 4) / k * (N * i + j + 1));
        }
        for (std::size_t j(i + 1); j < std::min(N, i + k + 1); ++j) {
            matrix->set_element(i, j, -(1. / 4) / k * (N * i + j + 1));
        }
    }

    solve_and_validate(*matrix);
}

TEST(Matrix, 3x3Blocks)
{
    std::size_t const N = 10;
    std::size_t const k = 10;
    std::unique_ptr<ddc::detail::SplinesLinearProblem<Kokkos::DefaultExecutionSpace>> center_block
            = std::make_unique<
                    ddc::detail::SplinesLinearProblemDense<Kokkos::DefaultExecutionSpace>>(N - 5);
    std::unique_ptr<ddc::detail::SplinesLinearProblem<Kokkos::DefaultExecutionSpace>> matrix
            = std::make_unique<ddc::detail::SplinesLinearProblem3x3Blocks<
                    Kokkos::DefaultExecutionSpace>>(N, 2, std::move(center_block));

    // Build a non-symmetric full-rank matrix (without zero)
    for (std::size_t i(0); i < N; ++i) {
        std::cout << i;
        matrix->set_element(i, i, 3. / 4 * ((N + 1) * i + 1));
        for (std::size_t j(std::max(0, int(i) - int(k))); j < i; ++j) {
            matrix->set_element(i, j, -(1. / 4) / k * (N * i + j + 1));
        }
        for (std::size_t j(i + 1); j < std::min(N, i + k + 1); ++j) {
            matrix->set_element(i, j, -(1. / 4) / k * (N * i + j + 1));
        }
    }

    solve_and_validate(*matrix);
}

class MatrixSizesFixture : public testing::TestWithParam<std::tuple<std::size_t, std::size_t>>
{
};

TEST_P(MatrixSizesFixture, NonSymmetric)
{
    auto const [N, k] = GetParam();
    std::unique_ptr<ddc::detail::SplinesLinearProblem<Kokkos::DefaultExecutionSpace>> matrix
            = ddc::detail::SplinesLinearProblemMaker::make_new_band<
                    Kokkos::DefaultExecutionSpace>(N, k, k, false);

    // Build a non-symmetric full-rank band matrix
    for (std::size_t i(0); i < N; ++i) {
        matrix->set_element(i, i, 3. / 4 * ((N + 1) * i + 1));
        for (std::size_t j(std::max(0, int(i) - int(k))); j < i; ++j) {
            matrix->set_element(i, j, -(1. / 4) / k * (N * i + j + 1));
        }
        for (std::size_t j(i + 1); j < std::min(N, i + k + 1); ++j) {
            matrix->set_element(i, j, -(1. / 4) / k * (N * i + j + 1));
        }
    }

    solve_and_validate(*matrix);
}

TEST_P(MatrixSizesFixture, PositiveDefiniteSymmetric)
{
    auto const [N, k] = GetParam();
    std::unique_ptr<ddc::detail::SplinesLinearProblem<Kokkos::DefaultExecutionSpace>> matrix
            = ddc::detail::SplinesLinearProblemMaker::make_new_band<
                    Kokkos::DefaultExecutionSpace>(N, k, k, true);

    // Build a positive-definite symmetric full-rank band matrix
    for (std::size_t i(0); i < N; ++i) {
        matrix->set_element(i, i, 2.0 * k + 1);
        for (std::size_t j(std::max(0, int(i) - int(k))); j < i; ++j) {
            matrix->set_element(i, j, -1.0);
        }
        for (std::size_t j(i + 1); j < std::min(N, i + k + 1); ++j) {
            matrix->set_element(i, j, -1.0);
        }
    }

    solve_and_validate(*matrix);
}

TEST_P(MatrixSizesFixture, OffsetBanded)
{
    auto const [N, k] = GetParam();
    std::unique_ptr<ddc::detail::SplinesLinearProblem<Kokkos::DefaultExecutionSpace>> matrix
            = ddc::detail::SplinesLinearProblemMaker::make_new_band<
                    Kokkos::DefaultExecutionSpace>(N, 0, 2 * k, true);

    // Build a positive-definite symmetric full-rank band matrix permuted in such a way the band is shifted
    for (std::size_t i(0); i < N; ++i) {
        for (std::size_t j(i); j < std::min(N, i + k); ++j) {
            matrix->set_element(i, i, -1.0);
        }
        if (i + k < N) {
            matrix->set_element(i, i + k, 2.0 * k + 1);
        }
        for (std::size_t j(i + k + 1); j < std::min(N, i + k + 1); ++j) {
            matrix->set_element(i, j, -1.0);
        }
    }

    solve_and_validate(*matrix);
}

TEST_P(MatrixSizesFixture, 2x2Blocks)
{
    auto const [N, k] = GetParam();
    std::unique_ptr<ddc::detail::SplinesLinearProblem<Kokkos::DefaultExecutionSpace>> matrix
            = ddc::detail::SplinesLinearProblemMaker::make_new_block_matrix_with_band_main_block<
                    Kokkos::DefaultExecutionSpace>(N, k, k, false, 3);

    // Build a non-symmetric full-rank band matrix
    for (std::size_t i(0); i < N; ++i) {
        matrix->set_element(i, i, 3. / 4 * ((N + 1) * i + 1));
        for (std::size_t j(std::max(0, int(i) - int(k))); j < i; ++j) {
            matrix->set_element(i, j, -(1. / 4) / k * (N * i + j + 1));
        }
        for (std::size_t j(i + 1); j < std::min(N, i + k + 1); ++j) {
            matrix->set_element(i, j, -(1. / 4) / k * (N * i + j + 1));
        }
    }

    solve_and_validate(*matrix);
}

TEST_P(MatrixSizesFixture, 3x3Blocks)
{
    auto const [N, k] = GetParam();
    std::unique_ptr<ddc::detail::SplinesLinearProblem<Kokkos::DefaultExecutionSpace>> matrix
            = ddc::detail::SplinesLinearProblemMaker::make_new_block_matrix_with_band_main_block<
                    Kokkos::DefaultExecutionSpace>(N, k, k, false, 3, 2);

    // Build a non-symmetric full-rank band matrix
    for (std::size_t i(0); i < N; ++i) {
        matrix->set_element(i, i, 3. / 4 * ((N + 1) * i + 1));
        for (std::size_t j(std::max(0, int(i) - int(k))); j < i; ++j) {
            matrix->set_element(i, j, -(1. / 4) / k * (N * i + j + 1));
        }
        for (std::size_t j(i + 1); j < std::min(N, i + k + 1); ++j) {
            matrix->set_element(i, j, -(1. / 4) / k * (N * i + j + 1));
        }
    }

    solve_and_validate(*matrix);
}

TEST_P(MatrixSizesFixture, PeriodicBand)
{
    auto const [N, k] = GetParam();

    // Build a full-rank periodic band matrix permuted in such a way the band is shifted
    for (std::ptrdiff_t s(-k + k / 2 + 1); s < static_cast<std::ptrdiff_t>(k - k / 2); ++s) {
        std::unique_ptr<ddc::detail::SplinesLinearProblem<Kokkos::DefaultExecutionSpace>> matrix
                = ddc::detail::SplinesLinearProblemMaker::make_new_periodic_band_matrix<
                        Kokkos::DefaultExecutionSpace>(
                        N,
                        static_cast<std::ptrdiff_t>(k - s),
                        k + s,
                        false);
        for (std::size_t i(0); i < N; ++i) {
            for (std::size_t j(0); j < N; ++j) {
                std::ptrdiff_t diag = ddc::detail::
                        modulo(static_cast<std::ptrdiff_t>(j - i), static_cast<std::ptrdiff_t>(N));
                if (diag == s || diag == N + s) {
                    matrix->set_element(i, j, 2.0 * k + 1);
                } else if (diag <= s + k || diag >= N + s - k) {
                    matrix->set_element(i, j, -1.);
                }
            }
        }

        solve_and_validate(*matrix);
    }
}

TEST_P(MatrixSizesFixture, Sparse)
{
    auto const [N, k] = GetParam();
    std::unique_ptr<ddc::detail::SplinesLinearProblem<Kokkos::DefaultExecutionSpace>> matrix
            = ddc::detail::SplinesLinearProblemMaker::make_new_sparse<
                    Kokkos::DefaultExecutionSpace>(N);

    // Build a positive-definite symmetric diagonal-dominant band matrix (stored as a sparse matrix)
    for (std::size_t i(0); i < N; ++i) {
        for (std::size_t j(0); j < N; ++j) {
            if (i == j) {
                matrix->set_element(i, j, 3. / 4);
            } else if (
                    std::abs(static_cast<std::ptrdiff_t>(j - i))
                    <= static_cast<std::ptrdiff_t>(k)) {
                matrix->set_element(i, j, -(1. / 4) / k);
            }
        }
    }

    solve_and_validate(*matrix);
}

INSTANTIATE_TEST_SUITE_P(
        MyGroup,
        MatrixSizesFixture,
        testing::Combine(testing::Values<std::size_t>(10, 20), testing::Range<std::size_t>(1, 7)));
