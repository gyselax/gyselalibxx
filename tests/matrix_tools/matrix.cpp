#include <algorithm>
#include <cmath>
#include <memory>

#include <gtest/gtest.h>

#include "math_tools.hpp"
#include "matrix.hpp"
#include "view.hpp"

namespace {

void fill_identity(DSpan2D mat)
{
    assert(mat.extent(0) == mat.extent(1));
    for (std::size_t i(0); i < mat.extent(0); ++i) {
        for (std::size_t j(0); j < mat.extent(1); ++j) {
            mat(i, j) = int(i == j);
        }
    }
}

void copy_matrix(DSpan2D copy, std::unique_ptr<Matrix>& mat)
{
    assert(mat->get_size() == int(copy.extent(0)));
    assert(mat->get_size() == int(copy.extent(1)));

    for (std::size_t i(0); i < copy.extent(0); ++i) {
        for (std::size_t j(0); j < copy.extent(1); ++j) {
            copy(i, j) = mat->get_element(i, j);
        }
    }
}

void check_inverse(DSpan2D matrix, DSpan2D inv)
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

void check_inverse_transpose(DSpan2D matrix, DSpan2D inv)
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

TEST_P(MatrixSizesFixture, PositiveDefiniteSymmetric)
{
    auto const [N, k] = GetParam();
    std::unique_ptr<Matrix> matrix = Matrix::make_new_banded(N, k, k, true);

    for (std::size_t i(0); i < N; ++i) {
        matrix->set_element(i, i, 2.0 * k);
        for (std::size_t j(std::max(0, int(i) - int(k))); j < i; ++j) {
            matrix->set_element(i, j, -1.0);
        }
        for (std::size_t j(i + 1); j < std::min(N, i + k + 1); ++j) {
            matrix->set_element(i, j, -1.0);
        }
    }
    std::vector<double> val_ptr(N * N);
    DSpan2D val(val_ptr.data(), N, N);
    copy_matrix(val, matrix);

    std::vector<double> inv_ptr(N * N);
    DSpan2D inv(inv_ptr.data(), N, N);
    fill_identity(inv);
    matrix->factorise();
    matrix->solve_multiple_inplace(inv);
    check_inverse(val, inv);
}

TEST_P(MatrixSizesFixture, OffsetBanded)
{
    auto const [N, k] = GetParam();
    std::unique_ptr<Matrix> matrix = Matrix::make_new_banded(N, 0, 2 * k, true);

    for (std::size_t i(0); i < N; ++i) {
        for (std::size_t j(i); j < std::min(N, i + k); ++j) {
            matrix->set_element(i, i, -1.0);
        }
        if (i + k < N) {
            matrix->set_element(i, i + k, 2.0 * k);
        }
        for (std::size_t j(i + k + 1); j < std::min(N, i + k + 1); ++j) {
            matrix->set_element(i, j, -1.0);
        }
    }
    std::vector<double> val_ptr(N * N);
    DSpan2D val(val_ptr.data(), N, N);
    copy_matrix(val, matrix);

    std::vector<double> inv_ptr(N * N);
    DSpan2D inv(inv_ptr.data(), N, N);
    fill_identity(inv);
    matrix->factorise();
    matrix->solve_multiple_inplace(inv);
    check_inverse(val, inv);
}

TEST_P(MatrixSizesFixture, PeriodicBanded)
{
    auto const [N, k] = GetParam();

    for (int s(-k); s < k + 1; ++s) {
        if (s == 0)
            continue;

        std::unique_ptr<Matrix> matrix = Matrix::make_new_periodic_banded(N, k - s, k + s, false);
        for (int i(0); i < N; ++i) {
            for (int j(0); j < N; ++j) {
                int diag = modulo(j - i, int(N));
                if (diag == s || diag == N + s) {
                    matrix->set_element(i, j, 0.5);
                } else if (diag <= s + k || diag >= N + s - k) {
                    matrix->set_element(i, j, -1.0 / k);
                }
            }
        }
        std::vector<double> val_ptr(N * N);
        DSpan2D val(val_ptr.data(), N, N);
        copy_matrix(val, matrix);

        std::vector<double> inv_ptr(N * N);
        DSpan2D inv(inv_ptr.data(), N, N);
        fill_identity(inv);
        matrix->factorise();
        matrix->solve_multiple_inplace(inv);
        check_inverse(val, inv);
    }
}

TEST_P(MatrixSizesFixture, PositiveDefiniteSymmetricTranspose)
{
    auto const [N, k] = GetParam();
    std::unique_ptr<Matrix> matrix = Matrix::make_new_banded(N, k, k, true);

    for (std::size_t i(0); i < N; ++i) {
        matrix->set_element(i, i, 2.0 * k);
        for (std::size_t j(std::max(0, int(i) - int(k))); j < i; ++j) {
            matrix->set_element(i, j, -1.0);
        }
        for (std::size_t j(i + 1); j < std::min(N, i + k + 1); ++j) {
            matrix->set_element(i, j, -1.0);
        }
    }
    std::vector<double> val_ptr(N * N);
    DSpan2D val(val_ptr.data(), N, N);
    copy_matrix(val, matrix);

    std::vector<double> inv_ptr(N * N);
    DSpan2D inv(inv_ptr.data(), N, N);
    fill_identity(inv);
    matrix->factorise();
    for (std::size_t i(0); i < N; ++i) {
        DSpan1D inv_line(inv_ptr.data() + i * N, N);
        matrix->solve_transpose_inplace(inv_line);
    }
    check_inverse_transpose(val, inv);
}

TEST_P(MatrixSizesFixture, OffsetBandedTranspose)
{
    auto const [N, k] = GetParam();
    std::unique_ptr<Matrix> matrix = Matrix::make_new_banded(N, 0, 2 * k, true);

    for (std::size_t i(0); i < N; ++i) {
        for (std::size_t j(i); j < std::min(N, i + k); ++j) {
            matrix->set_element(i, i, -1.0);
        }
        if (i + k < N) {
            matrix->set_element(i, i + k, 2.0 * k);
        }
        for (std::size_t j(i + k + 1); j < std::min(N, i + k + 1); ++j) {
            matrix->set_element(i, j, -1.0);
        }
    }
    std::vector<double> val_ptr(N * N);
    DSpan2D val(val_ptr.data(), N, N);
    copy_matrix(val, matrix);

    std::vector<double> inv_ptr(N * N);
    DSpan2D inv(inv_ptr.data(), N, N);
    fill_identity(inv);
    matrix->factorise();
    for (std::size_t i(0); i < N; ++i) {
        DSpan1D inv_line(inv_ptr.data() + i * N, N);
        matrix->solve_transpose_inplace(inv_line);
    }
    check_inverse_transpose(val, inv);
}

TEST_P(MatrixSizesFixture, PeriodicBandedTranspose)
{
    auto const [N, k] = GetParam();

    for (int s(-k); s < k + 1; ++s) {
        if (s == 0)
            continue;

        std::unique_ptr<Matrix> matrix = Matrix::make_new_periodic_banded(N, k - s, k + s, false);
        for (int i(0); i < N; ++i) {
            for (int j(0); j < N; ++j) {
                int diag = modulo(j - i, int(N));
                if (diag == s || diag == N + s) {
                    matrix->set_element(i, j, 0.5);
                } else if (diag <= s + k || diag >= N + s - k) {
                    matrix->set_element(i, j, -1.0 / k);
                }
            }
        }
        std::vector<double> val_ptr(N * N);
        DSpan2D val(val_ptr.data(), N, N);
        copy_matrix(val, matrix);

        std::vector<double> inv_ptr(N * N);
        DSpan2D inv(inv_ptr.data(), N, N);
        fill_identity(inv);
        matrix->factorise();
        for (int i(0); i < N; ++i) {
            DSpan1D inv_line(inv_ptr.data() + i * N, N);
            matrix->solve_transpose_inplace(inv_line);
        }
        check_inverse_transpose(val, inv);
    }
}

INSTANTIATE_TEST_SUITE_P(
        MyGroup,
        MatrixSizesFixture,
        testing::Combine(testing::Values<std::size_t>(10, 20), testing::Range<std::size_t>(1, 7)));
