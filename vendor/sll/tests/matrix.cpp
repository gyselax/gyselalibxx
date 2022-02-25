#include <algorithm>
#include <cmath>
#include <memory>

#include <sll/math_tools.hpp>
#include <sll/matrix.hpp>

#include <gtest/gtest.h>

#include "sll/view.hpp"

#include "test_utils.hpp"

template <class T>
struct MatrixSizesFixture;

template <std::size_t N, std::size_t k>
struct MatrixSizesFixture<
        std::tuple<std::integral_constant<std::size_t, N>, std::integral_constant<std::size_t, k>>>
    : public testing::Test
{
    static constexpr std::size_t matrix_size = N;
    static constexpr std::size_t ndiags = k;
};

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
    assert(mat->get_size() == copy.extent(0));
    assert(mat->get_size() == copy.extent(1));

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
            if (i == j) {
                ASSERT_TRUE(fabs(id_val - 1.0) < TOL);
            } else {
                ASSERT_TRUE(fabs(id_val) < TOL);
            }
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
            if (i == j) {
                ASSERT_TRUE(fabs(id_val - 1.0) < TOL);
            } else {
                ASSERT_TRUE(fabs(id_val) < TOL);
            }
        }
    }
}

using sizes = std::integer_sequence<std::size_t, 10, 20>;
using diagonals = std::integer_sequence<std::size_t, 1, 2, 3, 4, 5, 6>;

using Cases = tuple_to_types_t<cartesian_product_t<sizes, diagonals>>;

TYPED_TEST_SUITE(MatrixSizesFixture, Cases);

TYPED_TEST(MatrixSizesFixture, PositiveDefiniteSymmetric)
{
    constexpr std::size_t N = TestFixture::matrix_size;
    constexpr std::size_t k = TestFixture::ndiags;
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
    double val_ptr[N * N];
    DSpan2D val(val_ptr, N, N);
    copy_matrix(val, matrix);

    double inv_ptr[N * N];
    DSpan2D inv(inv_ptr, N, N);
    fill_identity(inv);
    matrix->factorize();
    matrix->solve_multiple_inplace(inv);
    check_inverse(val, inv);
}

TYPED_TEST(MatrixSizesFixture, OffsetBanded)
{
    constexpr std::size_t N = TestFixture::matrix_size;
    constexpr std::size_t k = TestFixture::ndiags;
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
    double val_ptr[N * N];
    DSpan2D val(val_ptr, N, N);
    copy_matrix(val, matrix);

    double inv_ptr[N * N];
    DSpan2D inv(inv_ptr, N, N);
    fill_identity(inv);
    matrix->factorize();
    matrix->solve_multiple_inplace(inv);
    check_inverse(val, inv);
}

TYPED_TEST(MatrixSizesFixture, PeriodicBanded)
{
    constexpr int N = TestFixture::matrix_size;
    constexpr int k = TestFixture::ndiags;

    for (int s(-k); s < k + 1; ++s) {
        if (s == 0)
            continue;

        std::unique_ptr<Matrix> matrix = Matrix::make_new_periodic_banded(N, k - s, k + s, false);
        for (int i(0); i < N; ++i) {
            for (int j(0); j < N; ++j) {
                int diag = modulo(j - i, int(N));
                if (diag == s or diag == N + s) {
                    matrix->set_element(i, j, 0.5);
                } else if (diag <= s + k or diag >= N + s - k) {
                    matrix->set_element(i, j, -1.0 / k);
                }
            }
        }
        double val_ptr[N * N];
        DSpan2D val(val_ptr, N, N);
        copy_matrix(val, matrix);

        double inv_ptr[N * N];
        DSpan2D inv(inv_ptr, N, N);
        fill_identity(inv);
        matrix->factorize();
        matrix->solve_multiple_inplace(inv);
        check_inverse(val, inv);
    }
}

TYPED_TEST(MatrixSizesFixture, PositiveDefiniteSymmetricTranspose)
{
    constexpr std::size_t N = TestFixture::matrix_size;
    constexpr std::size_t k = TestFixture::ndiags;
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
    double val_ptr[N * N];
    DSpan2D val(val_ptr, N, N);
    copy_matrix(val, matrix);

    double inv_ptr[N * N];
    DSpan2D inv(inv_ptr, N, N);
    fill_identity(inv);
    matrix->factorize();
    for (int i(0); i < N; ++i) {
        DSpan1D inv_line(inv_ptr + i * N, N);
        matrix->solve_transpose_inplace(inv_line);
    }
    check_inverse_transpose(val, inv);
}

TYPED_TEST(MatrixSizesFixture, OffsetBandedTranspose)
{
    constexpr std::size_t N = TestFixture::matrix_size;
    constexpr std::size_t k = TestFixture::ndiags;
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
    double val_ptr[N * N];
    DSpan2D val(val_ptr, N, N);
    copy_matrix(val, matrix);

    double inv_ptr[N * N];
    DSpan2D inv(inv_ptr, N, N);
    fill_identity(inv);
    matrix->factorize();
    for (int i(0); i < N; ++i) {
        DSpan1D inv_line(inv_ptr + i * N, N);
        matrix->solve_transpose_inplace(inv_line);
    }
    check_inverse_transpose(val, inv);
}

TYPED_TEST(MatrixSizesFixture, PeriodicBandedTranspose)
{
    constexpr int N = TestFixture::matrix_size;
    constexpr int k = TestFixture::ndiags;

    for (int s(-k); s < k + 1; ++s) {
        if (s == 0)
            continue;

        std::unique_ptr<Matrix> matrix = Matrix::make_new_periodic_banded(N, k - s, k + s, false);
        for (int i(0); i < N; ++i) {
            for (int j(0); j < N; ++j) {
                int diag = modulo(j - i, int(N));
                if (diag == s or diag == N + s) {
                    matrix->set_element(i, j, 0.5);
                } else if (diag <= s + k or diag >= N + s - k) {
                    matrix->set_element(i, j, -1.0 / k);
                }
            }
        }
        double val_ptr[N * N];
        DSpan2D val(val_ptr, N, N);
        copy_matrix(val, matrix);

        double inv_ptr[N * N];
        DSpan2D inv(inv_ptr, N, N);
        fill_identity(inv);
        matrix->factorize();
        for (int i(0); i < N; ++i) {
            DSpan1D inv_line(inv_ptr + i * N, N);
            matrix->solve_transpose_inplace(inv_line);
        }
        check_inverse_transpose(val, inv);
    }
}
