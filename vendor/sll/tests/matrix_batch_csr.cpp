#include <ddc/ddc.hpp>

#include <sll/matrix_batch_csr.hpp>

#include <gtest/gtest.h>

#include "sll/view.hpp"

#include "test_utils.hpp"

Kokkos::View<double**, Kokkos::LayoutRight, Kokkos::DefaultHostExecutionSpace> solve_system(
        Kokkos::View<double**, Kokkos::LayoutRight, Kokkos::DefaultHostExecutionSpace>
                values_view_host,
        Kokkos::View<int*, Kokkos::LayoutRight, Kokkos::DefaultHostExecutionSpace> idx_view_host,
        Kokkos::View<int*, Kokkos::LayoutRight, Kokkos::DefaultHostExecutionSpace>
                nnz_per_row_view_host,
        Kokkos::View<double**, Kokkos::LayoutRight, Kokkos::DefaultHostExecutionSpace> res_host)
{
    int const mat_size = nnz_per_row_view_host.size() - 1;
    int const batch_size = values_view_host.extent(0);
    int const non_zero_per_system = values_view_host.extent(1);
    Kokkos::View<double**, Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace>
            values_view("values", batch_size, non_zero_per_system);
    Kokkos::View<int*, Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace>
            idx_view("col_idxs", non_zero_per_system);
    Kokkos::View<int*, Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace>
            nnz_per_row_view("nnz_per_row", mat_size + 1);
    Kokkos::View<double**, Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace>
            res_view("res", batch_size, mat_size);

    Kokkos::deep_copy(values_view, values_view_host);
    Kokkos::deep_copy(idx_view, idx_view_host);
    Kokkos::deep_copy(nnz_per_row_view, nnz_per_row_view_host);
    Kokkos::deep_copy(res_view, res_host);

    MatrixBatchCsr<Kokkos::DefaultExecutionSpace>
            test_instance(values_view, idx_view, nnz_per_row_view, 100, 1e-12);

    Kokkos::deep_copy(res_view, res_host);
    test_instance.factorize();
    test_instance.solve_inplace(res_view);
    Kokkos::deep_copy(res_host, res_view);
    return res_host;
}

class MatrixBatchCsrFixture : public ::testing::Test
{
protected:
    MatrixBatchCsr<Kokkos::DefaultExecutionSpace> matrix_batch_test;
};

TEST(MatrixBatchCsrFixture, Get_batch_csr)
{
    int const batch_size = 2;
    int const mat_size = 4;
    int const non_zero_per_system = 4;
    double values[] = {2.0, 4.0, 6.0, 8.0, 3.0, 5.0, 7.0, 9.0};
    int col_idxs[] = {0, 1, 2, 3};
    int nnz_per_row[] = {0, 1, 2, 3, 4};
    // rhs and solution
    double res[] = {1., 3.5, 2., 3., 42., 17., 0.5, 1.};
    double solution[]
            = {1. / 2., 3.5 / 4., 2. / 6., 3. / 8., 42. / 3., 17. / 5., 1. / 14., 1. / 9.};

    Kokkos::View<double**, Kokkos::LayoutRight, Kokkos::DefaultHostExecutionSpace>
            values_view_host(values, batch_size, non_zero_per_system);
    Kokkos::View<int*, Kokkos::LayoutRight, Kokkos::DefaultHostExecutionSpace>
            col_idx_view_host(col_idxs, non_zero_per_system);
    Kokkos::View<int*, Kokkos::LayoutRight, Kokkos::DefaultHostExecutionSpace>
            nnz_per_row_view_host(nnz_per_row, mat_size + 1);

    Kokkos::View<double**, Kokkos::LayoutRight, Kokkos::DefaultHostExecutionSpace>
            res_host(res, batch_size, mat_size);

    Kokkos::View<double**, Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace>
            res_view("res", batch_size, mat_size);
    MatrixBatchCsr<Kokkos::DefaultExecutionSpace>
            test_instance(batch_size, mat_size, non_zero_per_system, 100, 1e-6);
    auto [values_view, col_idx_view, nnz_per_row_view] = test_instance.get_batch_csr();

    Kokkos::deep_copy(values_view, values_view_host);
    Kokkos::deep_copy(col_idx_view, col_idx_view_host);
    Kokkos::deep_copy(nnz_per_row_view, nnz_per_row_view_host);
    Kokkos::deep_copy(res_view, res_host);

    test_instance.factorize();
    test_instance.solve_inplace(res_view);

    Kokkos::deep_copy(res_host, res_view);
    ASSERT_EQ(test_instance.norm(0), 8);
    ASSERT_EQ(test_instance.norm(1), 9);
    for (int batch_idx = 0; batch_idx < batch_size; batch_idx++) {
        for (int i = 0; i < mat_size; i++) {
            ASSERT_FLOAT_EQ(res_host(batch_idx, i), solution[batch_idx * mat_size + i]);
        }
    }
}

TEST(MatrixBatchCsrFixture, SolveDiagonal)
{
    int const batch_size = 2;
    int const mat_size = 4;
    int const non_zero_per_system = 4;
    double values[] = {2.0, 4.0, 6.0, 8.0, 3.0, 5.0, 7.0, 9.0};
    int col_idxs[] = {0, 1, 2, 3};
    int nnz_per_row[] = {0, 1, 2, 3, 4};
    // rhs and solution
    double res[] = {0, 3.5, 2., 3., 42., 17., 0.5, 1.};
    double solution[] = {0, 3.5 / 4., 2. / 6., 3. / 8., 42. / 3., 17. / 5., 1. / 14., 1. / 9.};

    Kokkos::View<double**, Kokkos::LayoutRight, Kokkos::DefaultHostExecutionSpace>
            values_view_host(values, batch_size, non_zero_per_system);
    Kokkos::View<int*, Kokkos::LayoutRight, Kokkos::DefaultHostExecutionSpace>
            idx_view_host(col_idxs, non_zero_per_system);
    Kokkos::View<int*, Kokkos::LayoutRight, Kokkos::DefaultHostExecutionSpace>
            nnz_per_row_view_host(nnz_per_row, mat_size + 1);
    Kokkos::View<double**, Kokkos::LayoutRight, Kokkos::DefaultHostExecutionSpace>
            res_view_host(res, batch_size, mat_size);

    solve_system(values_view_host, idx_view_host, nnz_per_row_view_host, res_view_host);

    for (int batch_idx = 0; batch_idx < batch_size; batch_idx++) {
        for (int i = 0; i < mat_size; i++) {
            ASSERT_FLOAT_EQ(res_view_host(batch_idx, i), solution[batch_idx * mat_size + i]);
        }
    }
}

TEST(MatrixBatchCsrFixture, SolveSparse)
{
    int const batch_size = 2;
    int const mat_size = 5;
    int const non_zero_per_system = 16;
    int nnz_per_row[] = {0, 4, 8, 10, 13, 16};

    double solution[]
            = {-107. / 109.,
               114. / 109.,
               -119. / 109.,
               -101. / 109.,
               284. / 109.,
               2346624. / 234305.,
               -138722. / 46861.,
               324305. / 46861.,
               -30223 / 93722.,
               -1519579. / 187444.};
    double matvalues[2][5][5]
            = {{{1., 4., 0, 8, 2},
                {2., 3, 0, 3., 1},
                {0, 2., 1., 0, 0},
                {0, 0, 2., 5., 3},
                {0, 0, 9., 8, 7.}},
               {{11., 3.1, 0, 8.4, 12},
                {1.5, 0.5, 0, 3.7, 1.4},
                {0, 2., 1., 0, 0},
                {0, 0, 2.5, 5.3, 1.8},
                {0, 0, 9.42, 8, 7.6}}};

    Kokkos::View<double**, Kokkos::LayoutRight, Kokkos::DefaultHostExecutionSpace>
            values_view_host("values_host", batch_size, non_zero_per_system);
    Kokkos::View<int*, Kokkos::LayoutRight, Kokkos::DefaultHostExecutionSpace>
            idx_view_host("col_idxs_host", non_zero_per_system);
    Kokkos::View<int*, Kokkos::LayoutRight, Kokkos::DefaultHostExecutionSpace>
            nnz_per_row_view_host(nnz_per_row, mat_size + 1);
    int cpt;
    for (int batch_idx = 0; batch_idx < batch_size; batch_idx++) {
        cpt = 0;
        for (int i = 0; i < mat_size; i++) {
            for (int j = 0; j < mat_size; j++) {
                if (abs(matvalues[batch_idx][i][j]) > 1e-16) {
                    if (batch_idx == 0) {
                        idx_view_host(cpt) = j;
                    }
                    values_view_host(batch_idx, cpt) = matvalues[batch_idx][i][j];
                    cpt++;
                }
            }
        }
    }

    Kokkos::View<double**, Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace>
            values_view("values", batch_size, non_zero_per_system);
    Kokkos::View<int*, Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace>
            idx_view("col_idxs", non_zero_per_system);
    Kokkos::View<int*, Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace>
            nnz_per_row_view("nnz_per_row", mat_size + 1);

    Kokkos::View<double**, Kokkos::LayoutRight, Kokkos::DefaultHostExecutionSpace>
            res_view_host("res_host", batch_size, mat_size);
    Kokkos::View<double**, Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace>
            res_view("res", batch_size, mat_size);
    Kokkos::deep_copy(res_view_host, 1.);

    solve_system(values_view_host, idx_view_host, nnz_per_row_view_host, res_view_host);

    for (int batch_idx = 0; batch_idx < batch_size; batch_idx++) {
        for (int i = 0; i < mat_size; i++) {
            ASSERT_FLOAT_EQ(res_view_host(batch_idx, i), solution[batch_idx * mat_size + i]);
        }
    }
}