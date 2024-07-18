#include <ddc/ddc.hpp>

#include <sll/matrix_batch_ell.hpp>

#include <gtest/gtest.h>

#include "sll/view.hpp"

#include "test_utils.hpp"



void fill_values(Kokkos::View<int**, Kokkos::LayoutLeft, Kokkos::DefaultExecutionSpace> idx_view)
{
    Kokkos::parallel_for(
            "idx_inner_loop",
            Kokkos::RangePolicy<Kokkos::DefaultExecutionSpace>(0, 4),
            KOKKOS_LAMBDA(const int i) { idx_view(i, 0) = i; });
}



class MatrixBatchEllFixture : public ::testing::Test
{
protected:
    MatrixBatchEll<Kokkos::DefaultExecutionSpace> matrix_batch_test;
};

TEST(MatrixBatchEllFixture, Init)
{
    int const batch_size = 2;
    int const mat_size = 4;
    int const non_zero_per_col = 1;
    double values[] = {2.0, 2.0, 2.0, 2.0, 3.0, 3.0, 3.0, 3.0};
    int col_idxs[] = {0, 1, 2, 3};

    Kokkos::LayoutStride values_layout(
            batch_size,
            non_zero_per_col * mat_size,
            mat_size,
            1,
            non_zero_per_col,
            mat_size);
    Kokkos::View<double***, Kokkos::LayoutStride, Kokkos::DefaultExecutionSpace>
            values_view(values, values_layout);
    Kokkos::View<int**, Kokkos::LayoutLeft, Kokkos::DefaultExecutionSpace>
            idx_view(col_idxs, mat_size, non_zero_per_col);
    MatrixBatchEll<Kokkos::DefaultExecutionSpace> test_instance(idx_view, values_view, 1000, 0.001);
    ASSERT_EQ(test_instance.get_batch_size(), batch_size);
    ASSERT_EQ(test_instance.get_size(), mat_size);
}


TEST(MatrixBatchEllFixture, SetGetElement)
{
    int const batch_size = 2;
    int const mat_size = 4;
    int const non_zero_per_col = 1;
    double values[] = {2.0, 2.0, 2.0, 2.0, 3.0, 3.0, 3.0, 3.0};
    int col_idx[] = {0, 1, 2, 3};
    Kokkos::LayoutStride values_layout(
            batch_size,
            non_zero_per_col * mat_size,
            mat_size,
            1,
            non_zero_per_col,
            mat_size);
    Kokkos::View<double***, Kokkos::LayoutStride, Kokkos::DefaultExecutionSpace>
            values_view(values, values_layout);
    Kokkos::View<int**, Kokkos::LayoutLeft, Kokkos::DefaultExecutionSpace>
            idx_view(col_idx, mat_size, non_zero_per_col);
    MatrixBatchEll<Kokkos::DefaultExecutionSpace> test_instance(idx_view, values_view, 1000, 1e-8);
    auto [idx, vals] = test_instance.get_batch_idx_and_vals();

    test_instance.set_ell_element(0, 3, 0, 42);
    ASSERT_EQ(test_instance.get_ell_element(0, 3, 0), 42);
}


TEST(MatrixBatchEllFixture, GetIdxAndVals)
{
    int const batch_size = 2;
    int const mat_size = 4;
    int const non_zero_per_col = 1;
    double values[] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0};
    int col_idx[] = {0, 1, 2, 3};
    Kokkos::LayoutStride values_layout(
            batch_size,
            non_zero_per_col * mat_size,
            mat_size,
            1,
            non_zero_per_col,
            mat_size);
    Kokkos::View<double***, Kokkos::LayoutStride, Kokkos::DefaultExecutionSpace>
            values_view("values", values_layout);
    Kokkos::View<int**, Kokkos::LayoutLeft, Kokkos::DefaultExecutionSpace>
            idx_view("idx", mat_size, non_zero_per_col);
    Kokkos::View<double***, Kokkos::LayoutStride, Kokkos::DefaultHostExecutionSpace>
            values_host(values, values_layout);
    Kokkos::View<int**, Kokkos::LayoutLeft, Kokkos::DefaultHostExecutionSpace>
            idx_host(col_idx, mat_size, non_zero_per_col);
    Kokkos::deep_copy(idx_view, idx_host);
    Kokkos::deep_copy(values_view, values_host);
    MatrixBatchEll<Kokkos::DefaultExecutionSpace> test_instance(idx_view, values_view, 1000, 1e-8);
    auto [idx, vals] = test_instance.get_batch_idx_and_vals();

    Kokkos::deep_copy(idx_host, idx);
    Kokkos::deep_copy(values_host, vals);

    ASSERT_EQ(idx_host(0, 0), 0);
    ASSERT_EQ(idx_host(1, 0), 1);
    ASSERT_EQ(idx_host(2, 0), 2);
    ASSERT_EQ(idx_host(3, 0), 3);
    ASSERT_EQ(values_host(0, 0, 0), 1.0);
    ASSERT_EQ(values_host(0, 1, 0), 2.0);
    ASSERT_EQ(values_host(0, 2, 0), 3.0);
    ASSERT_EQ(values_host(0, 3, 0), 4.0);
    ASSERT_EQ(values_host(1, 0, 0), 5.0);
    ASSERT_EQ(values_host(1, 1, 0), 6.0);
    ASSERT_EQ(values_host(1, 2, 0), 7.0);
    ASSERT_EQ(values_host(1, 3, 0), 8.0);
}

TEST(MatrixBatchEllFixture, SolveDiagonal)
{
    int const batch_size = 2;
    int const mat_size = 4;
    int const non_zero_per_col = 1;

    double values[] = {2.0, 4.0, 6.0, 8.0, 3.0, 5.0, 7.0, 9.0};
    int col_idx[] = {0, 1, 2, 3};
    double res[] = {0, 3.5, 2., 3., 42., 17., 0.5, 1.};
    double solution[] = {0, 3.5 / 4., 2. / 6., 3. / 8., 42. / 3., 17. / 5., 1. / 14., 1. / 9.};

    Kokkos::LayoutStride values_layout(
            batch_size,
            non_zero_per_col * mat_size,
            mat_size,
            1,
            non_zero_per_col,
            mat_size);
    Kokkos::View<double***, Kokkos::LayoutStride, Kokkos::DefaultHostExecutionSpace>
            values_host(values, values_layout);
    Kokkos::View<int**, Kokkos::LayoutLeft, Kokkos::DefaultHostExecutionSpace>
            idx_host(col_idx, mat_size, non_zero_per_col);
    Kokkos::View<double**, Kokkos::LayoutRight, Kokkos::DefaultHostExecutionSpace>
            res_host(res, batch_size, mat_size);

    Kokkos::View<double***, Kokkos::LayoutStride, Kokkos::DefaultExecutionSpace>
            values_view("values", values_layout);
    Kokkos::View<int**, Kokkos::LayoutLeft, Kokkos::DefaultExecutionSpace>
            idx_view("col_idx", mat_size, non_zero_per_col);
    Kokkos::View<double**, Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace>
            res_view("res", batch_size, mat_size);

    fill_values(idx_view);
    Kokkos::deep_copy(values_view, values_host);
    Kokkos::deep_copy(res_view, res_host);

    MatrixBatchEll<Kokkos::DefaultExecutionSpace> test_instance(idx_view, values_view, 1000, 1e-6);
    test_instance.factorize();
    test_instance.solve_inplace(res_view);

    Kokkos::deep_copy(res_host, res_view);
    ASSERT_EQ(res_host(0, 0), 0.);
    ASSERT_FLOAT_EQ(res_host(0, 1), solution[1]);
    ASSERT_FLOAT_EQ(res_host(0, 2), solution[2]);
    ASSERT_FLOAT_EQ(res_host(0, 3), solution[3]);
    ASSERT_FLOAT_EQ(res_host(1, 0), solution[4]);
    ASSERT_FLOAT_EQ(res_host(1, 1), solution[5]);
    ASSERT_FLOAT_EQ(res_host(1, 2), solution[6]);
    ASSERT_FLOAT_EQ(res_host(1, 3), solution[7]);
}

TEST(MatrixBatchEllFixture, SolveSparse)
{
    int const batch_size = 2;
    int const mat_size = 5;
    int const non_zero_per_col = 2;

    double solution[]
            = {2. / 3.,
               1. / 9.,
               7. / 9.,
               -1. / 9.,
               2. / 9.,
               2. / 3.,
               1. / 9.,
               7. / 9.,
               -1. / 9.,
               2. / 9.};
    double matvalues[2][5][5]
            = {{{1., 3., 0, 0, 0},
                {2., 0, 0, 3., 0},
                {0, 2., 1., 0, 0},
                {0, 0, 2., 5., 0},
                {0, 0, 1., 0, 1.}},
               {{1., 3., 0, 0, 0},
                {2., 0, 0, 3., 0},
                {0, 2., 1., 0, 0},
                {0, 0, 2., 5., 0},
                {0, 0, 1., 0, 1.}}};
    Kokkos::LayoutStride values_layout(
            batch_size,
            non_zero_per_col * mat_size,
            mat_size,
            1,
            non_zero_per_col,
            mat_size);

    Kokkos::View<double***, Kokkos::LayoutStride, Kokkos::DefaultHostExecutionSpace>
            values_host("test_vals", values_layout);
    Kokkos::View<int**, Kokkos::LayoutLeft, Kokkos::DefaultHostExecutionSpace>
            idx_host("col_idx", mat_size, non_zero_per_col);
    int cpt;
    for (int batch_idx = 0; batch_idx < batch_size; batch_idx++) {
        for (int i = 0; i < mat_size; i++) {
            cpt = 0;
            for (int j = 0; j < mat_size; j++) {
                if (abs(matvalues[batch_idx][i][j]) > 1e-16 && cpt < non_zero_per_col) {
                    idx_host(i, cpt) = j;
                    values_host(batch_idx, i, cpt) = matvalues[batch_idx][i][j];
                    cpt++;
                }
            }
        }
    }

    Kokkos::View<double**, Kokkos::LayoutRight, Kokkos::DefaultHostExecutionSpace>
            res_host("h_res", batch_size, mat_size);
    Kokkos::View<double***, Kokkos::LayoutStride, Kokkos::DefaultExecutionSpace>
            values_view("values", values_layout);
    Kokkos::View<int**, Kokkos::LayoutLeft, Kokkos::DefaultExecutionSpace>
            idx_view("col_idx", mat_size, non_zero_per_col);
    Kokkos::View<double**, Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace>
            res_view("res", batch_size, mat_size);

    Kokkos::deep_copy(idx_view, idx_host);
    Kokkos::deep_copy(values_view, values_host);
    Kokkos::deep_copy(res_view, 1.);
    MatrixBatchEll<Kokkos::DefaultExecutionSpace> test_instance(idx_view, values_view, 1000, 1e-12);
    test_instance.factorize();
    test_instance.solve_inplace(res_view);
    ASSERT_EQ(test_instance.norm(0), 7);
    Kokkos::deep_copy(res_host, res_view);
    ASSERT_FLOAT_EQ(res_host(0, 0), solution[0]);
    ASSERT_FLOAT_EQ(res_host(0, 1), solution[1]);
    ASSERT_FLOAT_EQ(res_host(0, 2), solution[2]);
    ASSERT_FLOAT_EQ(res_host(0, 3), solution[3]);
    ASSERT_FLOAT_EQ(res_host(0, 4), solution[4]);
}
