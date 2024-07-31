#include <ddc/ddc.hpp>

#include <sll/matrix_batch_csr.hpp>

#include <gtest/gtest.h>

#include "sll/view.hpp"

#include "test_utils.hpp"

template <MatrixBatchCsrSolver Solver>
void solve_system(
        Kokkos::View<double**, Kokkos::LayoutRight, Kokkos::DefaultHostExecutionSpace>
                values_view_host,
        Kokkos::View<int*, Kokkos::LayoutRight, Kokkos::DefaultHostExecutionSpace> idx_view_host,
        Kokkos::View<int*, Kokkos::LayoutRight, Kokkos::DefaultHostExecutionSpace>
                nnz_per_row_view_host,
        Kokkos::View<double*, Kokkos::LayoutRight, Kokkos::DefaultHostExecutionSpace> solution)
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
    Kokkos::View<double**, Kokkos::LayoutRight, Kokkos::DefaultHostExecutionSpace>
            res_host("res_host", batch_size, mat_size);

    Kokkos::deep_copy(values_view, values_view_host);
    Kokkos::deep_copy(idx_view, idx_view_host);
    Kokkos::deep_copy(nnz_per_row_view, nnz_per_row_view_host);
    Kokkos::deep_copy(res_view, 1.);

    MatrixBatchCsr<Kokkos::DefaultExecutionSpace, Solver>
            test_instance(values_view, idx_view, nnz_per_row_view);

    test_instance.setup_solver();
    test_instance.solve(res_view);
    Kokkos::deep_copy(res_host, res_view);

    for (int batch_idx = 0; batch_idx < batch_size; batch_idx++) {
        for (int i = 0; i < mat_size; i++) {
            ASSERT_FLOAT_EQ(res_host(batch_idx, i), solution(batch_idx * mat_size + i));
        }
    }
}


class MatrixBatchCsrFixture : public ::testing::Test
{
protected:
    MatrixBatchCsr<Kokkos::DefaultExecutionSpace> matrix_batch_test;
};

TEST(MatrixBatchCsrFixture, Coo_to_Csr)
{
    int const batch_size = 2;
    int const mat_size = 5;
    int const non_zero_per_system = 16;
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

    Kokkos::View<double*, Kokkos::LayoutRight, Kokkos::DefaultHostExecutionSpace>
            vals_coo_host("vals_coo", batch_size * non_zero_per_system);
    Kokkos::View<int*, Kokkos::LayoutRight, Kokkos::DefaultHostExecutionSpace>
            row_coo_host("row_coo", batch_size * non_zero_per_system);
    Kokkos::View<int*, Kokkos::LayoutRight, Kokkos::DefaultHostExecutionSpace>
            col_coo_host("col_coo", batch_size * non_zero_per_system);

    int elem_idx = 0;
    for (int batch_idx = 0; batch_idx < batch_size; batch_idx++) {
        for (int row = 0; row < mat_size; row++) {
            for (int col = 0; col < mat_size; col++) {
                if (matvalues[batch_idx][row][col] != 0.) {
                    vals_coo_host(elem_idx) = matvalues[batch_idx][row][col];
                    row_coo_host(elem_idx) = row;
                    col_coo_host(elem_idx) = col;
                    elem_idx++;
                }
            }
        }
    }
    ASSERT_EQ(elem_idx, batch_size * non_zero_per_system);

    std::unique_ptr<MatrixBatchCsr<Kokkos::DefaultExecutionSpace>> test_instance
            = std::make_unique<MatrixBatchCsr<
                    Kokkos::DefaultExecutionSpace>>(batch_size, mat_size, non_zero_per_system);
    convert_coo_to_csr(test_instance, vals_coo_host, row_coo_host, col_coo_host);
    test_instance->setup_solver();
    Kokkos::View<double**, Kokkos::LayoutRight, Kokkos::DefaultHostExecutionSpace>
            res_host("h_res", batch_size, mat_size);
    Kokkos::View<double**, Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace>
            res_view("res", batch_size, mat_size);
    Kokkos::deep_copy(res_view, 1.);

    test_instance->solve(res_view);
    Kokkos::deep_copy(res_host, res_view);

    for (int batch_idx = 0; batch_idx < batch_size; batch_idx++) {
        for (int row = 0; row < mat_size; row++) {
            ASSERT_FLOAT_EQ(res_host(batch_idx, row), solution[batch_idx * mat_size + row]);
        }
    }
}

template <MatrixBatchCsrSolver Solver>
void solve_diagonal_system()
{
    int const batch_size = 2;
    int const mat_size = 4;
    int const non_zero_per_system = 4;
    double values[] = {2.0, 4.0, 6.0, 8.0, 3.0, 5.0, 7.0, 9.0};
    int col_idxs[] = {0, 1, 2, 3};
    int nnz_per_row[] = {0, 1, 2, 3, 4};
    // rhs and solution
    double res[] = {0, 3.5, 2., 3., 42., 17., 0.5, 1.};
    double solution[] = {1. / 2., 1. / 4., 1. / 6., 1. / 8., 1. / 3., 1. / 5., 1. / 7., 1. / 9.};

    Kokkos::View<double**, Kokkos::LayoutRight, Kokkos::DefaultHostExecutionSpace>
            values_view_host(values, batch_size, non_zero_per_system);
    Kokkos::View<int*, Kokkos::LayoutRight, Kokkos::DefaultHostExecutionSpace>
            idx_view_host(col_idxs, non_zero_per_system);
    Kokkos::View<int*, Kokkos::LayoutRight, Kokkos::DefaultHostExecutionSpace>
            nnz_per_row_view_host(nnz_per_row, mat_size + 1);
    Kokkos::View<double**, Kokkos::LayoutRight, Kokkos::DefaultHostExecutionSpace>
            res_view_host(res, batch_size, mat_size);
    Kokkos::View<double*, Kokkos::LayoutRight, Kokkos::DefaultHostExecutionSpace>
            solution_view_host(solution, batch_size * mat_size);

    solve_system<
            Solver>(values_view_host, idx_view_host, nnz_per_row_view_host, solution_view_host);
}

template <MatrixBatchCsrSolver Solver>
void solve_sparse_system()
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
    Kokkos::View<double*, Kokkos::LayoutRight, Kokkos::DefaultHostExecutionSpace>
            solution_view_host(solution, batch_size * mat_size);
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
    solve_system<
            Solver>(values_view_host, idx_view_host, nnz_per_row_view_host, solution_view_host);
}

template <MatrixBatchCsrSolver Solver>
void solve_pds_system()
{
    {
        int const batch_size = 2;
        int const mat_size = 5;
        int const non_zero_per_system = 25;
        int nnz_per_row[] = {0, 5, 10, 15, 20, 25};
        //first matrix eigenvalues [7.13490895 1.12269422 0.2177588  0.07179924 0.02072621]
        //second matrix eigenvalues  [9.38119529e+00 9.57862018e-01 5.16641290e-01 5.50717953e-02 4.69290376e-03]
        double solution[]
                = {-2.37358759,
                   2.41177772,
                   1.09774264,
                   1.10503935,
                   -0.16610463,
                   -19.035696,
                   27.611921,
                   -34.03087866,
                   30.813921,
                   11.64331918};
        double matvalues[2][5][5]
                = {{{2.85679107, 1.86116608, 1.55237696, 1.67953204, 1.6130445},
                    {1.86116608, 1.39832276, 1.07542783, 0.98889335, 1.37330181},
                    {1.55237696, 1.07542783, 1.06195425, 1.00023903, 1.08393101},
                    {1.67953204, 0.98889335, 1.00023903, 1.44583069, 0.56698767},
                    {1.6130445, 1.3733018, 1.08393101, 0.56698767, 1.80498864}},
                   {{1.67297988, 1.33253271, 2.15265423, 1.76538283, 1.28063759},
                    {1.33253271, 1.26654954, 1.75269646, 1.18861185, 1.23794696},
                    {2.15265423, 1.75269646, 3.36746571, 2.65006331, 2.27777785},
                    {1.76538283, 1.18861185, 2.65006331, 2.4502457, 1.41434842},
                    {1.28063759, 1.23794696, 2.27777785, 1.41434842, 2.15822247}}};

        Kokkos::View<double**, Kokkos::LayoutRight, Kokkos::DefaultHostExecutionSpace>
                values_view_host("values_host", batch_size, non_zero_per_system);
        Kokkos::View<int*, Kokkos::LayoutRight, Kokkos::DefaultHostExecutionSpace>
                idx_view_host("col_idxs_host", non_zero_per_system);
        Kokkos::View<int*, Kokkos::LayoutRight, Kokkos::DefaultHostExecutionSpace>
                nnz_per_row_view_host(nnz_per_row, mat_size + 1);
        Kokkos::View<double*, Kokkos::LayoutRight, Kokkos::DefaultHostExecutionSpace>
                solution_view_host(solution, batch_size * mat_size);
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
        solve_system<
                Solver>(values_view_host, idx_view_host, nnz_per_row_view_host, solution_view_host);
    }
}

TEST(MatrixBatchCsrFixture, SolveDiagonalCg)
{
    solve_diagonal_system<MatrixBatchCsrSolver::CG>();
}

TEST(MatrixBatchCsrFixture, SolveDiagonalBatchCg)
{
    solve_diagonal_system<MatrixBatchCsrSolver::BATCH_CG>();
}

/* Disabled because of Ginkgo issue #1563 (OpenMP-specific)
TEST(MatrixBatchCsrFixture, SolveDiagonalBicgstab)
{
    solve_diagonal_system<MatrixBatchCsrSolver::BICGSTAB>();
}
*/

TEST(MatrixBatchCsrFixture, SolveDiagonalBatchBicgstab)
{
    solve_diagonal_system<MatrixBatchCsrSolver::BATCH_BICGSTAB>();
}

TEST(MatrixBatchCsrFixture, SolvePDSCg)
{
    solve_pds_system<MatrixBatchCsrSolver::CG>();
}

TEST(MatrixBatchCsrFixture, SolvePDSBatchCg)
{
    solve_pds_system<MatrixBatchCsrSolver::BATCH_CG>();
}

TEST(MatrixBatchCsrFixture, SolveSparseBicgstab)
{
    solve_sparse_system<MatrixBatchCsrSolver::BICGSTAB>();
}

TEST(MatrixBatchCsrFixture, SolveSparseBatchBicgstab)
{
    solve_sparse_system<MatrixBatchCsrSolver::BATCH_BICGSTAB>();
}
