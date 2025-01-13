
#include <gtest/gtest.h>

#include "matrix_batch_tridiag.hpp"

using ConstField2d = Kokkos::View<double**, Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace>;

void solve_batched_tridiag_system(
        int const batch_size,
        int const mat_size,
        double const sub_diag,
        double const diag,
        double const up_diag,
        Kokkos::View<double**, Kokkos::LayoutRight, Kokkos::DefaultHostExecutionSpace>
                Rhs_view_host)
{
    ConstField2d A_view("A", batch_size, mat_size);
    ConstField2d B_view("B", batch_size, mat_size);
    ConstField2d C_view("C", batch_size, mat_size);
    ConstField2d Rhs_view("R", batch_size, mat_size);
    Kokkos::deep_copy(A_view, sub_diag);
    Kokkos::deep_copy(B_view, diag);
    Kokkos::deep_copy(C_view, up_diag);
    Kokkos::deep_copy(Rhs_view, Rhs_view_host);
    MatrixBatchTridiag<Kokkos::DefaultExecutionSpace>
            matrix(batch_size, mat_size, A_view, B_view, C_view);
    matrix.setup_solver();
    matrix.solve(Rhs_view);
    Kokkos::deep_copy(Rhs_view_host, Rhs_view);
}


TEST(MatrixBatchTridiag, Symmetric)
{
    int const batch_size = 2;
    int const mat_size = 4;

    Kokkos::View<double**, Kokkos::LayoutRight, Kokkos::DefaultHostExecutionSpace>
            Res_view_host("R_host", batch_size, mat_size);
    Kokkos::deep_copy(Res_view_host, 1.);
    //Cond number for the matrix 2.72
    solve_batched_tridiag_system(batch_size, mat_size, 0.5, 2, 0.5, Res_view_host);

    ASSERT_DOUBLE_EQ(Res_view_host(0, 0), 8. / 19.);
    ASSERT_DOUBLE_EQ(Res_view_host(0, 1), 6. / 19.);
    ASSERT_DOUBLE_EQ(Res_view_host(0, 2), 6. / 19.);
    ASSERT_DOUBLE_EQ(Res_view_host(0, 3), 8. / 19.);
}

TEST(MatrixBatchTridiag, DominantDiag)
{
    int const batch_size = 1;
    int const mat_size = 4;

    Kokkos::View<double**, Kokkos::LayoutRight, Kokkos::DefaultHostExecutionSpace>
            Res_view_host("R_host", batch_size, mat_size);
    Kokkos::deep_copy(Res_view_host, 1.);
    //Cond number for the matrix 2.56
    solve_batched_tridiag_system(batch_size, mat_size, 0.5, 1.1, -0.25, Res_view_host);

    ASSERT_DOUBLE_EQ(Res_view_host(0, 0), 1.0468844955326548);
    ASSERT_DOUBLE_EQ(Res_view_host(0, 1), 0.6062917803436817);
    ASSERT_DOUBLE_EQ(Res_view_host(0, 2), 0.7614528245775094);
    ASSERT_DOUBLE_EQ(Res_view_host(0, 3), 0.5629759888284048);
}

TEST(MatrixBatchTridiag, NotValidMatrix)
{
    int const batch_size = 2;
    int const mat_size = 4;

    ConstField2d A_view("A", batch_size, mat_size);
    ConstField2d B_view("B", batch_size, mat_size);
    ConstField2d C_view("C", batch_size, mat_size);
    //neither positive definite symmetric nor diagonal dominant
    Kokkos::deep_copy(A_view, 0.93);
    Kokkos::deep_copy(B_view, -0.367);
    Kokkos::deep_copy(C_view, 1.42);
    MatrixBatchTridiag<Kokkos::DefaultExecutionSpace>
            matrix(batch_size, mat_size, A_view, B_view, C_view);
    EXPECT_FALSE(matrix.check_stability());
}

TEST(MatrixBatchTridiag, GeneralValidMatrix)
{
    int const batch_size = 1;
    int const mat_size = 4;

    Kokkos::View<double**, Kokkos::LayoutRight, Kokkos::DefaultHostExecutionSpace>
            Res_view_host("R_host", batch_size, mat_size);
    double subdiag[] = {0., 2.0, 0.9, -0.8};
    double diag[] = {2.7, 5.0, 2.0, 3.3};
    double uppdiag[] = {1.8, 2.4, -0.7, 0.};
    double r[] = {7.8, 11.5, 42., 0.5};


    Kokkos::View<double**, Kokkos::LayoutRight, Kokkos::DefaultHostExecutionSpace>
            A_view_host(subdiag, batch_size, mat_size);
    Kokkos::View<double**, Kokkos::LayoutRight, Kokkos::DefaultHostExecutionSpace>
            B_view_host(diag, batch_size, mat_size);
    Kokkos::View<double**, Kokkos::LayoutRight, Kokkos::DefaultHostExecutionSpace>
            C_view_host(uppdiag, batch_size, mat_size);
    Kokkos::View<double**, Kokkos::LayoutRight, Kokkos::DefaultHostExecutionSpace>
            Rhs_view_host(r, batch_size, mat_size);

    ConstField2d A_view("A", batch_size, mat_size);
    ConstField2d B_view("B", batch_size, mat_size);
    ConstField2d C_view("C", batch_size, mat_size);
    ConstField2d Rhs_view("R", batch_size, mat_size);
    Kokkos::deep_copy(A_view, A_view_host);
    Kokkos::deep_copy(B_view, B_view_host);
    Kokkos::deep_copy(C_view, C_view_host);
    Kokkos::deep_copy(Rhs_view, Rhs_view_host);
    MatrixBatchTridiag<Kokkos::DefaultExecutionSpace>
            matrix(batch_size, mat_size, A_view, B_view, C_view);
    matrix.solve(Rhs_view);
    Kokkos::deep_copy(Res_view_host, Rhs_view);

    ASSERT_DOUBLE_EQ(Res_view_host(0, 0), 272999. / 16896.);
    ASSERT_DOUBLE_EQ(Res_view_host(0, 1), -672565. / 33792.);
    ASSERT_DOUBLE_EQ(Res_view_host(0, 2), 134315. / 4096.);
    ASSERT_DOUBLE_EQ(Res_view_host(0, 3), 45625. / 5632.);
}
