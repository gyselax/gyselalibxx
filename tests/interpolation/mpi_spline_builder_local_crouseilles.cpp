// SPDX-License-Identifier: MIT

#include <mpi.h>

#include <algorithm>
#include <cmath>

#include <ddc/ddc.hpp>
#include <ddc/kernels/splines.hpp>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "ddc_alias_inline_functions.hpp"
#include "ddc_aliases.hpp"
#include "mpi_scope_guard.hpp"
#include "mpi_spline_builder_local_vidal.hpp"

namespace {

/// @brief Build a rank-local, non-periodic uniform cubic B-spline basis + interpolation grid
/// covering [x_min, x_min + n_cells * dx], and return the resulting local interpolation index
/// range. Every rank calls this with a distinct Dim/Grid/BSplines type so that each type's
/// (process-local, per-MPI-rank) ddc discrete space is only ever initialised once.
template <class Dim, class Grid, class BSplines>
IdxRange<Grid> init_local_mesh(
        MPI_Comm comm,
        Coord<Dim> x_min,
        Coord<Dim> x_max,
        std::size_t n_cells)
{
    int mpi_size, mpi_rank;
    MPI_Comm_size(comm, &mpi_size);
    MPI_Comm_rank(comm, &mpi_rank);

    assert(n_cells % mpi_size == 0);

    ddc::init_discrete_space<Grid>(
            Grid::template init<Grid>(x_min, x_max, IdxStep<Grid>(n_cells + 1)));

    ddc::init_discrete_space<BSplines>(x_min, x_max, IdxStep<Grid>(n_cells));

    std::size_t n_local_cells(n_cells / mpi_size);
    return IdxRange<Grid>(
            ddc::discrete_space<Grid>().front() + n_local_cells * mpi_rank,
            IdxStep<Grid>(n_local_cells + 1));
}

} // namespace

struct X
{
    static bool constexpr PERIODIC = false;
};
struct GridX : UniformGridBase<X>
{
};
struct LocalBSplinesXLower : ddc::UniformBSplines<X, 3, false>
{
};
struct LocalBSplinesX : ddc::UniformBSplines<X, 3, false>
{
};
struct LocalBSplinesXUpper : ddc::UniformBSplines<X, 3, false>
{
};
struct BSplinesX : ddc::UniformBSplines<X, 3, false>
{
};
struct TestGridX : UniformGridBase<X>
{
};

void test_MPISplineBuilderLocalVidal_ReproducesLinearFunction()
{
    using Builder = MPISplineBuilderLocalVidal<
            Kokkos::DefaultExecutionSpace,
            Kokkos::DefaultExecutionSpace::memory_space,
            LocalBSplinesXLower,
            LocalBSplinesX,
            LocalBSplinesXUpper,
            BSplinesX,
            GridX,
            ddc::SplineBuilderClosure::HERMITE,
            ddc::SplineBuilderClosure::HERMITE>;

    std::size_t constexpr n_cells = 20;
    Coord<X> const x_min(0.0);
    Coord<X> const x_max(1.0);

    IdxRange<GridX> const local_idx_range_x = init_local_mesh<
            X,
            GridX,
            BSplinesX>(MPI_COMM_WORLD, x_min, x_max, IdxStep<GridX>(n_cells));
    IdxRange<GridX> const global_idx_range_x(Idx<GridX> {0}, IdxStep<GridX> {n_cells + 1});

    construct_local_bsplines_and_neighbours<
            BSplinesX,
            LocalBSplinesXLower,
            LocalBSplinesX,
            LocalBSplinesXUpper>(MPI_COMM_WORLD);

    Coord<X> const local_x_min(ddc::coordinate(local_idx_range_x.front()));
    Coord<X> const local_x_max(ddc::coordinate(local_idx_range_x.back()));
    IdxRange<TestGridX> const idx_range_test_x = ddc::init_discrete_space<TestGridX>(
            TestGridX::template init<
                    TestGridX>(local_x_min, local_x_max, IdxStep<TestGridX>(n_cells)));

    IdxRange<LocalBSplinesX> const coeff_idx_range(
            ddc::discrete_space<LocalBSplinesX>().full_domain());
    DFieldMem<IdxRange<LocalBSplinesX>> coeffs_alloc(coeff_idx_range);
    DField<IdxRange<LocalBSplinesX>> coeffs(coeffs_alloc);

    Builder builder(local_idx_range_x, global_idx_range_x, MPI_COMM_WORLD, 10);

    double constexpr slope = 1.3;
    double constexpr intercept = 0.7;

    DFieldMem<IdxRange<GridX>> vals_alloc(local_idx_range_x);
    DField<IdxRange<GridX>> vals(vals_alloc);

    IdxRange<ddc::Deriv<X>> const deriv_idx_range(
            builder.batched_derivs_xmin_domain(local_idx_range_x));
    DFieldMem<IdxRange<ddc::Deriv<X>>> derivs_alloc(deriv_idx_range);
    DField<IdxRange<ddc::Deriv<X>>> derivs(derivs_alloc);

    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            local_idx_range_x,
            KOKKOS_LAMBDA(Idx<GridX> const idx) {
                double const x = ddc::coordinate(idx);
                vals(idx) = slope * x + intercept;
            });
    Idx<ddc::Deriv<X>> first_deriv(1);
    Kokkos::parallel_for(
            "Fill deriv",
            Kokkos::RangePolicy<Kokkos::DefaultExecutionSpace>(0, 1),
            KOKKOS_LAMBDA(const int) { derivs(first_deriv) = slope; });

    builder(coeffs,
            get_const_field(vals),
            std::optional(get_const_field(derivs)),
            std::optional(get_const_field(derivs)));

    ddc::NullExtrapolationRule extrapolation;
    ddc::SplineEvaluator<
            Kokkos::DefaultExecutionSpace,
            Kokkos::DefaultExecutionSpace::memory_space,
            LocalBSplinesX,
            GridX,
            ddc::NullExtrapolationRule,
            ddc::NullExtrapolationRule>
            evaluator(extrapolation, extrapolation);

    double max_error = ddc::parallel_transform_reduce(
            Kokkos::DefaultExecutionSpace(),
            idx_range_test_x,
            0.,
            ddc::reducer::max<double>(),
            KOKKOS_LAMBDA(Idx<TestGridX> const idx) {
                Coord<X> const coord(ddc::coordinate(idx));
                double const expected = slope * double(coord) + intercept;
                double const actual = evaluator(coord, get_const_field(coeffs));
                return Kokkos::abs(expected - actual);
            });

    EXPECT_LT(max_error, 5e-7);
}

TEST(MPISplineBuilderLocalVidal, ReproducesLinearFunction)
{
    test_MPISplineBuilderLocalVidal_ReproducesLinearFunction();
}

namespace {

template <
        class Grid,
        class BSplines,
        class LocalBSplinesLower,
        class LocalBSplines,
        class LocalBSplinesUpper>
double run_cosine_case(std::size_t n_cells, IdxRange<TestGridX> idx_range_test)
{
    using Builder = MPISplineBuilderLocalVidal<
            Kokkos::DefaultExecutionSpace,
            Kokkos::DefaultExecutionSpace::memory_space,
            LocalBSplinesLower,
            LocalBSplines,
            LocalBSplinesUpper,
            BSplines,
            Grid,
            ddc::SplineBuilderClosure::PERIODIC,
            ddc::SplineBuilderClosure::PERIODIC>;

    IdxRange<Grid> const local_idx_range_x = init_local_mesh<X, Grid, BSplines>(
            MPI_COMM_WORLD,
            ddc::coordinate(idx_range_test.front()),
            ddc::coordinate(idx_range_test.back()),
            n_cells);
    IdxRange<Grid> const global_idx_range_x(Idx<Grid> {0}, IdxStep<Grid> {n_cells + 1});

    construct_local_bsplines_and_neighbours<
            BSplines,
            LocalBSplinesLower,
            LocalBSplines,
            LocalBSplinesUpper>(MPI_COMM_WORLD);

    Coord<X> const local_x_min(ddc::coordinate(local_idx_range_x.front()));
    Coord<X> const local_x_max(ddc::coordinate(local_idx_range_x.back()));

    DFieldMem<IdxRange<Grid>> vals_alloc(local_idx_range_x);
    DField<IdxRange<Grid>> vals(vals_alloc);
    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            local_idx_range_x,
            KOKKOS_LAMBDA(Idx<Grid> const idx) {
                double const x = ddc::coordinate(idx);
                vals(idx) = Kokkos::cos(x);
            });

    IdxRange<LocalBSplines> const coeff_idx_range(
            ddc::discrete_space<LocalBSplines>().full_domain());
    DFieldMem<IdxRange<LocalBSplines>> coeffs_alloc(coeff_idx_range);
    DField<IdxRange<LocalBSplines>> coeffs(coeffs_alloc);

    Builder builder(local_idx_range_x, global_idx_range_x, MPI_COMM_WORLD, 10);

    builder(coeffs, get_const_field(vals));

    ddc::NullExtrapolationRule extrapolation;
    ddc::SplineEvaluator<
            Kokkos::DefaultExecutionSpace,
            Kokkos::DefaultExecutionSpace::memory_space,
            LocalBSplines,
            Grid,
            ddc::NullExtrapolationRule,
            ddc::NullExtrapolationRule>
            evaluator(extrapolation, extrapolation);

    double max_error = ddc::parallel_transform_reduce(
            Kokkos::DefaultExecutionSpace(),
            idx_range_test,
            0.,
            ddc::reducer::max<double>(),
            KOKKOS_LAMBDA(Idx<TestGridX> const idx) {
                Coord<X> const coord(ddc::coordinate(idx));
                double const expected = Kokkos::cos(double(coord));
                if (coord >= local_x_min && coord <= local_x_max) {
                    double const actual = evaluator(coord, get_const_field(coeffs));
                    return Kokkos::abs(actual - expected);
                } else {
                    return 0.0;
                }
            });
    return max_error;
}

} // namespace

struct GridXCoarse : UniformGridBase<X>
{
};
struct BSplinesXCoarse : ddc::UniformBSplines<X, 3, true>
{
};
struct LocalBSplinesXCoarseLower : ddc::UniformBSplines<X, 3, false>
{
};
struct LocalBSplinesXCoarse : ddc::UniformBSplines<X, 3, false>
{
};
struct LocalBSplinesXCoarseUpper : ddc::UniformBSplines<X, 3, false>
{
};

struct GridXFine : UniformGridBase<X>
{
};
struct BSplinesXFine : ddc::UniformBSplines<X, 3, true>
{
};
struct LocalBSplinesXFineLower : ddc::UniformBSplines<X, 3, false>
{
};
struct LocalBSplinesXFine : ddc::UniformBSplines<X, 3, false>
{
};
struct LocalBSplinesXFineUpper : ddc::UniformBSplines<X, 3, false>
{
};

TEST(MPISplineBuilderLocalVidal, ConvergesForCosine)
{
    Coord<X> const x_min(0.0);
    Coord<X> const x_max(2.0 * M_PI);
    IdxRange<TestGridX> const idx_range_test_x = ddc::init_discrete_space<TestGridX>(
            TestGridX::init<TestGridX>(x_min, x_max, IdxStep<TestGridX>(242)));

    double const error_coarse = run_cosine_case<
            GridXCoarse,
            BSplinesXCoarse,
            LocalBSplinesXCoarseLower,
            LocalBSplinesXCoarse,
            LocalBSplinesXCoarseUpper>(40, idx_range_test_x);
    double const error_fine = run_cosine_case<
            GridXFine,
            BSplinesXFine,
            LocalBSplinesXFineLower,
            LocalBSplinesXFine,
            LocalBSplinesXFineUpper>(80, idx_range_test_x);
    double order = std::log(error_coarse / error_fine) / std::log(2.0);

    // The placeholder weights (see MPISplineBuilderLocalVidal' class documentation) only
    // give a low-order-consistent boundary derivative estimate, not yet the full cubic-spline
    // accuracy of the real Vidal weights, so we only check the error goes down as the
    // resolution increases, not a specific convergence order.
    EXPECT_LT(error_fine, error_coarse);
    EXPECT_NEAR(order, 4.0, 0.5);
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleMock(&argc, argv);
    MpiScopeGuard mpi_scope(argc, argv);
    ::Kokkos::ScopeGuard kokkos_scope(argc, argv);
    ::ddc::ScopeGuard ddc_scope(argc, argv);
    ::testing::TestEventListeners& listeners = ::testing::UnitTest::GetInstance()->listeners();
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank != 0) {
        delete listeners.Release(listeners.default_result_printer());
    }
    return RUN_ALL_TESTS();
}
