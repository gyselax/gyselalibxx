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
#include "mpi_spline_builder_local_crouseilles.hpp"

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

    double dx = ddc::discrete_space<Grid>().step();

    std::size_t n_local_cells(n_cells / mpi_size);
    double rank_len = dx * n_local_cells;
    Coord<Dim> local_x_min(x_min + rank_len * mpi_rank);
    Coord<Dim> local_x_max(x_min + rank_len * (mpi_rank + 1));
    ddc::init_discrete_space<BSplines>(local_x_min, local_x_max, IdxStep<Grid>(n_local_cells));

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
struct BSplinesX : ddc::UniformBSplines<X, 3>
{
};
struct TestGridX : UniformGridBase<X>
{
};

void test_MPISplineBuilderLocalCrouseilles_ReproducesLinearFunction()
{
    using Builder = MPISplineBuilderLocalCrouseilles<
            Kokkos::DefaultExecutionSpace,
            Kokkos::DefaultExecutionSpace::memory_space,
            BSplinesX,
            IdxRange<GridX>,
            ddc::SplineBuilderClosure::HERMITE,
            ddc::SplineBuilderClosure::HERMITE>;

    std::size_t constexpr n_cells = 20;
    Coord<X> const x_min(0.0);
    Coord<X> const x_max(1.0);

    IdxRange<GridX> const idx_range_x = init_local_mesh<
            X,
            GridX,
            BSplinesX>(MPI_COMM_WORLD, x_min, x_max, IdxStep<GridX>(n_cells));

    Coord<X> const local_x_min(ddc::coordinate(idx_range_x.front()));
    Coord<X> const local_x_max(ddc::coordinate(idx_range_x.back()));
    IdxRange<TestGridX> const idx_range_test_x = ddc::init_discrete_space<TestGridX>(
            TestGridX::template init<
                    TestGridX>(local_x_min, local_x_max, IdxStep<TestGridX>(n_cells * 3 + 2)));

    IdxRange<BSplinesX> const coeff_idx_range(ddc::discrete_space<BSplinesX>().full_domain());
    DFieldMem<IdxRange<BSplinesX>> coeffs_alloc(coeff_idx_range);
    DField<IdxRange<BSplinesX>> coeffs(coeffs_alloc);

    Builder builder(idx_range_x, MPI_COMM_WORLD);

    double constexpr slope = 1.3;
    double constexpr intercept = 0.7;

    DFieldMem<IdxRange<GridX>> vals_alloc(idx_range_x);
    DField<IdxRange<GridX>> vals(vals_alloc);

    IdxRange<ddc::Deriv<X>> const deriv_idx_range(builder.batched_derivs_xmin_domain(idx_range_x));
    DFieldMem<IdxRange<ddc::Deriv<X>>> derivs_alloc(deriv_idx_range);
    DField<IdxRange<ddc::Deriv<X>>> derivs(derivs_alloc);

    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            idx_range_x,
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
            BSplinesX,
            GridX,
            ddc::NullExtrapolationRule,
            ddc::NullExtrapolationRule>
            evaluator(extrapolation, extrapolation);

    ddc::host_for_each(idx_range_test_x, [&](Idx<TestGridX> const idx) {
        Coord<X> const coord(ddc::coordinate(idx));
        double const expected = slope * double(coord) + intercept;
        double const actual = evaluator(coord, get_const_field(coeffs));
        EXPECT_NEAR(actual, expected, 1e-8);
    });
}

TEST(MPISplineBuilderLocalCrouseilles, ReproducesLinearFunction)
{
    test_MPISplineBuilderLocalCrouseilles_ReproducesLinearFunction();
}

namespace {

template <class Grid, class BSplines>
double run_cosine_case(std::size_t n_cells, IdxRange<TestGridX> idx_range_test)
{
    using Builder = MPISplineBuilderLocalCrouseilles<
            Kokkos::DefaultExecutionSpace,
            Kokkos::DefaultExecutionSpace::memory_space,
            BSplines,
            IdxRange<Grid>,
            ddc::SplineBuilderClosure::PERIODIC,
            ddc::SplineBuilderClosure::PERIODIC>;

    IdxRange<Grid> const idx_range_x = init_local_mesh<X, Grid, BSplines>(
            MPI_COMM_WORLD,
            ddc::coordinate(idx_range_test.front()),
            ddc::coordinate(idx_range_test.back()),
            n_cells);

    Coord<X> const local_x_min(ddc::coordinate(idx_range_x.front()));
    Coord<X> const local_x_max(ddc::coordinate(idx_range_x.back()));

    DFieldMem<IdxRange<Grid>> vals_alloc(idx_range_x);
    DField<IdxRange<Grid>> vals(vals_alloc);
    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            idx_range_x,
            KOKKOS_LAMBDA(Idx<Grid> const idx) {
                double const x = ddc::coordinate(idx);
                vals(idx) = Kokkos::cos(x);
            });

    IdxRange<BSplines> const coeff_idx_range(ddc::discrete_space<BSplines>().full_domain());
    DFieldMem<IdxRange<BSplines>> coeffs_alloc(coeff_idx_range);
    DField<IdxRange<BSplines>> coeffs(coeffs_alloc);

    Builder builder(idx_range_x, MPI_COMM_WORLD);
    builder(coeffs, get_const_field(vals));

    ddc::NullExtrapolationRule extrapolation;
    ddc::SplineEvaluator<
            Kokkos::DefaultExecutionSpace,
            Kokkos::DefaultExecutionSpace::memory_space,
            BSplines,
            Grid,
            ddc::NullExtrapolationRule,
            ddc::NullExtrapolationRule>
            evaluator(extrapolation, extrapolation);

    double max_error = 0.0;
    ddc::host_for_each(idx_range_test, [&](Idx<TestGridX> const idx) {
        Coord<X> const coord(ddc::coordinate(idx));
        double const expected = std::cos(double(coord));
        if (coord >= local_x_min && coord <= local_x_max) {
            double const actual = evaluator(coord, get_const_field(coeffs));
            max_error = std::max(max_error, std::abs(actual - expected));
        }
    });
    return max_error;
}

} // namespace

struct GridXCoarse : UniformGridBase<X>
{
};
struct BSplinesXCoarse : ddc::UniformBSplines<X, 3>
{
};

struct GridXFine : UniformGridBase<X>
{
};
struct BSplinesXFine : ddc::UniformBSplines<X, 3>
{
};

TEST(MPISplineBuilderLocalCrouseilles, ConvergesForCosine)
{
    Coord<X> const x_min(0.0);
    Coord<X> const x_max(2.0 * M_PI);
    IdxRange<TestGridX> const idx_range_test_x = ddc::init_discrete_space<TestGridX>(
            TestGridX::init<TestGridX>(x_min, x_max, IdxStep<TestGridX>(242)));

    double const error_coarse = run_cosine_case<GridXCoarse, BSplinesXCoarse>(40, idx_range_test_x);
    double const error_fine = run_cosine_case<GridXFine, BSplinesXFine>(80, idx_range_test_x);
    double order = std::log(error_coarse / error_fine) / std::log(2.0);

    // The placeholder weights (see MPISplineBuilderLocalCrouseilles' class documentation) only
    // give a low-order-consistent boundary derivative estimate, not yet the full cubic-spline
    // accuracy of the real Crouseilles weights, so we only check the error goes down as the
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
