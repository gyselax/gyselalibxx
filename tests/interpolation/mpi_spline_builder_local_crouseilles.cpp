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

// These tests are run with 2 MPI ranks (see tests/interpolation/CMakeLists.txt), each rank
// owning one half of a globally periodic domain. Every dimension tag below is used to build a
// *local*, non-periodic B-spline basis for a single rank (X::PERIODIC = false): the wrapped
// internal ddc::SplineBuilder inside MPISplineBuilderLocalCrouseilles always uses HERMITE
// closure. "PERIODIC" as a MPISplineBuilderLocalCrouseilles template argument only selects the
// wrap-around MPI rank topology (rank 0's lower neighbour is the last rank), independent of the
// (always non-periodic) local B-spline basis.

// Run entirely on the host execution space: this keeps the test focused on the algorithm
// (index/MPI logic) rather than device-memory plumbing, and matches the host-only
// SplineRThetaBuilder_host precedent in spline_interpolation_centre.cpp.
using ExecSpace = Kokkos::DefaultHostExecutionSpace;
using MemorySpace = Kokkos::HostSpace;

namespace {

/// @brief Build a rank-local, non-periodic uniform cubic B-spline basis + interpolation grid
/// covering [x_min, x_min + n_cells * dx], and return the resulting local interpolation index
/// range. Every rank calls this with a distinct Dim/Grid/BSplines type so that each type's
/// (process-local, per-MPI-rank) ddc discrete space is only ever initialised once.
template <class Dim, class Grid, class BSplines>
IdxRange<Grid> init_local_mesh(Coord<Dim> x_min, double dx, std::size_t n_cells)
{
    Coord<Dim> const x_max(double(x_min) + dx * n_cells);
    ddc::init_discrete_space<BSplines>(x_min, x_max, IdxStep<Grid>(n_cells));

    using SplineInterpPoints = ddc::GrevilleInterpolationPoints<
            BSplines,
            ddc::SplineBuilderClosure::HERMITE,
            ddc::SplineBuilderClosure::HERMITE>;
    ddc::init_discrete_space<Grid>(SplineInterpPoints::template get_sampling<Grid>());
    return SplineInterpPoints::template get_domain<Grid>();
}

int get_rank()
{
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    return rank;
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

TEST(MPISplineBuilderLocalCrouseilles, ReproducesLinearFunction)
{
    int const rank = get_rank();

    using Builder = MPISplineBuilderLocalCrouseilles<
            ExecSpace,
            MemorySpace,
            BSplinesX,
            IdxRange<GridX>,
            ddc::SplineBuilderClosure::HERMITE,
            ddc::SplineBuilderClosure::HERMITE>;

    double constexpr dx = 0.05;
    std::size_t constexpr n_cells = 20;
    // Rank 0 covers [0, n_cells*dx]; rank 1 covers [n_cells*dx, 2*n_cells*dx], so the two local
    // domains meet exactly at the shared boundary point.
    Coord<X> const local_x_min(rank * n_cells * dx);

    IdxRange<GridX> const idx_range_x
            = init_local_mesh<X, GridX, BSplinesX>(local_x_min, dx, n_cells);

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
            ExecSpace(),
            idx_range_x,
            KOKKOS_LAMBDA(Idx<GridX> const idx) {
                double const x = ddc::coordinate(idx);
                vals(idx) = slope * x + intercept;
            });
    Idx<ddc::Deriv<X>> first_deriv(1);
    derivs(first_deriv) = slope;

    std::cout << "Ready to build" << std::endl;

    builder(coeffs,
            get_const_field(vals),
            std::optional(get_const_field(derivs)),
            std::optional(get_const_field(derivs)));

    std::cout << "Built" << std::endl;

    ddc::NullExtrapolationRule extrapolation;
    ddc::SplineEvaluator<
            ExecSpace,
            MemorySpace,
            BSplinesX,
            GridX,
            ddc::NullExtrapolationRule,
            ddc::NullExtrapolationRule>
            evaluator(extrapolation, extrapolation);

    ddc::host_for_each(idx_range_x, [&](Idx<GridX> const idx) {
        Coord<X> const coord(ddc::coordinate(idx));
        double const expected = slope * double(coord) + intercept;
        double const actual = evaluator(coord, get_const_field(coeffs));
        EXPECT_NEAR(actual, expected, 1e-8);
    });
}

namespace {

template <class Dim, class Grid, class BSplines>
double run_cosine_case(std::size_t n_cells)
{
    int const rank = get_rank();

    using Builder = MPISplineBuilderLocalCrouseilles<
            ExecSpace,
            MemorySpace,
            BSplines,
            IdxRange<Grid>,
            ddc::SplineBuilderClosure::PERIODIC,
            ddc::SplineBuilderClosure::PERIODIC>;

    // Each rank covers half of one full period of cos, so the wrap-around MPI boundary sees
    // physically continuous (periodic) data.
    double constexpr period = 2.0 * M_PI;
    double const dx = period / (2 * n_cells);
    Coord<Dim> const local_x_min(rank * n_cells * dx);

    IdxRange<Grid> const idx_range_x
            = init_local_mesh<Dim, Grid, BSplines>(local_x_min, dx, n_cells);

    DFieldMem<IdxRange<Grid>> vals_alloc(idx_range_x);
    DField<IdxRange<Grid>> vals(vals_alloc);
    ddc::parallel_for_each(
            ExecSpace(),
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
            ExecSpace,
            MemorySpace,
            BSplines,
            Grid,
            ddc::NullExtrapolationRule,
            ddc::NullExtrapolationRule>
            evaluator(extrapolation, extrapolation);

    double max_error = 0.0;
    ddc::host_for_each(idx_range_x, [&](Idx<Grid> const idx) {
        Coord<Dim> const coord(ddc::coordinate(idx));
        double const expected = std::cos(double(coord));
        double const actual = evaluator(coord, get_const_field(coeffs));
        max_error = std::max(max_error, std::abs(actual - expected));
    });
    return max_error;
}

} // namespace

struct XCoarse
{
    static bool constexpr PERIODIC = false;
};
struct GridXCoarse : UniformGridBase<XCoarse>
{
};
struct BSplinesXCoarse : ddc::UniformBSplines<XCoarse, 3>
{
};

struct XFine
{
    static bool constexpr PERIODIC = false;
};
struct GridXFine : UniformGridBase<XFine>
{
};
struct BSplinesXFine : ddc::UniformBSplines<XFine, 3>
{
};

TEST(MPISplineBuilderLocalCrouseilles, ConvergesForCosine)
{
    double const error_coarse = run_cosine_case<XCoarse, GridXCoarse, BSplinesXCoarse>(20);
    double const error_fine = run_cosine_case<XFine, GridXFine, BSplinesXFine>(80);

    // The placeholder weights (see MPISplineBuilderLocalCrouseilles' class documentation) only
    // give a low-order-consistent boundary derivative estimate, not yet the full cubic-spline
    // accuracy of the real Crouseilles weights, so we only check the error goes down as the
    // resolution increases, not a specific convergence order.
    EXPECT_LT(error_fine, error_coarse);
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleMock(&argc, argv);
    MpiScopeGuard mpi_scope(argc, argv);
    ::Kokkos::ScopeGuard kokkos_scope(argc, argv);
    ::ddc::ScopeGuard ddc_scope(argc, argv);
    ::testing::TestEventListeners& listeners = ::testing::UnitTest::GetInstance()->listeners();
    if (get_rank() != 0) {
        delete listeners.Release(listeners.default_result_printer());
    }
    return RUN_ALL_TESTS();
}
