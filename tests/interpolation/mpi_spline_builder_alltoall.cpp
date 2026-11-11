// SPDX-License-Identifier: MIT
#include <cmath>

#include <ddc/ddc.hpp>
#include <ddc/kernels/splines.hpp>

#include <gtest/gtest.h>

#include "ddc_alias_inline_functions.hpp"
#include "ddc_aliases.hpp"
#include "i_interpolation_builder.hpp"
#include "mpi_builder_alltoall.hpp"

namespace {

struct X
{
    static bool constexpr PERIODIC = true;
};

int constexpr Degree = 3;

struct BSplinesX : ddc::UniformBSplines<X, Degree, X::PERIODIC>
{
};

using SplineInterpPointsX = ddc::GrevilleInterpolationPoints<
        BSplinesX,
        ddc::SplineBuilderClosure::PERIODIC,
        ddc::SplineBuilderClosure::PERIODIC>;

struct GridX : SplineInterpPointsX::interpolation_discrete_dimension_type
{
};

// Batch dimensions (no coordinates are needed)
struct GridY : ddc::DiscreteDimension
{
};
struct GridZ : ddc::DiscreteDimension
{
};

using SplineXBuilder = ddc::SplineBuilder<
        Kokkos::DefaultExecutionSpace,
        Kokkos::DefaultExecutionSpace::memory_space,
        BSplinesX,
        GridX,
        ddc::SplineBuilderClosure::PERIODIC,
        ddc::SplineBuilderClosure::PERIODIC,
        ddc::SplineSolver::LAPACK>;

using IdxRangeX = IdxRange<GridX>;
using IdxRangeY = IdxRange<GridY>;
using IdxRangeZ = IdxRange<GridZ>;
using IdxRangeXY = IdxRange<GridX, GridY>;
using IdxRangeYXZ = IdxRange<GridY, GridX, GridZ>;

int constexpr NCells = 16;

/// A phase which differs for each batch index.
template <class... BatchDims>
double batch_phase(IdxStep<BatchDims...> ib)
{
    double phase = 0.0;
    ((phase = 1.7 * phase + IdxStep<BatchDims>(ib).value()), ...);
    return 0.1 * phase;
}

using MPIBuilder = MPIBuilderAlltoall<SplineXBuilder>;

template <class IdxRangeBatchedInterp>
void test_mpi_spline_builder(
        SplineXBuilder const& builder,
        MPIBuilder const& mpi_builder,
        IdxRangeBatchedInterp global_idx_range)
{
    using IdxBatchedInterp = typename IdxRangeBatchedInterp::discrete_element_type;
    using IdxRangeBatchedSpline =
            typename MPIBuilder::template batched_spline_domain_type<IdxRangeBatchedInterp>;
    using IdxBatchedSpline = typename IdxRangeBatchedSpline::discrete_element_type;
    using IdxRangeBatch = ddc::remove_dims_of_t<IdxRangeBatchedInterp, GridX>;
    using IdxBatch = typename IdxRangeBatch::discrete_element_type;

    int rank;
    int comm_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

    IdxRangeX const idx_range_x(global_idx_range);
    IdxRangeBatch const idx_range_batch(global_idx_range);

    // Serial reference computed on every rank
    DFieldMem<IdxRangeBatchedInterp> global_vals_alloc(global_idx_range);
    DField<IdxRangeBatchedInterp> global_vals(global_vals_alloc);
    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            global_idx_range,
            KOKKOS_LAMBDA(IdxBatchedInterp idx) {
                double const x = ddc::coordinate(Idx<GridX>(idx));
                double const phase = batch_phase(IdxBatch(idx) - idx_range_batch.front());
                global_vals(idx) = Kokkos::cos(x + phase) + 0.5 * Kokkos::sin(3 * x - phase);
            });
    DFieldMem<IdxRangeBatchedSpline> global_coeffs_alloc(
            builder.batched_spline_domain(global_idx_range));
    builder(get_field(global_coeffs_alloc), get_const_field(global_vals));
    auto global_coeffs = ddc::create_mirror_view_and_copy(get_field(global_coeffs_alloc));

    // Distributed computation
    IdxRangeX const local_idx_range_x = mpi_builder.local_interpolation_domain();
    IdxRangeBatchedInterp const local_idx_range(local_idx_range_x, idx_range_batch);
    IdxRangeBatchedSpline const local_spline_idx_range
            = mpi_builder.batched_spline_domain(local_idx_range);

    // Periodic splines: there is one interpolation point per cell
    int const nloc = local_idx_range_x.size();
    EXPECT_EQ(nloc * comm_size, NCells);
    EXPECT_EQ(local_idx_range_x.front(), idx_range_x.front() + rank * nloc);
    EXPECT_EQ(ddc::select<BSplinesX>(local_spline_idx_range).size(), nloc + Degree);

    DFieldMem<IdxRangeBatchedInterp> local_vals_alloc(local_idx_range);
    ddc::parallel_deepcopy(get_field(local_vals_alloc), global_vals[local_idx_range]);

    DFieldMem<IdxRangeBatchedSpline> local_coeffs_alloc(local_spline_idx_range);
    mpi_builder(get_field(local_coeffs_alloc), get_const_field(local_vals_alloc));

    auto local_coeffs_host = ddc::create_mirror_view_and_copy(get_field(local_coeffs_alloc));

    // Compare against the serial reference
    Idx<BSplinesX> const spline_front
            = ddc::select<BSplinesX>(get_idx_range(global_coeffs)).front();
    int const nbasis = ddc::discrete_space<BSplinesX>().nbasis();
    EXPECT_EQ(ddc::select<BSplinesX>(local_spline_idx_range).front(), spline_front + rank * nloc);
    ddc::host_for_each(local_spline_idx_range, [&](IdxBatchedSpline idx) {
        Idx<BSplinesX> const ib_local(idx);
        Idx<BSplinesX> const ib_global(spline_front + (ib_local - spline_front).value() % nbasis);
        double const expected = global_coeffs(IdxBatchedSpline(ib_global, IdxBatch(idx)));
        EXPECT_NEAR(local_coeffs_host(idx), expected, 1e-14 * std::max(1.0, std::abs(expected)));
    });
}

} // namespace

TEST(MPIBuilderAlltoall, Concept)
{
    static_assert(concepts::InterpolationBuilder1D<MPIBuilder>);
}

TEST(MPIBuilderAlltoall, Build2D)
{
    ddc::init_discrete_space<BSplinesX>(Coord<X>(0.0), Coord<X>(2 * M_PI), NCells);
    ddc::init_discrete_space<GridX>(SplineInterpPointsX::get_sampling<GridX>());
    IdxRangeX const idx_range_x(SplineInterpPointsX::get_domain<GridX>());
    IdxRangeY const idx_range_y(Idx<GridY>(0), IdxStep<GridY>(4));
    SplineXBuilder const builder(idx_range_x);
    MPIBuilder const mpi_builder(builder, MPI_COMM_WORLD);
    test_mpi_spline_builder(builder, mpi_builder, IdxRangeXY(idx_range_x, idx_range_y));
}

TEST(MPIBuilderAlltoall, Build3D)
{
    ddc::init_discrete_space<BSplinesX>(Coord<X>(0.0), Coord<X>(2 * M_PI), NCells);
    ddc::init_discrete_space<GridX>(SplineInterpPointsX::get_sampling<GridX>());
    IdxRangeX const idx_range_x(SplineInterpPointsX::get_domain<GridX>());
    // With 2 ranks Y (size 3) cannot be distributed so Z is distributed instead
    IdxRangeY const idx_range_y(Idx<GridY>(0), IdxStep<GridY>(3));
    IdxRangeZ const idx_range_z(Idx<GridZ>(0), IdxStep<GridZ>(4));
    SplineXBuilder const builder(idx_range_x);
    MPIBuilder const mpi_builder(builder, MPI_COMM_WORLD);
    test_mpi_spline_builder(
            builder,
            mpi_builder,
            IdxRangeYXZ(idx_range_y, idx_range_x, idx_range_z));
}

TEST(MPIBuilderAlltoall, ReuseForDifferentBatches)
{
    ddc::init_discrete_space<BSplinesX>(Coord<X>(0.0), Coord<X>(2 * M_PI), NCells);
    ddc::init_discrete_space<GridX>(SplineInterpPointsX::get_sampling<GridX>());
    IdxRangeX const idx_range_x(SplineInterpPointsX::get_domain<GridX>());
    IdxRangeY const idx_range_y(Idx<GridY>(0), IdxStep<GridY>(6));
    IdxRangeZ const idx_range_z(Idx<GridZ>(0), IdxStep<GridZ>(2));
    SplineXBuilder const builder(idx_range_x);
    MPIBuilder const mpi_builder(builder, MPI_COMM_WORLD);
    test_mpi_spline_builder(builder, mpi_builder, IdxRangeXY(idx_range_x, idx_range_y));
    test_mpi_spline_builder(
            builder,
            mpi_builder,
            IdxRangeYXZ(idx_range_y, idx_range_x, idx_range_z));
}
