// SPDX-License-Identifier: MIT

#include <ddc/ddc.hpp>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "ddc_helper.hpp"
#include "fem_1d_poisson_solver.hpp"
#include "neumann_spline_quadrature.hpp"
#include "quadrature.hpp"
#include "species_info.hpp"

namespace {

struct X
{
    /// @brief A boolean indicating if the dimension is periodic.
    static bool constexpr PERIODIC = true;
};

struct Batch
{
};

using CoordX = Coord<X>;
using CoordBatch = Coord<Batch>;

struct BSplinesX : ddc::UniformBSplines<X, 3>
{
};

using SplineInterpPointsX = ddc::
        GrevilleInterpolationPoints<BSplinesX, ddc::BoundCond::PERIODIC, ddc::BoundCond::PERIODIC>;

struct GridX : SplineInterpPointsX::interpolation_discrete_dimension_type
{
};

struct GridBatch : UniformGridBase<Batch>
{
};

using IdxX = Idx<GridX>;
using IdxStepX = IdxStep<GridX>;
using IdxRangeX = IdxRange<GridX>;

using IdxBatch = Idx<GridBatch>;
using IdxStepBatch = IdxStep<GridBatch>;
using IdxRangeBatch = IdxRange<GridBatch>;

using IdxRangeBatchX = IdxRange<GridBatch, GridX>;

using SplineXBuilder_1d = ddc::SplineBuilder<
        Kokkos::DefaultExecutionSpace,
        Kokkos::DefaultExecutionSpace::memory_space,
        BSplinesX,
        GridX,
        ddc::BoundCond::PERIODIC,
        ddc::BoundCond::PERIODIC,
        ddc::SplineSolver::LAPACK,
        GridX>;

using SplineXEvaluator_1d = ddc::SplineEvaluator<
        Kokkos::DefaultExecutionSpace,
        Kokkos::DefaultExecutionSpace::memory_space,
        BSplinesX,
        GridX,
        ddc::PeriodicExtrapolationRule<X>,
        ddc::PeriodicExtrapolationRule<X>,
        GridX>;

using BatchedSplineXBuilder_1d = ddc::SplineBuilder<
        Kokkos::DefaultExecutionSpace,
        Kokkos::DefaultExecutionSpace::memory_space,
        BSplinesX,
        GridX,
        ddc::BoundCond::PERIODIC,
        ddc::BoundCond::PERIODIC,
        ddc::SplineSolver::LAPACK,
        GridBatch,
        GridX>;

using BatchedSplineXEvaluator_1d = ddc::SplineEvaluator<
        Kokkos::DefaultExecutionSpace,
        Kokkos::DefaultExecutionSpace::memory_space,
        BSplinesX,
        GridX,
        ddc::PeriodicExtrapolationRule<X>,
        ddc::PeriodicExtrapolationRule<X>,
        GridBatch,
        GridX>;

using DFieldMemX = DFieldMem<IdxRangeX>;
using DFieldMemBatchX = DFieldMem<IdxRangeBatchX>;

TEST(FemPeriodicPoissonSolver, CosineSource)
{
    CoordX const x_min(0.0);
    CoordX const x_max(2.0 * M_PI);
    IdxStepX const x_size(100);

    // Creating mesh & supports
    ddc::init_discrete_space<BSplinesX>(x_min, x_max, x_size);
    ddc::init_discrete_space<GridX>(SplineInterpPointsX::get_sampling<GridX>());
    IdxRangeX gridx(SplineInterpPointsX::get_domain<GridX>());

    SplineXBuilder_1d const builder_x(gridx);

    ddc::PeriodicExtrapolationRule<X> x_extrapolation_rule_min;
    ddc::PeriodicExtrapolationRule<X> x_extrapolation_rule_max;
    SplineXEvaluator_1d const
            spline_x_evaluator(x_extrapolation_rule_min, x_extrapolation_rule_max);

    FEM1DPoissonSolver poisson(builder_x, spline_x_evaluator);

    host_t<DFieldMemX> electrostatic_potential_host(gridx);
    host_t<DFieldMemX> electric_field_host(gridx);
    host_t<DFieldMemX> rhs_host(gridx);

    // Initialization of the distribution function --> fill values
    for (IdxX const ix : gridx) {
        rhs_host(ix) = cos(ddc::coordinate(ix));
    }
    DFieldMemX electrostatic_potential(gridx);
    DFieldMemX electric_field(gridx);
    DFieldMemX rhs(gridx);

    ddc::parallel_deepcopy(rhs, rhs_host);
    poisson(electrostatic_potential, electric_field, rhs);
    ddc::parallel_deepcopy(electric_field_host, electric_field);
    ddc::parallel_deepcopy(electrostatic_potential_host, electrostatic_potential);

    double error_pot = 0.0;
    double error_field = 0.0;


    for (IdxX const ix : gridx) {
        double const exact_pot = cos(ddc::coordinate(ix));
        error_pot = fmax(fabs(electrostatic_potential_host(ix) - exact_pot), error_pot);
        double const exact_field = sin(ddc::coordinate(ix));
        error_field = fmax(fabs(electric_field_host(ix) - exact_field), error_field);
    }
    EXPECT_LE(error_pot, 1e-8);
    EXPECT_LE(error_field, 1e-6);
}

TEST(FemPeriodicPoissonSolver, BatchedCosineSource)
{
    CoordX const x_min(0.0);
    CoordX const x_max(2.0 * M_PI);
    IdxStepX const x_size(100);

    // Creating mesh & supports
    ddc::init_discrete_space<BSplinesX>(x_min, x_max, x_size);
    ddc::init_discrete_space<GridX>(SplineInterpPointsX::get_sampling<GridX>());
    IdxRangeX gridx(SplineInterpPointsX::get_domain<GridX>());

    CoordBatch const b_min(1.0);
    CoordBatch const b_max(3.0);
    IdxStepBatch const b_size(3);

    IdxRangeBatch gridb
            = ddc::init_discrete_space<GridBatch>(GridBatch::init(b_min, b_max, b_size));

    IdxRangeBatchX gridbx(gridb, gridx);

    BatchedSplineXBuilder_1d const builder_x(gridbx);

    ddc::PeriodicExtrapolationRule<X> x_extrapolation_rule_min;
    ddc::PeriodicExtrapolationRule<X> x_extrapolation_rule_max;
    BatchedSplineXEvaluator_1d const
            spline_x_evaluator(x_extrapolation_rule_min, x_extrapolation_rule_max);

    FEM1DPoissonSolver poisson(builder_x, spline_x_evaluator);

    host_t<DFieldMemBatchX> electrostatic_potential_host(gridbx);
    host_t<DFieldMemBatchX> electric_field_host(gridbx);
    host_t<DFieldMemBatchX> rhs_host(gridbx);

    // Initialization of the distribution function --> fill values
    for (IdxBatch const ib : gridb) {
        for (IdxX const ix : gridx) {
            rhs_host(ib, ix) = cos(ddc::coordinate(ib) * ddc::coordinate(ix));
        }
    }
    DFieldMemBatchX electrostatic_potential(gridbx);
    DFieldMemBatchX electric_field(gridbx);
    DFieldMemBatchX rhs(gridbx);

    ddc::parallel_deepcopy(rhs, rhs_host);
    poisson(electrostatic_potential, electric_field, rhs);
    ddc::parallel_deepcopy(electric_field_host, electric_field);
    ddc::parallel_deepcopy(electrostatic_potential_host, electrostatic_potential);


    for (IdxBatch const ib : gridb) {
        double mult = ddc::coordinate(ib);
        double error_field = 0.0;
        double error_pot = 0.0;
        for (IdxX const ix : gridx) {
            double const exact_pot = cos(mult * ddc::coordinate(ix)) / (mult * mult);
            error_pot = fmax(fabs(electrostatic_potential_host(ib, ix) - exact_pot), error_pot);
            double const exact_field = sin(mult * ddc::coordinate(ix)) / mult;
            error_field = fmax(fabs(electric_field_host(ib, ix) - exact_field), error_field);
        }
        EXPECT_LE(error_pot, mult * 1e-8);
        EXPECT_LE(error_field, mult * 1e-6);
    }
}

} // namespace
