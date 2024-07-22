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

struct RDimX
{
    /// @brief A boolean indicating if the dimension is periodic.
    static bool constexpr PERIODIC = true;
};

struct RDimBatch
{
};

using CoordX = ddc::Coordinate<RDimX>;
using CoordB = ddc::Coordinate<RDimBatch>;

struct BSplinesX : ddc::UniformBSplines<RDimX, 3>
{
};

using SplineInterpPointsX = ddc::
        GrevilleInterpolationPoints<BSplinesX, ddc::BoundCond::PERIODIC, ddc::BoundCond::PERIODIC>;

struct IDimX : SplineInterpPointsX::interpolation_discrete_dimension_type
{
};

struct IDimBatch : ddc::UniformPointSampling<RDimBatch>
{
};

using IndexX = ddc::DiscreteElement<IDimX>;
using IVectX = ddc::DiscreteVector<IDimX>;
using IDomainX = ddc::DiscreteDomain<IDimX>;

using IndexB = ddc::DiscreteElement<IDimBatch>;
using IVectB = ddc::DiscreteVector<IDimBatch>;
using IDomainB = ddc::DiscreteDomain<IDimBatch>;

using IndexBX = ddc::DiscreteElement<IDimBatch, IDimX>;
using IVectBX = ddc::DiscreteVector<IDimBatch, IDimX>;
using IDomainBX = ddc::DiscreteDomain<IDimBatch, IDimX>;

using SplineXBuilder_1d = ddc::SplineBuilder<
        Kokkos::DefaultExecutionSpace,
        Kokkos::DefaultExecutionSpace::memory_space,
        BSplinesX,
        IDimX,
        ddc::BoundCond::PERIODIC,
        ddc::BoundCond::PERIODIC,
        ddc::SplineSolver::LAPACK,
        IDimX>;

using SplineXEvaluator_1d = ddc::SplineEvaluator<
        Kokkos::DefaultExecutionSpace,
        Kokkos::DefaultExecutionSpace::memory_space,
        BSplinesX,
        IDimX,
        ddc::PeriodicExtrapolationRule<RDimX>,
        ddc::PeriodicExtrapolationRule<RDimX>,
        IDimX>;

using BatchedSplineXBuilder_1d = ddc::SplineBuilder<
        Kokkos::DefaultExecutionSpace,
        Kokkos::DefaultExecutionSpace::memory_space,
        BSplinesX,
        IDimX,
        ddc::BoundCond::PERIODIC,
        ddc::BoundCond::PERIODIC,
        ddc::SplineSolver::LAPACK,
        IDimBatch,
        IDimX>;

using BatchedSplineXEvaluator_1d = ddc::SplineEvaluator<
        Kokkos::DefaultExecutionSpace,
        Kokkos::DefaultExecutionSpace::memory_space,
        BSplinesX,
        IDimX,
        ddc::PeriodicExtrapolationRule<RDimX>,
        ddc::PeriodicExtrapolationRule<RDimX>,
        IDimBatch,
        IDimX>;

using DFieldX = device_t<ddc::Chunk<double, IDomainX>>;
using DFieldBX = device_t<ddc::Chunk<double, IDomainBX>>;

TEST(FemPeriodicPoissonSolver, CosineSource)
{
    CoordX const x_min(0.0);
    CoordX const x_max(2.0 * M_PI);
    IVectX const x_size(100);

    // Creating mesh & supports
    ddc::init_discrete_space<BSplinesX>(x_min, x_max, x_size);
    ddc::init_discrete_space<IDimX>(SplineInterpPointsX::get_sampling<IDimX>());
    IDomainX gridx(SplineInterpPointsX::get_domain<IDimX>());

    SplineXBuilder_1d const builder_x(gridx);

    ddc::PeriodicExtrapolationRule<RDimX> x_extrapolation_rule_min;
    ddc::PeriodicExtrapolationRule<RDimX> x_extrapolation_rule_max;
    SplineXEvaluator_1d const
            spline_x_evaluator(x_extrapolation_rule_min, x_extrapolation_rule_max);

    FEM1DPoissonSolver poisson(builder_x, spline_x_evaluator);

    host_t<DFieldX> electrostatic_potential_host(gridx);
    host_t<DFieldX> electric_field_host(gridx);
    host_t<DFieldX> rhs_host(gridx);

    // Initialization of the distribution function --> fill values
    for (IndexX const ix : gridx) {
        rhs_host(ix) = cos(ddc::coordinate(ix));
    }
    DFieldX electrostatic_potential(gridx);
    DFieldX electric_field(gridx);
    DFieldX rhs(gridx);

    ddc::parallel_deepcopy(rhs, rhs_host);
    poisson(electrostatic_potential, electric_field, rhs);
    ddc::parallel_deepcopy(electric_field_host, electric_field);
    ddc::parallel_deepcopy(electrostatic_potential_host, electrostatic_potential);

    double error_pot = 0.0;
    double error_field = 0.0;


    for (IndexX const ix : gridx) {
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
    IVectX const x_size(100);

    // Creating mesh & supports
    ddc::init_discrete_space<BSplinesX>(x_min, x_max, x_size);
    ddc::init_discrete_space<IDimX>(SplineInterpPointsX::get_sampling<IDimX>());
    IDomainX gridx(SplineInterpPointsX::get_domain<IDimX>());

    CoordB const b_min(1.0);
    CoordB const b_max(3.0);
    IVectB const b_size(3);

    IDomainB gridb = ddc::init_discrete_space<IDimBatch>(IDimBatch::init(b_min, b_max, b_size));

    IDomainBX gridbx(gridb, gridx);

    BatchedSplineXBuilder_1d const builder_x(gridbx);

    ddc::PeriodicExtrapolationRule<RDimX> x_extrapolation_rule_min;
    ddc::PeriodicExtrapolationRule<RDimX> x_extrapolation_rule_max;
    BatchedSplineXEvaluator_1d const
            spline_x_evaluator(x_extrapolation_rule_min, x_extrapolation_rule_max);

    FEM1DPoissonSolver poisson(builder_x, spline_x_evaluator);

    host_t<DFieldBX> electrostatic_potential_host(gridbx);
    host_t<DFieldBX> electric_field_host(gridbx);
    host_t<DFieldBX> rhs_host(gridbx);

    // Initialization of the distribution function --> fill values
    for (IndexB const ib : gridb) {
        for (IndexX const ix : gridx) {
            rhs_host(ib, ix) = cos(ddc::coordinate(ib) * ddc::coordinate(ix));
        }
    }
    DFieldBX electrostatic_potential(gridbx);
    DFieldBX electric_field(gridbx);
    DFieldBX rhs(gridbx);

    ddc::parallel_deepcopy(rhs, rhs_host);
    poisson(electrostatic_potential, electric_field, rhs);
    ddc::parallel_deepcopy(electric_field_host, electric_field);
    ddc::parallel_deepcopy(electrostatic_potential_host, electrostatic_potential);


    for (IndexB const ib : gridb) {
        double mult = ddc::coordinate(ib);
        double error_field = 0.0;
        double error_pot = 0.0;
        for (IndexX const ix : gridx) {
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
