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
    static bool constexpr PERIODIC = false;
};

using CoordX = Coord<X>;

struct BSplinesX : ddc::UniformBSplines<X, 3>
{
};

using SplineInterpPointsX = ddc::
        GrevilleInterpolationPoints<BSplinesX, ddc::BoundCond::GREVILLE, ddc::BoundCond::GREVILLE>;

struct GridX : SplineInterpPointsX::interpolation_discrete_dimension_type
{
};

using IdxX = Idx<GridX>;
using IdxStepX = IdxStep<GridX>;
using IdxRangeX = IdxRange<GridX>;

using SplineXBuilder_1d = ddc::SplineBuilder<
        Kokkos::DefaultExecutionSpace,
        Kokkos::DefaultExecutionSpace::memory_space,
        BSplinesX,
        GridX,
        ddc::BoundCond::GREVILLE,
        ddc::BoundCond::GREVILLE,
        ddc::SplineSolver::LAPACK,
        GridX>;

using SplineXEvaluator_1d = ddc::SplineEvaluator<
        Kokkos::DefaultExecutionSpace,
        Kokkos::DefaultExecutionSpace::memory_space,
        BSplinesX,
        GridX,
        ddc::NullExtrapolationRule,
        ddc::NullExtrapolationRule,
        GridX>;

using DFieldMemX = DFieldMem<IdxRangeX>;

TEST(FemNonPeriodicPoissonSolver, Ordering)
{
    CoordX const x_min(0.0);
    CoordX const x_max(M_PI);
    IdxStepX const x_size(100);

    // Creating mesh & supports
    ddc::init_discrete_space<BSplinesX>(x_min, x_max, x_size);

    ddc::init_discrete_space<GridX>(SplineInterpPointsX::get_sampling<GridX>());
    IdxRange<GridX> gridx(SplineInterpPointsX::get_domain<GridX>());

    SplineXBuilder_1d const builder_x(gridx);

    ddc::NullExtrapolationRule x_extrapolation_rule_min;
    ddc::NullExtrapolationRule x_extrapolation_rule_max;

    SplineXEvaluator_1d const
            spline_x_evaluator(x_extrapolation_rule_min, x_extrapolation_rule_max);

    FEM1DPoissonSolver poisson(builder_x, spline_x_evaluator);

    host_t<DFieldMemX> electrostatic_potential_host(gridx);
    host_t<DFieldMemX> electric_field_host(gridx);
    host_t<DFieldMemX> rhs_host(gridx);

    // Initialization of the distribution function --> fill values
    for (IdxX const ix : gridx) {
        rhs_host(ix) = sin(ddc::coordinate(ix));
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
        double const exact_pot = sin(ddc::coordinate(ix));
        error_pot = fmax(fabs(electrostatic_potential_host(ix) - exact_pot), error_pot);
        double const exact_field = -cos(ddc::coordinate(ix));
        error_field = fmax(fabs(electric_field_host(ix) - exact_field), error_field);
    }
    EXPECT_LE(error_pot, 1e-2);
    EXPECT_LE(error_field, 1e-1);
}

} // namespace
