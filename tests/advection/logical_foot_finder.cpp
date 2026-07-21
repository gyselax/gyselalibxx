// SPDX-License-Identifier: MIT
#include <gtest/gtest.h>

#include "../test_utils.hpp"

#include "circular_to_cartesian.hpp"
#include "ddc_aliases.hpp"
#include "mesh_builder.hpp"
#include "polar_foot_finder.hpp"
#include "rk4.hpp"
#include "species_info.hpp"
#include "vector_field.hpp"
#include "vector_field_mem.hpp"
#include "vector_index_tools.hpp"

namespace {

struct R_cov;
struct Theta_cov;
struct R
{
    static constexpr bool PERIODIC = false;
    static constexpr bool IS_CONTRAVARIANT = true;
    static constexpr bool IS_COVARIANT = false;
    using Dual = R_cov;
};
struct Theta
{
    static constexpr bool PERIODIC = true;
    static constexpr bool IS_CONTRAVARIANT = true;
    static constexpr bool IS_COVARIANT = false;
    using Dual = Theta_cov;
};
struct R_cov
{
    static constexpr bool IS_CONTRAVARIANT = false;
    static constexpr bool IS_COVARIANT = true;
    using Dual = R;
};
struct Theta_cov
{
    static constexpr bool IS_CONTRAVARIANT = false;
    static constexpr bool IS_COVARIANT = true;
    using Dual = Theta;
};
struct X
{
    static constexpr bool PERIODIC = false;
    static constexpr bool IS_CONTRAVARIANT = true;
    static constexpr bool IS_COVARIANT = true;
    using Dual = X;
};
struct Y
{
    static constexpr bool PERIODIC = false;
    static constexpr bool IS_CONTRAVARIANT = true;
    static constexpr bool IS_COVARIANT = true;
    using Dual = Y;
};

struct GridR : NonUniformGridBase<R>
{
};
struct GridTheta : NonUniformGridBase<Theta>
{
};

static constexpr int BSDegree = 3;
static constexpr ddc::SplineBuilderClosure SplineRClosure = ddc::SplineBuilderClosure::GREVILLE;
static constexpr ddc::SplineBuilderClosure SplineThetaClosure = ddc::SplineBuilderClosure::PERIODIC;

struct BSplinesR : ddc::NonUniformBSplines<R, BSDegree>
{
};
struct BSplinesTheta : ddc::NonUniformBSplines<Theta, BSDegree>
{
};

using SplineRThetaBuilder = ddc::SplineBuilder2D<
        Kokkos::DefaultExecutionSpace,
        typename Kokkos::DefaultExecutionSpace::memory_space,
        BSplinesR,
        BSplinesTheta,
        GridR,
        GridTheta,
        SplineRClosure,
        SplineRClosure,
        SplineThetaClosure,
        SplineThetaClosure,
        ddc::SplineSolver::LAPACK>;

using SplineRThetaEvaluator = ddc::SplineEvaluator2D<
        Kokkos::DefaultExecutionSpace,
        typename Kokkos::DefaultExecutionSpace::memory_space,
        BSplinesR,
        BSplinesTheta,
        GridR,
        GridTheta,
        ddc::NullExtrapolationRule,
        ddc::ConstantExtrapolationRule<R, Theta>,
        ddc::PeriodicExtrapolationRule<Theta>,
        ddc::PeriodicExtrapolationRule<Theta>>;

using SplineInterpPointsR
        = ddc::GrevilleInterpolationPoints<BSplinesR, SplineRClosure, SplineRClosure>;
using SplineInterpPointsTheta
        = ddc::GrevilleInterpolationPoints<BSplinesTheta, SplineThetaClosure, SplineThetaClosure>;

using IdxRangeR = IdxRange<GridR>;
using IdxRangeTheta = IdxRange<GridTheta>;
using IdxRangeRTheta = IdxRange<GridR, GridTheta>;
using IdxRangeSpRTheta = IdxRange<Species, GridR, GridTheta>;
using IdxRTheta = Idx<GridR, GridTheta>;
using IdxSpRTheta = Idx<Species, GridR, GridTheta>;

using LogicalBasis = VectorIndexSet<R, Theta>;
using CoordRTheta = Coord<R, Theta>;

static constexpr Coord<R> r_min(0.1);
static constexpr Coord<R> r_max(1.0);
static constexpr Coord<Theta> theta_min(0.0);
static constexpr Coord<Theta> theta_max(2.0 * M_PI);

} // namespace

/**
 * @brief Test SplineLogicalFootFinder with a pure poloidal rotation.
 *
 * The advection field is @f$ (A^r, A^\theta) = (0, \omega) @f$ (constant everywhere).
 * For a time step @f$ dt @f$, the exact characteristic foot from @f$ (r, \theta) @f$ is:
 * @f$ (r, \theta - \omega \, dt) @f$.
 */
TEST(LogicalFootFinder, PureRotation)
{
    IdxStep<GridR> const nr_cells(25);
    IdxStep<GridTheta> const ntheta_cells(50);

    IdxStepSp const nb_kinspecies(2);
    IdxRangeSp const idx_range_sp(IdxSp(0), nb_kinspecies);

    ddc::init_discrete_space<BSplinesR>(
            build_random_non_uniform_break_points(r_min, r_max, nr_cells, 0.5));
    ddc::init_discrete_space<BSplinesTheta>(
            build_random_non_uniform_break_points(theta_min, theta_max, ntheta_cells, 0.5));

    ddc::init_discrete_space<GridR>(SplineInterpPointsR::template get_sampling<GridR>());
    ddc::init_discrete_space<GridTheta>(
            SplineInterpPointsTheta::template get_sampling<GridTheta>());

    IdxRangeR r_idx_range(SplineInterpPointsR::template get_domain<GridR>());
    IdxRangeTheta theta_idx_range(SplineInterpPointsTheta::template get_domain<GridTheta>());
    IdxRangeRTheta idx_range(r_idx_range, theta_idx_range);
    IdxRangeSpRTheta batched_idx_range(idx_range_sp, idx_range);

    ddc::NullExtrapolationRule r_min_extrap;
    ddc::PeriodicExtrapolationRule<Theta> theta_extrap;
    SplineRThetaBuilder builder(idx_range);
    ddc::ConstantExtrapolationRule<R, Theta> r_max_extrap(r_max);
    SplineRThetaEvaluator evaluator(r_min_extrap, r_max_extrap, theta_extrap, theta_extrap);

    RK4Builder time_stepper;

    CircularToCartesian<R, Theta, X, Y> mapping;

    PolarFootFinder<
            FootFindingSpace::LOGICAL,
            AdvectionFieldSpace::LOGICAL,
            CircularToCartesian<R, Theta, X, Y>,
            IdxRangeSpRTheta,
            RK4Builder,
            decltype(builder),
            decltype(evaluator)>
            foot_finder(time_stepper, mapping, builder, evaluator);

    const double omega = 2 * M_PI;
    const double dt = 0.001;

    // Build constant logical advection field (0, omega).
    DVectorFieldMem<IdxRangeSpRTheta, LogicalBasis> adv_field_alloc(batched_idx_range);
    DVectorField<IdxRangeSpRTheta, LogicalBasis> adv_field = get_field(adv_field_alloc);
    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            batched_idx_range,
            KOKKOS_LAMBDA(IdxSpRTheta idx) {
                ddcHelper::get<R>(adv_field)(idx) = 0.0;
                ddcHelper::get<Theta>(adv_field)(idx) = omega;
            });

    // Initialise feet to mesh coordinates and compute exact feet.
    FieldMem<CoordRTheta, IdxRangeSpRTheta> feet_alloc(batched_idx_range);
    FieldMem<CoordRTheta, IdxRangeSpRTheta> exact_feet_alloc(batched_idx_range);
    Field<CoordRTheta, IdxRangeSpRTheta> feet = get_field(feet_alloc);
    Field<CoordRTheta, IdxRangeSpRTheta> exact_feet = get_field(exact_feet_alloc);

    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            batched_idx_range,
            KOKKOS_LAMBDA(IdxSpRTheta idx) {
                CoordRTheta const coord = ddc::coordinate(IdxRTheta(idx));
                feet(idx) = coord;
                exact_feet(idx) = coord - dt * DVector<R, Theta>(0.0, omega);
            });

    foot_finder(feet, get_const_field(adv_field), dt);

    double const error = error_norm_inf(
            Kokkos::DefaultExecutionSpace(),
            get_const_field(feet),
            get_const_field(exact_feet));

    EXPECT_NEAR(error, 0.0, 1e-4);
}

/**
 * @brief Test SplineLogicalFootFinder O-point reflection.
 *
 * With a large outward-pointing advection field @f$ A^r > 0 @f$, points near @f$ r_{\min} @f$
 * are driven to @f$ r < 0 @f$ and must be reflected through the O-point.
 * After foot-finding, all radial coordinates must be non-negative.
 */
TEST(LogicalFootFinder, OPointReflection)
{
    IdxStep<GridR> const nr_cells(25);
    IdxStep<GridTheta> const ntheta_cells(50);

    ddc::init_discrete_space<BSplinesR>(
            build_random_non_uniform_break_points(r_min, r_max, nr_cells, 0.5));
    ddc::init_discrete_space<BSplinesTheta>(
            build_random_non_uniform_break_points(theta_min, theta_max, ntheta_cells, 0.5));

    ddc::init_discrete_space<GridR>(SplineInterpPointsR::template get_sampling<GridR>());
    ddc::init_discrete_space<GridTheta>(
            SplineInterpPointsTheta::template get_sampling<GridTheta>());

    IdxRangeR r_idx_range(SplineInterpPointsR::template get_domain<GridR>());
    IdxRangeTheta theta_idx_range(SplineInterpPointsTheta::template get_domain<GridTheta>());
    IdxRangeRTheta idx_range(r_idx_range, theta_idx_range);

    ddc::NullExtrapolationRule r_min_extrap;
    ddc::PeriodicExtrapolationRule<Theta> theta_extrap;
    SplineRThetaBuilder builder(idx_range);
    ddc::ConstantExtrapolationRule<R, Theta> r_max_extrap(r_max);
    SplineRThetaEvaluator evaluator(r_min_extrap, r_max_extrap, theta_extrap, theta_extrap);

    RK4Builder time_stepper;

    CircularToCartesian<R, Theta, X, Y> mapping;

    PolarFootFinder<
            FootFindingSpace::LOGICAL,
            AdvectionFieldSpace::LOGICAL,
            CircularToCartesian<R, Theta, X, Y>,
            IdxRangeRTheta,
            RK4Builder,
            decltype(builder),
            decltype(evaluator)>
            foot_finder(time_stepper, mapping, builder, evaluator);

    // Large positive A^r: foot_r = r - dt * vr < 0 for small r.
    // With r_min = 0.1, dt = 0.1, vr = 5.0: foot_r = r - 0.5 < 0 for r < 0.5.
    const double vr = 5.0;
    const double dt = 0.1;

    DVectorFieldMem<IdxRangeRTheta, LogicalBasis> adv_field_alloc(idx_range);
    DVectorField<IdxRangeRTheta, LogicalBasis> adv_field = get_field(adv_field_alloc);
    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            idx_range,
            KOKKOS_LAMBDA(IdxRTheta idx) {
                ddcHelper::get<R>(adv_field)(idx) = vr;
                ddcHelper::get<Theta>(adv_field)(idx) = 0.0;
            });

    FieldMem<CoordRTheta, IdxRangeRTheta> feet_alloc(idx_range);
    Field<CoordRTheta, IdxRangeRTheta> feet = get_field(feet_alloc);
    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            idx_range,
            KOKKOS_LAMBDA(IdxRTheta idx) { feet(idx) = ddc::coordinate(idx); });

    foot_finder(feet, get_const_field(adv_field), dt);

    // Mirror to host to check all r values are non-negative.
    auto feet_host = ddc::create_mirror_view_and_copy(get_const_field(feet));
    ddc::host_for_each(idx_range, [&](IdxRTheta idx) { EXPECT_GE(ddc::get<R>(feet_host(idx)), 0.0); });
}
