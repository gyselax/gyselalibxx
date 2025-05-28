// SPDX-License-Identifier: MIT
#include <gtest/gtest.h>

#include "../test_utils.hpp"

#include "cartesian_to_circular.hpp"
#include "cartesian_to_czarny.hpp"
#include "circular_to_cartesian.hpp"
#include "crank_nicolson.hpp"
#include "czarny_to_cartesian.hpp"
#include "ddc_aliases.hpp"
#include "euler.hpp"
#include "geometry_pseudo_cartesian.hpp"
#include "mesh_builder.hpp"
#include "r_theta_test_cases.hpp"
#include "rk2.hpp"
#include "rk3.hpp"
#include "rk4.hpp"
#include "species_info.hpp"
#include "spline_polar_foot_finder.hpp"
#include "vector_field.hpp"
#include "vector_field_mem.hpp"
#include "vector_index_tools.hpp"

namespace {
struct X
{
    static constexpr bool IS_CONTRAVARIANT = true;
    static constexpr bool IS_COVARIANT = true;
    using Dual = X;
};
struct Y
{
    static constexpr bool IS_CONTRAVARIANT = true;
    static constexpr bool IS_COVARIANT = true;
    using Dual = Y;
};

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

struct GridR : NonUniformGridBase<R>
{
};
struct GridTheta : NonUniformGridBase<Theta>
{
};

static constexpr int BSDegree = 3;
static constexpr ddc::BoundCond SplineRBoundary = ddc::BoundCond::GREVILLE;
static constexpr ddc::BoundCond SplineThetaBoundary = ddc::BoundCond::PERIODIC;


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
        SplineRBoundary, // boundary at r=0
        SplineRBoundary, // boundary at rmax
        SplineThetaBoundary,
        SplineThetaBoundary,
        ddc::SplineSolver::LAPACK>;

using SplineRThetaEvaluator = ddc::SplineEvaluator2D<
        Kokkos::DefaultExecutionSpace,
        typename Kokkos::DefaultExecutionSpace::memory_space,
        BSplinesR,
        BSplinesTheta,
        GridR,
        GridTheta,
        ddc::NullExtrapolationRule, // boundary at r=0
        ddc::ConstantExtrapolationRule<R, Theta>, // boundary at rmax
        ddc::PeriodicExtrapolationRule<Theta>,
        ddc::PeriodicExtrapolationRule<Theta>>;

using SplineInterpPointsR
        = ddc::GrevilleInterpolationPoints<BSplinesR, SplineRBoundary, SplineRBoundary>;
using SplineInterpPointsTheta
        = ddc::GrevilleInterpolationPoints<BSplinesTheta, SplineThetaBoundary, SplineThetaBoundary>;

using IdxRangeR = IdxRange<GridR>;
using IdxRangeTheta = IdxRange<GridTheta>;
using IdxRangeRTheta = IdxRange<GridR, GridTheta>;
using IdxRangeSpRTheta = IdxRange<Species, GridR, GridTheta>;
using IdxRTheta = Idx<GridR, GridTheta>;
using IdxSpRTheta = Idx<Species, GridR, GridTheta>;

using CartBasis = VectorIndexSet<X, Y>;

using CoordXY = Coord<X, Y>;
using CoordRTheta = Coord<R, Theta>;

struct AnalyticalCircular
{
    using LogicalToPhysicalMapping = CircularToCartesian<R, Theta, X, Y>;
    using LogicalToPseudoPhysicalMapping = CircularToCartesian<R, Theta, X, Y>;
};

struct AnalyticalCzarny
{
    using LogicalToPhysicalMapping = CzarnyToCartesian<R, Theta, X, Y>;
    using LogicalToPseudoPhysicalMapping = CzarnyToCartesian<R, Theta, X, Y>;
};

struct PseudoCartCzarny
{
    using LogicalToPhysicalMapping = CzarnyToCartesian<R, Theta, X, Y>;
    using LogicalToPseudoPhysicalMapping = CircularToCartesian<R, Theta, X_pC, Y_pC>;
};

template <class T>
struct PolarAdvectionFixture;

template <class TimeStepperBuilderType, class Mappings, class AdvectionFieldType>
struct PolarAdvectionFixture<std::tuple<TimeStepperBuilderType, Mappings, AdvectionFieldType>>
    : public testing::Test
{
    using LogicalToPhysicalMapping = typename Mappings::LogicalToPhysicalMapping;
    using LogicalToPseudoPhysicalMapping = typename Mappings::LogicalToPseudoPhysicalMapping;
    using X_adv = typename LogicalToPseudoPhysicalMapping::cartesian_tag_x;
    using Y_adv = typename LogicalToPseudoPhysicalMapping::cartesian_tag_y;
    using TimeStepperBuilder = TimeStepperBuilderType;
    using AdvectionField = AdvectionFieldType;
};

template <class LogicalToOtherMapping>
LogicalToOtherMapping init_mapping()
{
    using OtherX = typename LogicalToOtherMapping::cartesian_tag_x;
    using OtherY = typename LogicalToOtherMapping::cartesian_tag_y;
    // At x0,y0 to match rotation centre
    double x0 = 0.0;
    double y0 = 0.0;
    if constexpr (std::is_same_v<
                          LogicalToOtherMapping,
                          CircularToCartesian<R, Theta, OtherX, OtherY>>) {
        return LogicalToOtherMapping(x0, y0);
    } else if constexpr (std::is_same_v<
                                 LogicalToOtherMapping,
                                 CzarnyToCartesian<R, Theta, OtherX, OtherY>>) {
        return LogicalToOtherMapping(0.3, 1.4, y0, y0);
    }
}

template <class AdvectionField>
AdvectionField init_field()
{
    static_assert(
            std::is_same_v<
                    AdvectionField,
                    AdvectionField_translation<
                            X,
                            Y>> || std::is_same_v<AdvectionField, AdvectionField_rotation<X, Y, R, Theta>> || std::is_same_v<AdvectionField, AdvectionField_decentred_rotation<X, Y>>);
    if constexpr (std::is_same_v<AdvectionField, AdvectionField_translation<X, Y>>) {
        return AdvectionField(DVector<X, Y>(
                std::cos(2 * M_PI * 511. / 4096.) / 2.,
                std::sin(2 * M_PI * 511. / 4096.) / 2.));
    } else if constexpr (std::is_same_v<AdvectionField, AdvectionField_rotation<X, Y, R, Theta>>) {
        return AdvectionField(DVector<R, Theta>(0., 2 * M_PI));
    } else if constexpr (std::is_same_v<AdvectionField, AdvectionField_decentred_rotation<X, Y>>) {
        return AdvectionField();
    }
}

using TimeSteppers = std::tuple<RK4Builder>;
using Mappings = std::tuple<AnalyticalCircular, AnalyticalCzarny, PseudoCartCzarny>;
using AdvectionFieldTypes = std::tuple<
        AdvectionField_translation<X, Y>,
        AdvectionField_rotation<X, Y, R, Theta>,
        AdvectionField_decentred_rotation<X, Y>>;

using Cases = tuple_to_types_t<cartesian_product_t<TimeSteppers, Mappings, AdvectionFieldTypes>>;


TYPED_TEST_SUITE(PolarAdvectionFixture, Cases);

TYPED_TEST(PolarAdvectionFixture, Analytical)
{
    using LogicalToPhysicalMapping = typename TestFixture::LogicalToPhysicalMapping;
    using LogicalToPseudoPhysicalMapping = typename TestFixture::LogicalToPseudoPhysicalMapping;
    using TimeStepperBuilder = typename TestFixture::TimeStepperBuilder;
    using AdvectionField = typename TestFixture::AdvectionField;

    Coord<R> const r_min(0.0);
    Coord<R> const r_max(1.0);
    IdxStep<GridR> const nr_cells(40);

    Coord<Theta> const theta_min(0.0);
    Coord<Theta> const theta_max(2.0 * M_PI);
    IdxStep<GridTheta> const ntheta_cells(30);

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

    LogicalToPhysicalMapping to_physical = init_mapping<LogicalToPhysicalMapping>();
    LogicalToPseudoPhysicalMapping to_pseudo_physical
            = init_mapping<LogicalToPseudoPhysicalMapping>();
    inverse_mapping_t<LogicalToPhysicalMapping> from_physical = to_physical.get_inverse_mapping();

    TimeStepperBuilder time_stepper;
    AdvectionField advection_field = init_field<AdvectionField>();

    SplinePolarFootFinder const batched_foot_finder(
            batched_idx_range,
            time_stepper,
            to_physical,
            to_pseudo_physical,
            builder,
            evaluator);

    const double t = 0.0;
    const double dt = 0.001;
    DVectorFieldMem<IdxRangeSpRTheta, CartBasis> adv_field_alloc(batched_idx_range);
    FieldMem<CoordRTheta, IdxRangeSpRTheta> feet_alloc(batched_idx_range);
    FieldMem<CoordRTheta, IdxRangeSpRTheta> exact_feet_alloc(batched_idx_range);
    DVectorField<IdxRangeSpRTheta, CartBasis> adv_field = get_field(adv_field_alloc);
    Field<CoordRTheta, IdxRangeSpRTheta> feet = get_field(feet_alloc);
    Field<CoordRTheta, IdxRangeSpRTheta> exact_feet = get_field(exact_feet_alloc);
    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            batched_idx_range,
            KOKKOS_LAMBDA(IdxSpRTheta idx) {
                IdxRTheta idx_rtheta(idx);
                CoordRTheta coord_rtheta = ddc::coordinate(idx_rtheta);
                CoordXY coord_xy = to_physical(coord_rtheta);
                ddcHelper::
                        assign_vector_field_element(adv_field, idx, advection_field(coord_xy, t));
                feet(idx) = coord_rtheta;
                exact_feet(idx) = from_physical(advection_field.exact_feet(coord_xy, dt));
            });
    batched_foot_finder(feet, adv_field, dt);

    double error = error_norm_inf(
            Kokkos::DefaultExecutionSpace(),
            get_const_field(feet),
            get_const_field(exact_feet));

    double TOL = 1e-4;
    EXPECT_NEAR(error, 0.0, TOL);
}

}; // namespace
