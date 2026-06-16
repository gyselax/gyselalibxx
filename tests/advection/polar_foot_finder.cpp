// SPDX-License-Identifier: MIT
#include <cmath>

#include <gtest/gtest.h>

#include "../test_utils.hpp"

#include "cartesian_to_circular.hpp"
#include "cartesian_to_czarny.hpp"
#include "circular_to_cartesian.hpp"
#include "czarny_to_cartesian.hpp"
#include "ddc_aliases.hpp"
#include "geometry_pseudo_cartesian.hpp"
#include "mesh_builder.hpp"
#include "polar_foot_finder.hpp"
#include "r_theta_test_cases.hpp"
#include "rk4.hpp"
#include "species_info.hpp"
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
using LogicalBasis = VectorIndexSet<R, Theta>;

using CoordXY = Coord<X, Y>;
using CoordRTheta = Coord<R, Theta>;

struct AnalyticalCircular
{
    using LogicalToPhysicalMapping = CircularToCartesian<R, Theta, X, Y>;
    static constexpr FootFindingSpace FFSpace = FootFindingSpace::PHYSICAL;
    static constexpr AdvectionFieldSpace AFSpace = AdvectionFieldSpace::PHYSICAL;
    using AdvectionBasis = CartBasis;
};

struct AnalyticalCzarny
{
    using LogicalToPhysicalMapping = CzarnyToCartesian<R, Theta, X, Y>;
    static constexpr FootFindingSpace FFSpace = FootFindingSpace::PHYSICAL;
    static constexpr AdvectionFieldSpace AFSpace = AdvectionFieldSpace::PHYSICAL;
    using AdvectionBasis = CartBasis;
};

struct PseudoCartCzarny
{
    using LogicalToPhysicalMapping = CzarnyToCartesian<R, Theta, X, Y>;
    static constexpr FootFindingSpace FFSpace = FootFindingSpace::PSEUDO_PHYSICAL;
    static constexpr AdvectionFieldSpace AFSpace = AdvectionFieldSpace::PHYSICAL;
    using AdvectionBasis = CartBasis;
};

struct CircularLogical
{
    using LogicalToPhysicalMapping = CircularToCartesian<R, Theta, X, Y>;
    static constexpr FootFindingSpace FFSpace = FootFindingSpace::LOGICAL;
    static constexpr AdvectionFieldSpace AFSpace = AdvectionFieldSpace::LOGICAL;
    using AdvectionBasis = LogicalBasis;
};

struct PseudoCartCircularLogical
{
    using LogicalToPhysicalMapping = CircularToCartesian<R, Theta, X, Y>;
    static constexpr FootFindingSpace FFSpace = FootFindingSpace::PSEUDO_PHYSICAL;
    static constexpr AdvectionFieldSpace AFSpace = AdvectionFieldSpace::LOGICAL;
    using AdvectionBasis = LogicalBasis;
};

template <class T>
struct PolarAdvectionFixture;

template <class TimeStepperBuilderType, class Mappings, class AdvectionFieldType>
struct PolarAdvectionFixture<std::tuple<TimeStepperBuilderType, Mappings, AdvectionFieldType>>
    : public testing::Test
{
    using LogicalToPhysicalMapping = typename Mappings::LogicalToPhysicalMapping;
    using TimeStepperBuilder = TimeStepperBuilderType;
    using AdvectionField = AdvectionFieldType;
    static constexpr FootFindingSpace FFSpace = Mappings::FFSpace;
    static constexpr AdvectionFieldSpace AFSpace = Mappings::AFSpace;
    using AdvectionBasis = typename Mappings::AdvectionBasis;
};

template <class LogicalToOtherMapping>
LogicalToOtherMapping init_mapping()
{
    using OtherX = typename LogicalToOtherMapping::cartesian_tag_x;
    using OtherY = typename LogicalToOtherMapping::cartesian_tag_y;
    // At x0,y0 to match rotation centre
    Coord<OtherX, OtherY> origin_point(0.0, 0.0);
    if constexpr (std::is_same_v<
                          LogicalToOtherMapping,
                          CircularToCartesian<R, Theta, OtherX, OtherY>>) {
        return LogicalToOtherMapping(origin_point);
    } else if constexpr (std::is_same_v<
                                 LogicalToOtherMapping,
                                 CzarnyToCartesian<R, Theta, OtherX, OtherY>>) {
        return LogicalToOtherMapping(0.3, 1.4, origin_point);
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

template <class LogicalToPhysicalMapping, class AdvectionField, class Basis>
void fill_feet_and_advection_field(
        Field<CoordRTheta, IdxRangeSpRTheta> feet,
        Field<CoordRTheta, IdxRangeSpRTheta> exact_feet,
        DVectorField<IdxRangeSpRTheta, Basis> adv_field,
        IdxRangeSpRTheta batched_idx_range,
        LogicalToPhysicalMapping const& to_physical,
        AdvectionField const& advection_field,
        double t,
        double dt)
{
    inverse_mapping_t<LogicalToPhysicalMapping> from_physical = to_physical.get_inverse_mapping();
    // This function is required for GPU compilation
    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            batched_idx_range,
            KOKKOS_LAMBDA(IdxSpRTheta idx) {
                IdxRTheta idx_rtheta(idx);
                CoordRTheta coord_rtheta = ddc::coordinate(idx_rtheta);
                CoordXY coord_xy = to_physical(coord_rtheta);
                ddcHelper::assign_vector_field_element(
                        adv_field,
                        idx,
                        to_vector_space<
                                Basis>(from_physical, coord_xy, advection_field(coord_xy, t)));
                feet(idx) = coord_rtheta;
                exact_feet(idx) = from_physical(advection_field.exact_feet(coord_xy, dt));
            });
}

double calculate_periodic_error(
        Field<CoordRTheta, IdxRangeSpRTheta> feet,
        Field<CoordRTheta, IdxRangeSpRTheta> exact_feet)
{
    const std::source_location location = std::source_location::current();
    return ddc::parallel_transform_reduce(
            location.function_name(),
            Kokkos::DefaultExecutionSpace(),
            get_idx_range(feet),
            0.,
            ddc::reducer::max<double>(),
            KOKKOS_LAMBDA(IdxSpRTheta const idx) {
                Coord<R, Theta> error = feet(idx) - exact_feet(idx);
                if (ddc::get<Theta>(error) < -M_PI || ddc::get<Theta>(error) > M_PI) {
                    double theta_error = Kokkos::fmod(ddc::get<Theta>(error), 2 * M_PI);
                    // Put error on domain [-pi, pi]
                    if (theta_error > M_PI) {
                        theta_error -= 2 * M_PI;
                    }
                    if (theta_error < -M_PI) {
                        theta_error += 2 * M_PI;
                    }
                    error = Coord<R, Theta>(ddc::get<R>(error), theta_error);
                }
                return ::norm_inf(error);
            });
}

using TimeSteppers = std::tuple<RK4Builder>;
using Mappings = std::tuple<
        AnalyticalCircular,
        AnalyticalCzarny,
        PseudoCartCzarny,
        CircularLogical,
        PseudoCartCircularLogical>;
using AdvectionFieldTypes = std::tuple<
        AdvectionField_translation<X, Y>,
        AdvectionField_rotation<X, Y, R, Theta>,
        AdvectionField_decentred_rotation<X, Y>>;

using Cases = tuple_to_types_t<cartesian_product_t<TimeSteppers, Mappings, AdvectionFieldTypes>>;


TYPED_TEST_SUITE(PolarAdvectionFixture, Cases);

TYPED_TEST(PolarAdvectionFixture, Analytical)
{
    using LogicalToPhysicalMapping = typename TestFixture::LogicalToPhysicalMapping;
    using TimeStepperBuilder = typename TestFixture::TimeStepperBuilder;
    using AdvectionField = typename TestFixture::AdvectionField;

    Coord<R> const r_min(0.0);
    Coord<R> const r_max(1.0);
    IdxStep<GridR> nr_cells(25);

    Coord<Theta> const theta_min(0.0);
    Coord<Theta> const theta_max(2.0 * M_PI);
    IdxStep<GridTheta> ntheta_cells(50);

    IdxStepSp const nb_kinspecies(2);
    IdxRangeSp const idx_range_sp(IdxSp(0), nb_kinspecies);

    ddc::init_discrete_space<BSplinesR>(
            build_random_non_uniform_break_points(r_min, r_max, nr_cells, 0.5));
    ddc::init_discrete_space<BSplinesTheta>(
            build_random_non_uniform_break_points(theta_min, theta_max, ntheta_cells, 0.5));

    if (TestFixture::AFSpace == AdvectionFieldSpace::LOGICAL) {
        std::vector<Coord<R>> radial_points;
        NonUniformGridBase<R>::Impl<GridR, Kokkos::HostSpace> greville_sampling_r
                = SplineInterpPointsR::template get_sampling<GridR>();
        Idx<GridR> ir = greville_sampling_r.front();
        for (int i(0); i < greville_sampling_r.size(); ++i, ++ir) {
            radial_points.push_back(greville_sampling_r.coordinate(ir));
        }
        radial_points[0] = Coord<R>(radial_points[1] * 0.5);
        ddc::init_discrete_space<GridR>(radial_points);
    } else {
        ddc::init_discrete_space<GridR>(SplineInterpPointsR::template get_sampling<GridR>());
    }
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

    TimeStepperBuilder time_stepper;
    AdvectionField advection_field = init_field<AdvectionField>();

    PolarFootFinder const batched_foot_finder = make_polar_foot_finder<
            TestFixture::FFSpace,
            TestFixture::AFSpace>(time_stepper, to_physical, batched_idx_range, builder, evaluator);

    const double t = 0.0;
    double dt = 0.001;

    if ((TestFixture::AFSpace == AdvectionFieldSpace::LOGICAL
         or TestFixture::FFSpace == FootFindingSpace::LOGICAL)) {
        if (std::is_same_v<AdvectionField, AdvectionField_translation<X, Y>>) {
            dt /= 8;
        } else if (std::is_same_v<AdvectionField, AdvectionField_decentred_rotation<X, Y>>) {
            dt /= 32;
        }
    }

    DVectorFieldMem<IdxRangeSpRTheta, typename TestFixture::AdvectionBasis> adv_field_alloc(
            batched_idx_range);
    FieldMem<CoordRTheta, IdxRangeSpRTheta> feet_alloc(batched_idx_range);
    FieldMem<CoordRTheta, IdxRangeSpRTheta> exact_feet_alloc(batched_idx_range);

    DVectorField<IdxRangeSpRTheta, typename TestFixture::AdvectionBasis> adv_field
            = get_field(adv_field_alloc);
    Field<CoordRTheta, IdxRangeSpRTheta> feet = get_field(feet_alloc);
    Field<CoordRTheta, IdxRangeSpRTheta> exact_feet = get_field(exact_feet_alloc);

    fill_feet_and_advection_field(
            feet,
            exact_feet,
            adv_field,
            batched_idx_range,
            to_physical,
            advection_field,
            t,
            dt);

    batched_foot_finder(feet, adv_field, dt);

    double error = calculate_periodic_error(feet, exact_feet);

    double TOL = 1e-4;
    EXPECT_NEAR(error, 0.0, TOL);
}

}; // namespace
