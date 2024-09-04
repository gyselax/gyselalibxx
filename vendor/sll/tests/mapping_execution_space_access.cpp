#include <array>
#include <cassert>

#include <ddc/ddc.hpp>

#include "sll/mapping/analytical_invertible_curvilinear2d_to_cartesian.hpp"
#include "sll/mapping/circular_to_cartesian.hpp"
#include "sll/mapping/curvilinear2d_to_cartesian.hpp"
#include "sll/mapping/czarny_to_cartesian.hpp"
#include "sll/mapping/discrete_mapping_builder.hpp"
#include "sll/mapping/discrete_to_cartesian.hpp"

#include "test_utils.hpp"



namespace {
struct X
{
};
struct Y
{
};
struct R
{
    static bool constexpr PERIODIC = false;
};

struct Theta
{
    static bool constexpr PERIODIC = true;
};

using CoordR = ddc::Coordinate<R>;
using CoordTheta = ddc::Coordinate<Theta>;
using CoordRTheta = ddc::Coordinate<R, Theta>;

using CoordXY = ddc::Coordinate<X, Y>;

int constexpr BSDegree = 3;

struct BSplinesR : ddc::NonUniformBSplines<R, BSDegree>
{
};
struct BSplinesTheta : ddc::NonUniformBSplines<Theta, BSDegree>
{
};

using InterpPointsR = ddc::
        GrevilleInterpolationPoints<BSplinesR, ddc::BoundCond::GREVILLE, ddc::BoundCond::GREVILLE>;
using InterpPointsTheta = ddc::GrevilleInterpolationPoints<
        BSplinesTheta,
        ddc::BoundCond::PERIODIC,
        ddc::BoundCond::PERIODIC>;

struct GridR : InterpPointsR::interpolation_discrete_dimension_type
{
};
struct GridTheta : InterpPointsTheta::interpolation_discrete_dimension_type
{
};

using HostExecSpace = Kokkos::DefaultHostExecutionSpace;
using DeviceExecSpace = Kokkos::DefaultExecutionSpace;

template <class ExecSpace>
using SplineRThetaBuilder = ddc::SplineBuilder2D<
        ExecSpace,
        typename ExecSpace::memory_space,
        BSplinesR,
        BSplinesTheta,
        GridR,
        GridTheta,
        ddc::BoundCond::GREVILLE,
        ddc::BoundCond::GREVILLE,
        ddc::BoundCond::PERIODIC,
        ddc::BoundCond::PERIODIC,
        ddc::SplineSolver::LAPACK,
        GridR,
        GridTheta>;

template <class ExecSpace>
using SplineRThetaEvaluator = ddc::SplineEvaluator2D<
        ExecSpace,
        typename ExecSpace::memory_space,
        BSplinesR,
        BSplinesTheta,
        GridR,
        GridTheta,
        ddc::NullExtrapolationRule,
        ddc::NullExtrapolationRule,
        ddc::PeriodicExtrapolationRule<Theta>,
        ddc::PeriodicExtrapolationRule<Theta>,
        GridR,
        GridTheta>;

using IdxRangeBSR = ddc::DiscreteDomain<BSplinesR>;
using IdxRangeBSTheta = ddc::DiscreteDomain<BSplinesTheta>;
using IdxRangeBSRTheta = ddc::DiscreteDomain<BSplinesR, BSplinesTheta>;

using IdxRangeR = ddc::DiscreteDomain<GridR>;
using IdxRangeTheta = ddc::DiscreteDomain<GridTheta>;
using IdxRangeRTheta = ddc::DiscreteDomain<GridR, GridTheta>;

using IdxR = ddc::DiscreteElement<GridR>;
using IdxTheta = ddc::DiscreteElement<GridTheta>;
using IdxRTheta = ddc::DiscreteElement<GridR, GridTheta>;

using IdxStepR = ddc::DiscreteVector<GridR>;
using IdxStepTheta = ddc::DiscreteVector<GridTheta>;
using IdxStepRTheta = ddc::DiscreteVector<GridR, GridTheta>;

using IdxRangeRTheta = ddc::DiscreteDomain<GridR, GridTheta>;


class MappingMemoryAccess : public ::testing::Test
{
protected:
    static int constexpr npts_r = 32;
    static int constexpr npts_theta = 64;

    static constexpr CoordR r_min = CoordR(0.0);
    static constexpr CoordR r_max = CoordR(1.0);

    static constexpr CoordTheta theta_min = CoordTheta(0.0);
    static constexpr CoordTheta theta_max = CoordTheta(2.0 * M_PI);

    IdxRangeR const interpolation_idx_range_r;
    IdxRangeTheta const interpolation_idx_range_theta;
    IdxRangeRTheta const interpolation_idx_range_rtheta;

public:
    MappingMemoryAccess()
        : interpolation_idx_range_r(InterpPointsR::get_domain<GridR>())
        , interpolation_idx_range_theta(InterpPointsTheta::get_domain<GridTheta>())
        , interpolation_idx_range_rtheta(
                  interpolation_idx_range_r,
                  interpolation_idx_range_theta) {};

    static void SetUpTestSuite()
    {
        double const dr((r_max - r_min) / npts_r);
        double const dtheta((theta_max - theta_min) / npts_theta);

        std::vector<CoordR> r_break_points(npts_r + 1);
        std::vector<CoordTheta> theta_break_points(npts_theta + 1);

        for (int i(0); i < npts_r + 1; ++i) {
            r_break_points[i] = CoordR(r_min + i * dr);
        }
        r_break_points[npts_r] = CoordR(r_max);
        for (int i(0); i < npts_theta + 1; ++i) {
            theta_break_points[i] = CoordTheta(theta_min + i * dtheta);
        }

        ddc::init_discrete_space<BSplinesR>(r_break_points);
        ddc::init_discrete_space<BSplinesTheta>(theta_break_points);

        ddc::init_discrete_space<GridR>(InterpPointsR::get_sampling<GridR>());
        ddc::init_discrete_space<GridTheta>(InterpPointsTheta::get_sampling<GridTheta>());
    }
};

template <class Mapping>
double check_logical_to_physical_coord_converter(
        Mapping const& mapping,
        CoordRTheta const& coord_rtheta,
        CoordXY const& coord_xy)
{
    DeviceExecSpace exec;

    double max_error = 0;
    Kokkos::parallel_reduce(
            Kokkos::RangePolicy<DeviceExecSpace>(exec, 0, 1),
            KOKKOS_LAMBDA(int const i, double& err) {
                // Coord converter: logical -> physical.
                CoordXY diff_coord_xy = mapping(coord_rtheta) - coord_xy;
                err = Kokkos::max(ddc::get<X>(diff_coord_xy), err);
                err = Kokkos::max(ddc::get<Y>(diff_coord_xy), err);
            },
            Kokkos::Max<double>(max_error));
    return max_error;
};

template <class Mapping>
double check_physical_to_logical_coord_converter(
        Mapping const& mapping,
        CoordRTheta const& coord_rtheta,
        CoordXY const& coord_xy)
{
    DeviceExecSpace exec;

    double max_error = 0;
    Kokkos::parallel_reduce(
            Kokkos::RangePolicy<DeviceExecSpace>(exec, 0, 1),
            KOKKOS_LAMBDA(int const i, double& err) {
                // Coord converter: physical -> logical.
                CoordRTheta diff_coord_rtheta = mapping(coord_xy) - coord_rtheta;
                err = Kokkos::max(ddc::get<R>(diff_coord_rtheta), err);
                err = Kokkos::max(ddc::get<Theta>(diff_coord_rtheta), err);
            },
            Kokkos::Max<double>(max_error));
    return max_error;
};

} // namespace



TEST_F(MappingMemoryAccess, HostCircularCoordConverter)
{
    // Mapping
    CircularToCartesian<X, Y, R, Theta> const mapping;

    // Test coordinates
    double const r = .75;
    double const theta = 1 / 3. * M_PI;
    CoordRTheta const coord_rtheta(r, theta);

    double const x = r * Kokkos::cos(theta);
    double const y = r * Kokkos::sin(theta);
    CoordXY const coord_xy(x, y);

    // Coord converter: logical -> physical.
    CoordXY diff_coord_xy = mapping(coord_rtheta) - coord_xy;
    EXPECT_LE(ddc::get<X>(diff_coord_xy), 1e-15);
    EXPECT_LE(ddc::get<Y>(diff_coord_xy), 1e-15);

    // Coord converter: physical--> logical.
    CoordRTheta diff_coord_rtheta = mapping(coord_xy) - coord_rtheta;
    EXPECT_LE(ddc::get<R>(diff_coord_rtheta), 1e-15);
    EXPECT_LE(ddc::get<Theta>(diff_coord_rtheta), 1e-15);
}


TEST_F(MappingMemoryAccess, HostCzarnyCoordConverter)
{
    // Mapping
    double const epsilon = 0.3;
    double const e = 1.4;
    CzarnyToCartesian<X, Y, R, Theta> const mapping(epsilon, e);

    // Test coordinates
    double const r = .75;
    double const theta = 1 / 3. * M_PI;
    CoordRTheta const coord_rtheta(r, theta);

    const double tmp1 = Kokkos::sqrt(epsilon * (epsilon + 2.0 * r * Kokkos::cos(theta)) + 1.0);
    const double x = (1.0 - tmp1) / epsilon;
    const double y = e * r * Kokkos::sin(theta)
                     / (Kokkos::sqrt(1.0 - 0.25 * epsilon * epsilon) * (2.0 - tmp1));
    CoordXY const coord_xy(x, y);

    // Coord converter: logical -> physical.
    CoordXY diff_coord_xy = mapping(coord_rtheta) - coord_xy;
    EXPECT_LE(ddc::get<X>(diff_coord_xy), 1e-15);
    EXPECT_LE(ddc::get<Y>(diff_coord_xy), 1e-15);

    // Coord converter: physical -> logical.
    CoordRTheta diff_coord_rtheta = mapping(coord_xy) - coord_rtheta;
    EXPECT_LE(ddc::get<R>(diff_coord_rtheta), 1e-15);
    EXPECT_LE(ddc::get<Theta>(diff_coord_rtheta), 1e-15);
}


TEST_F(MappingMemoryAccess, HostDiscreteCoordConverter)
{
    // Mapping
    double const epsilon = 0.3;
    double const e = 1.4;
    CzarnyToCartesian<X, Y, R, Theta> const analytical_mapping(epsilon, e);

    SplineRThetaBuilder<HostExecSpace> builder(interpolation_idx_range_rtheta);

    ddc::NullExtrapolationRule r_extrapolation_rule;
    ddc::PeriodicExtrapolationRule<Theta> theta_extrapolation_rule;
    SplineRThetaEvaluator<HostExecSpace> evaluator(
            r_extrapolation_rule,
            r_extrapolation_rule,
            theta_extrapolation_rule,
            theta_extrapolation_rule);

    DiscreteToCartesianBuilder<
            X,
            Y,
            SplineRThetaBuilder<HostExecSpace>,
            SplineRThetaEvaluator<HostExecSpace>>
            mapping_builder(HostExecSpace(), analytical_mapping, builder, evaluator);
    DiscreteToCartesian mapping = mapping_builder();

    // Test coordinates
    double const r = .75;
    double const theta = 1 / 3. * M_PI;
    CoordRTheta const coord_rtheta(r, theta);

    const double tmp1 = Kokkos::sqrt(epsilon * (epsilon + 2.0 * r * Kokkos::cos(theta)) + 1.0);
    const double x = (1.0 - tmp1) / epsilon;
    const double y = e * r * Kokkos::sin(theta)
                     / (Kokkos::sqrt(1.0 - 0.25 * epsilon * epsilon) * (2.0 - tmp1));
    CoordXY const coord_xy(x, y);

    // Coord converter: logical -> physical.
    CoordXY diff_coord_xy = mapping(coord_rtheta) - coord_xy;
    EXPECT_LE(ddc::get<X>(diff_coord_xy), 1e-7);
    EXPECT_LE(ddc::get<Y>(diff_coord_xy), 1e-7);
}



TEST_F(MappingMemoryAccess, DeviceCircularCoordConverter)
{
    // Mapping
    CircularToCartesian<X, Y, R, Theta> const mapping;

    // Test coordinates
    double const r = .75;
    double const theta = 1 / 3. * M_PI;
    CoordRTheta const coord_rtheta(r, theta);

    double const x = r * Kokkos::cos(theta);
    double const y = r * Kokkos::sin(theta);
    CoordXY const coord_xy(x, y);

    // Coord converter: logical -> physical.
    double err = check_logical_to_physical_coord_converter(mapping, coord_rtheta, coord_xy);
    EXPECT_LE(err, 1e-15);

    // Coord converter: physical--> logical.
    err = check_physical_to_logical_coord_converter(mapping, coord_rtheta, coord_xy);
    EXPECT_LE(err, 1e-15);
}


TEST_F(MappingMemoryAccess, DeviceCzarnyCoordConverter)
{
    // Mapping
    double const epsilon = 0.3;
    double const e = 1.4;
    CzarnyToCartesian<X, Y, R, Theta> const mapping(epsilon, e);

    // Test coordinates
    double const r = .75;
    double const theta = 1 / 3. * M_PI;
    CoordRTheta const coord_rtheta(r, theta);

    const double tmp1 = Kokkos::sqrt(epsilon * (epsilon + 2.0 * r * Kokkos::cos(theta)) + 1.0);
    const double x = (1.0 - tmp1) / epsilon;
    const double y = e * r * Kokkos::sin(theta)
                     / (Kokkos::sqrt(1.0 - 0.25 * epsilon * epsilon) * (2.0 - tmp1));
    CoordXY const coord_xy(x, y);

    // Coord converter: logical -> physical.
    double err = check_logical_to_physical_coord_converter(mapping, coord_rtheta, coord_xy);
    EXPECT_LE(err, 1e-15);

    // Coord converter: physical--> logical.
    err = check_physical_to_logical_coord_converter(mapping, coord_rtheta, coord_xy);
    EXPECT_LE(err, 1e-15);
}


TEST_F(MappingMemoryAccess, DeviceDiscreteCoordConverter)
{
    // Mapping
    double const epsilon = 0.3;
    double const e = 1.4;
    CzarnyToCartesian<X, Y, R, Theta> const analytical_mapping(epsilon, e);

    SplineRThetaBuilder<DeviceExecSpace> builder(interpolation_idx_range_rtheta);

    ddc::NullExtrapolationRule r_extrapolation_rule;
    ddc::PeriodicExtrapolationRule<Theta> theta_extrapolation_rule;
    SplineRThetaEvaluator<DeviceExecSpace> evaluator(
            r_extrapolation_rule,
            r_extrapolation_rule,
            theta_extrapolation_rule,
            theta_extrapolation_rule);

    DiscreteToCartesianBuilder<
            X,
            Y,
            SplineRThetaBuilder<DeviceExecSpace>,
            SplineRThetaEvaluator<DeviceExecSpace>>
            mapping_builder(DeviceExecSpace(), analytical_mapping, builder, evaluator);
    DiscreteToCartesian mapping = mapping_builder();

    // Test coordinates
    double const r = .75;
    double const theta = 1 / 3. * M_PI;
    CoordRTheta const coord_rtheta(r, theta);

    const double tmp1 = Kokkos::sqrt(epsilon * (epsilon + 2.0 * r * Kokkos::cos(theta)) + 1.0);
    const double x = (1.0 - tmp1) / epsilon;
    const double y = e * r * Kokkos::sin(theta)
                     / (Kokkos::sqrt(1.0 - 0.25 * epsilon * epsilon) * (2.0 - tmp1));
    CoordXY const coord_xy(x, y);

    // Coord converter: logical -> physical.
    double err = check_logical_to_physical_coord_converter(mapping, coord_rtheta, coord_xy);
    EXPECT_LE(err, 1e-7);
}
