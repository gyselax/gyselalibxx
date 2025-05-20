#include <array>
#include <cassert>

#include <ddc/ddc.hpp>

#include <gtest/gtest.h>

#include "cartesian_to_circular.hpp"
#include "cartesian_to_czarny.hpp"
#include "circular_to_cartesian.hpp"
#include "czarny_to_cartesian.hpp"
#include "discrete_mapping_builder.hpp"
#include "discrete_to_cartesian.hpp"
#include "geometry_mapping_tests.hpp"



namespace {
using HostExecSpace = Kokkos::DefaultHostExecutionSpace;
using DeviceExecSpace = Kokkos::DefaultExecutionSpace;


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
        Mapping const& to_physical_mapping,
        CoordRTheta const& coord_rtheta,
        CoordXY const& coord_xy)
{
    DeviceExecSpace exec;
    static_assert(is_accessible_v<DeviceExecSpace, Mapping>);

    double max_error = 0;
    Kokkos::parallel_reduce(
            Kokkos::RangePolicy<DeviceExecSpace>(exec, 0, 1),
            KOKKOS_LAMBDA(int const i, double& err) {
                // Coord converter: logical -> physical.
                CoordXY diff_coord_xy = to_physical_mapping(coord_rtheta) - coord_xy;
                err = Kokkos::max(ddc::get<X>(diff_coord_xy), err);
                err = Kokkos::max(ddc::get<Y>(diff_coord_xy), err);
            },
            Kokkos::Max<double>(max_error));
    return max_error;
};

template <class Mapping>
double check_physical_to_logical_coord_converter(
        Mapping const& to_physical_mapping,
        CoordRTheta const& coord_rtheta,
        CoordXY const& coord_xy)
{
    DeviceExecSpace exec;
    static_assert(is_accessible_v<DeviceExecSpace, Mapping>);

    double max_error = 0;
    Kokkos::parallel_reduce(
            Kokkos::RangePolicy<DeviceExecSpace>(exec, 0, 1),
            KOKKOS_LAMBDA(int const i, double& err) {
                // Coord converter: physical -> logical.
                CoordRTheta diff_coord_rtheta = to_physical_mapping(coord_xy) - coord_rtheta;
                err = Kokkos::max(ddc::get<R>(diff_coord_rtheta), err);
                err = Kokkos::max(ddc::get<Theta>(diff_coord_rtheta), err);
            },
            Kokkos::Max<double>(max_error));
    return max_error;
};

} // namespace


TEST_F(MappingMemoryAccess, HostCircularCoordConverter)
{
    static_assert(is_mapping_v<CartesianToCircular<X, Y, R, Theta>>);
    static_assert(is_mapping_v<CircularToCartesian<R, Theta, X, Y>>);
    static_assert(is_analytical_mapping_v<CartesianToCircular<X, Y, R, Theta>>);
    static_assert(std::is_same_v<
                  inverse_mapping_t<CartesianToCircular<X, Y, R, Theta>>,
                  CircularToCartesian<R, Theta, X, Y>>);
    static_assert(is_accessible_v<
                  Kokkos::DefaultHostExecutionSpace,
                  CartesianToCircular<X, Y, R, Theta>>);
    static_assert(is_accessible_v<
                  Kokkos::DefaultHostExecutionSpace,
                  CircularToCartesian<R, Theta, X, Y>>);

    // Mapping
    CircularToCartesian<R, Theta, X, Y> const to_physical_mapping;
    CartesianToCircular<X, Y, R, Theta> const to_logical_mapping
            = to_physical_mapping.get_inverse_mapping();

    // Test coordinates
    double const r = .75;
    double const theta = 1 / 3. * M_PI;
    CoordRTheta const coord_rtheta(r, theta);

    double const x = r * std::cos(theta);
    double const y = r * std::sin(theta);
    CoordXY const coord_xy(x, y);

    // Coord converter: logical -> physical.
    CoordXY diff_coord_xy = to_physical_mapping(coord_rtheta) - coord_xy;
    EXPECT_LE(ddc::get<X>(diff_coord_xy), 1e-15);
    EXPECT_LE(ddc::get<Y>(diff_coord_xy), 1e-15);

    // Coord converter: physical--> logical.
    CoordRTheta diff_coord_rtheta = to_logical_mapping(coord_xy) - coord_rtheta;
    EXPECT_LE(ddc::get<R>(diff_coord_rtheta), 1e-15);
    EXPECT_LE(ddc::get<Theta>(diff_coord_rtheta), 1e-15);
}


TEST_F(MappingMemoryAccess, HostCzarnyCoordConverter)
{
    static_assert(is_mapping_v<CartesianToCzarny<X, Y, R, Theta>>);
    static_assert(is_mapping_v<CzarnyToCartesian<R, Theta, X, Y>>);
    static_assert(is_analytical_mapping_v<CartesianToCzarny<X, Y, R, Theta>>);
    static_assert(std::is_same_v<
                  inverse_mapping_t<CartesianToCzarny<X, Y, R, Theta>>,
                  CzarnyToCartesian<R, Theta, X, Y>>);
    static_assert(
            is_accessible_v<Kokkos::DefaultHostExecutionSpace, CartesianToCzarny<X, Y, R, Theta>>);
    static_assert(
            is_accessible_v<Kokkos::DefaultHostExecutionSpace, CzarnyToCartesian<R, Theta, X, Y>>);

    // Mapping
    double const epsilon = 0.3;
    double const e = 1.4;
    CzarnyToCartesian<R, Theta, X, Y> const to_physical_mapping(epsilon, e);
    CartesianToCzarny<X, Y, R, Theta> const to_logical_mapping
            = to_physical_mapping.get_inverse_mapping();

    // Test coordinates
    double const r = .75;
    double const theta = 1 / 3. * M_PI;
    CoordRTheta const coord_rtheta(r, theta);

    const double tmp1 = std::sqrt(epsilon * (epsilon + 2.0 * r * std::cos(theta)) + 1.0);
    const double x = (1.0 - tmp1) / epsilon;
    const double y
            = e * r * std::sin(theta) / (std::sqrt(1.0 - 0.25 * epsilon * epsilon) * (2.0 - tmp1));
    CoordXY const coord_xy(x, y);

    // Coord converter: logical -> physical.
    CoordXY diff_coord_xy = to_physical_mapping(coord_rtheta) - coord_xy;
    EXPECT_LE(ddc::get<X>(diff_coord_xy), 1e-15);
    EXPECT_LE(ddc::get<Y>(diff_coord_xy), 1e-15);

    // Coord converter: physical -> logical.
    CoordRTheta diff_coord_rtheta = to_logical_mapping(coord_xy) - coord_rtheta;
    EXPECT_LE(ddc::get<R>(diff_coord_rtheta), 1e-15);
    EXPECT_LE(ddc::get<Theta>(diff_coord_rtheta), 1e-15);
}


TEST_F(MappingMemoryAccess, HostDiscreteCoordConverter)
{
    // Mapping
    double const epsilon = 0.3;
    double const e = 1.4;
    CzarnyToCartesian<R, Theta, X, Y> const analytical_mapping(epsilon, e);
    static_assert(
            is_accessible_v<Kokkos::DefaultHostExecutionSpace, CzarnyToCartesian<R, Theta, X, Y>>);

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
    DiscreteToCartesian to_physical_mapping = mapping_builder();
    static_assert(
            is_accessible_v<Kokkos::DefaultHostExecutionSpace, decltype(to_physical_mapping)>);

    // Test coordinates
    double const r = .75;
    double const theta = 1 / 3. * M_PI;
    CoordRTheta const coord_rtheta(r, theta);

    const double tmp1 = std::sqrt(epsilon * (epsilon + 2.0 * r * std::cos(theta)) + 1.0);
    const double x = (1.0 - tmp1) / epsilon;
    const double y
            = e * r * std::sin(theta) / (std::sqrt(1.0 - 0.25 * epsilon * epsilon) * (2.0 - tmp1));
    CoordXY const coord_xy(x, y);

    // Coord converter: logical -> physical.
    CoordXY diff_coord_xy = to_physical_mapping(coord_rtheta) - coord_xy;
    EXPECT_LE(ddc::get<X>(diff_coord_xy), 1e-7);
    EXPECT_LE(ddc::get<Y>(diff_coord_xy), 1e-7);
}


TEST_F(MappingMemoryAccess, DeviceCircularCoordConverter)
{
    // Mapping
    CircularToCartesian<R, Theta, X, Y> const to_physical_mapping;
    CartesianToCircular<X, Y, R, Theta> const to_logical_mapping;
    static_assert(
            is_accessible_v<Kokkos::DefaultExecutionSpace, CircularToCartesian<R, Theta, X, Y>>);

    // Test coordinates
    double const r = .75;
    double const theta = 1 / 3. * M_PI;
    CoordRTheta const coord_rtheta(r, theta);

    double const x = r * std::cos(theta);
    double const y = r * std::sin(theta);
    CoordXY const coord_xy(x, y);

    // Coord converter: logical -> physical.
    double err = check_logical_to_physical_coord_converter(
            to_physical_mapping,
            coord_rtheta,
            coord_xy);
    EXPECT_LE(err, 1e-15);

    // Coord converter: physical--> logical.
    err = check_physical_to_logical_coord_converter(to_logical_mapping, coord_rtheta, coord_xy);
    EXPECT_LE(err, 1e-15);
}


TEST_F(MappingMemoryAccess, DeviceCzarnyCoordConverter)
{
    // Mapping
    double const epsilon = 0.3;
    double const e = 1.4;
    CzarnyToCartesian<R, Theta, X, Y> const to_physical_mapping(epsilon, e);
    CartesianToCzarny<X, Y, R, Theta> const to_logical_mapping(epsilon, e);

    static_assert(
            is_accessible_v<Kokkos::DefaultExecutionSpace, CzarnyToCartesian<R, Theta, X, Y>>);

    // Test coordinates
    double const r = .75;
    double const theta = 1 / 3. * M_PI;
    CoordRTheta const coord_rtheta(r, theta);

    const double tmp1 = std::sqrt(epsilon * (epsilon + 2.0 * r * std::cos(theta)) + 1.0);
    const double x = (1.0 - tmp1) / epsilon;
    const double y
            = e * r * std::sin(theta) / (std::sqrt(1.0 - 0.25 * epsilon * epsilon) * (2.0 - tmp1));
    CoordXY const coord_xy(x, y);

    // Coord converter: logical -> physical.
    double err = check_logical_to_physical_coord_converter(
            to_physical_mapping,
            coord_rtheta,
            coord_xy);
    EXPECT_LE(err, 1e-15);

    // Coord converter: physical--> logical.
    err = check_physical_to_logical_coord_converter(to_logical_mapping, coord_rtheta, coord_xy);
    EXPECT_LE(err, 1e-15);
}


TEST_F(MappingMemoryAccess, DeviceDiscreteCoordConverter)
{
    // Mapping
    double const epsilon = 0.3;
    double const e = 1.4;
    CzarnyToCartesian<R, Theta, X, Y> const analytical_mapping(epsilon, e);
    static_assert(
            is_accessible_v<Kokkos::DefaultExecutionSpace, CzarnyToCartesian<R, Theta, X, Y>>);

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
    DiscreteToCartesian to_physical_mapping = mapping_builder();
    static_assert(is_accessible_v<DeviceExecSpace, decltype(to_physical_mapping)>);

    // Test coordinates
    double const r = .75;
    double const theta = 1 / 3. * M_PI;
    CoordRTheta const coord_rtheta(r, theta);

    const double tmp1 = std::sqrt(epsilon * (epsilon + 2.0 * r * std::cos(theta)) + 1.0);
    const double x = (1.0 - tmp1) / epsilon;
    const double y
            = e * r * std::sin(theta) / (std::sqrt(1.0 - 0.25 * epsilon * epsilon) * (2.0 - tmp1));
    CoordXY const coord_xy(x, y);

    // Coord converter: logical -> physical.
    double err = check_logical_to_physical_coord_converter(
            to_physical_mapping,
            coord_rtheta,
            coord_xy);
    EXPECT_LE(err, 1e-7);
}
