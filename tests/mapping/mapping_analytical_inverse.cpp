// SPDX-License-Identifier: MIT
#include <cstdlib>
#include <ctime>

#include <ddc/ddc.hpp>

#include <gtest/gtest.h>

#include "cartesian_to_circular.hpp"
#include "cartesian_to_cylindrical.hpp"
#include "cartesian_to_czarny.hpp"
#include "circular_to_cartesian.hpp"
#include "cylindrical_to_cartesian.hpp"
#include "czarny_to_cartesian.hpp"
#include "geometry_mapping_tests.hpp"
#include "math_tools.hpp"

namespace {

inline void expect_near(double d1, double d2, double Tol)
{
    EXPECT_NEAR(d1, d2, Tol);
}

template <class... Dims>
void compare_coord(Coord<Dims...> const& coord1, Coord<Dims...> const& coord2, double Tol)
{
    ((expect_near(ddc::get<Dims>(coord1), ddc::get<Dims>(coord2), Tol)), ...);
}

template <class Mapping>
void test_analytical_inverse(Mapping map, typename Mapping::CoordArg coord_start)
{
    static_assert(is_analytical_mapping_v<Mapping>);
    using InverseMapping = inverse_mapping_t<Mapping>;
    using CoordArg = typename Mapping::CoordArg;
    using CoordResult = typename Mapping::CoordResult;
    InverseMapping inv_map = map.get_inverse_mapping();
    std::srand(std::time(nullptr)); // Seed with random value (the time)

    // Choose a theta value with no chance of a modulo issue
    CoordResult coord_xy = map(coord_start);
    CoordArg coord_end = inv_map(coord_xy);
    compare_coord(coord_start, coord_end, 1e-14);
}

template <class Mapping>
void test_analytical_inverse_jacobian(Mapping map, typename Mapping::CoordArg coord_arg)
{
    static_assert(is_analytical_mapping_v<Mapping>);
    using InverseMapping = inverse_mapping_t<Mapping>;
    using CoordResult = typename Mapping::CoordResult;
    InverseMapping inv_map = map.get_inverse_mapping();
    std::srand(std::time(nullptr)); // Seed with random value (the time)

    CoordResult coord_result = map(coord_arg);

    check_inverse_tensor(
            map.jacobian_matrix(coord_arg),
            inv_map.jacobian_matrix(coord_result),
            1e-14);
    EXPECT_NEAR(map.jacobian(coord_arg) * inv_map.jacobian(coord_result), 1.0, 1e-14);
    EXPECT_NEAR(map.jacobian(coord_arg), determinant(map.jacobian_matrix(coord_arg)), 1e-14);
    EXPECT_NEAR(
            inv_map.jacobian(coord_result),
            determinant(inv_map.jacobian_matrix(coord_result)),
            1e-14);
}

} // namespace

TEST(AnalyticalMappingTests, Circular)
{
    CircularToCartesian<R, Theta, X, Y> mapping;
    CoordRTheta coord(double(rand()) / RAND_MAX, double(rand()) / RAND_MAX * M_PI + M_PI * 0.5);
    test_analytical_inverse(mapping, coord);
    test_analytical_inverse_jacobian(mapping, coord);
}

TEST(AnalyticalMappingTests, Czarny)
{
    CzarnyToCartesian<R, Theta, X, Y> mapping(0.3, 0.4);
    CoordRTheta coord(double(rand()) / RAND_MAX, double(rand()) / RAND_MAX * M_PI + M_PI * 0.5);
    test_analytical_inverse(mapping, coord);
}

TEST(AnalyticalMappingTests, Cylindrical)
{
    CylindricalToCartesian<R, Z, Zeta, X, Y> mapping;
    CoordRZZeta
    coord(double(rand()) / RAND_MAX,
          double(rand()) / RAND_MAX,
          double(rand()) / RAND_MAX * M_PI + M_PI * 0.5);
    test_analytical_inverse(mapping, coord);
    test_analytical_inverse_jacobian(mapping, coord);
}
