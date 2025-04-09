// SPDX-License-Identifier: MIT
#include <cstdlib>
#include <ctime>

#include <ddc/ddc.hpp>

#include <gtest/gtest.h>

#include "circular_to_cartesian.hpp"
#include "czarny_to_cartesian.hpp"
#include "cartesian_to_circular.hpp"
#include "cartesian_to_czarny.hpp"
#include "geometry_mapping_tests.hpp"

template <class Mapping>
void test_analytical_inverse(Mapping map)
{
    static_assert(is_analytical_mapping_v<Mapping>);
    using InverseMapping = inverse_mapping_t<Mapping>;
    InverseMapping inv_map = map.get_inverse_mapping();
    std::srand(std::time(nullptr)); // Seed with random value (the time)

    // Choose a theta value with no chance of a modulo issue
    CoordRTheta
    coord_start(double(rand()) / RAND_MAX, double(rand()) / RAND_MAX * M_PI + M_PI * 0.5);
    CoordXY coord_xy = map(coord_start);
    CoordRTheta coord_end = inv_map(coord_xy);
    EXPECT_NEAR(ddc::get<R>(coord_start), ddc::get<R>(coord_end), 1e-14);
    EXPECT_NEAR(ddc::get<Theta>(coord_start), ddc::get<Theta>(coord_end), 1e-14);
}

template <class Mapping>
void test_analytical_inverse_jacobian(Mapping map)
{
    static_assert(is_analytical_mapping_v<Mapping>);
    using InverseMapping = inverse_mapping_t<Mapping>;
    InverseMapping inv_map = map.get_inverse_mapping();
    std::srand(std::time(nullptr)); // Seed with random value (the time)

    CoordRTheta coord_rtheta(double(rand()) / RAND_MAX, double(rand()) / RAND_MAX * 2.0 * M_PI);
    CoordXY coord_xy = map(coord_rtheta);

    check_inverse_tensor(map.jacobian_matrix(coord_rtheta), inv_map.jacobian_matrix(coord_xy), 1e-14);
    EXPECT_NEAR(map.jacobian(coord_rtheta) * inv_map.jacobian(coord_xy), 1.0, 1e-14);
}

TEST(AnalyticalMappingTests, Circular)
{
    CircularToCartesian<R, Theta, X, Y> mapping;
    test_analytical_inverse(mapping);
    test_analytical_inverse_jacobian(mapping);
}

TEST(AnalyticalMappingTests, Czarny)
{
    CzarnyToCartesian<R, Theta, X, Y> mapping(0.3, 0.4);
    test_analytical_inverse(mapping);
}
