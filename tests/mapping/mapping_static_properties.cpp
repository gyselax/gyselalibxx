#include <gtest/gtest.h>

#include "barycentric_to_cartesian.hpp"
#include "cartesian_to_barycentric.hpp"
#include "cartesian_to_circular.hpp"
#include "cartesian_to_cylindrical.hpp"
#include "cartesian_to_czarny.hpp"
#include "circular_to_cartesian.hpp"
#include "combined_mapping.hpp"
#include "cylindrical_to_cartesian.hpp"
#include "czarny_to_cartesian.hpp"
#include "discrete_to_cartesian.hpp"
#include "geometry_mapping_tests.hpp"
#include "identity_coordinate_change.hpp"
#include "mapping_tools.hpp"
#include "toroidal_to_cylindrical.hpp"


TEST(MappingStaticAsserts, CirctoCart)
{
    static_assert(is_mapping_v<CircularToCartesian<R, Theta, X, Y>>);
    static_assert(has_jacobian_v<CircularToCartesian<R, Theta, X, Y>>);
    static_assert(has_inv_jacobian_v<CircularToCartesian<R, Theta, X, Y>>);
    static_assert(is_curvilinear_2d_mapping_v<CircularToCartesian<R, Theta, X, Y>>);
    static_assert(is_analytical_mapping_v<CircularToCartesian<R, Theta, X, Y>>);
    static_assert(has_singular_o_point_inv_jacobian_v<CircularToCartesian<R, Theta, X, Y>>);
}

TEST(MappingStaticAsserts, CartToCirc)
{
    static_assert(is_mapping_v<CartesianToCircular<X, Y, R, Theta>>);
    static_assert(has_jacobian_v<CartesianToCircular<X, Y, R, Theta>>);
    static_assert(is_analytical_mapping_v<CartesianToCircular<X, Y, R, Theta>>);
}


TEST(MappingStaticAsserts, CzarnytoCart)
{
    static_assert(is_mapping_v<CzarnyToCartesian<R, Theta, X, Y>>);
    static_assert(has_jacobian_v<CzarnyToCartesian<R, Theta, X, Y>>);
    static_assert(has_inv_jacobian_v<CzarnyToCartesian<R, Theta, X, Y>>);
    static_assert(is_curvilinear_2d_mapping_v<CzarnyToCartesian<R, Theta, X, Y>>);
    static_assert(is_analytical_mapping_v<CzarnyToCartesian<R, Theta, X, Y>>);
    static_assert(has_singular_o_point_inv_jacobian_v<CzarnyToCartesian<R, Theta, X, Y>>);
}

TEST(MappingStaticAsserts, CartToCzarny)
{
    static_assert(is_mapping_v<CartesianToCzarny<X, Y, R, Theta>>);
    static_assert(is_analytical_mapping_v<CartesianToCzarny<X, Y, R, Theta>>);
}

TEST(MappingStaticAsserts, DiscToCart)
{
    using Mapping = DiscreteToCartesian<X, Y, SplineRThetaEvaluator_host>;
    static_assert(is_mapping_v<Mapping>);
    static_assert(has_jacobian_v<Mapping>);
    static_assert(is_curvilinear_2d_mapping_v<Mapping>);
    static_assert(!is_analytical_mapping_v<Mapping>);
    static_assert(has_singular_o_point_inv_jacobian_v<Mapping>);
}

TEST(MappingStaticAsserts, CombinedMapping)
{
    using Mapping = CombinedMapping<
            CircularToCartesian<R, Theta, X, Y>,
            CartesianToCzarny<X, Y, R, Theta>>;
    static_assert(is_mapping_v<Mapping>);
    static_assert(has_jacobian_v<Mapping>);
    static_assert(has_inv_jacobian_v<Mapping>);
}

TEST(MappingStaticAsserts, CylToCart)
{
    static_assert(is_mapping_v<CylindricalToCartesian<R, Z, Zeta, X, Y>>);
    static_assert(has_jacobian_v<CylindricalToCartesian<R, Z, Zeta, X, Y>>);
    static_assert(has_inv_jacobian_v<CylindricalToCartesian<R, Z, Zeta, X, Y>>);
    static_assert(is_analytical_mapping_v<CylindricalToCartesian<R, Z, Zeta, X, Y>>);
}

TEST(MappingStaticAsserts, CartToCyl)
{
    static_assert(is_mapping_v<CartesianToCylindrical<X, Y, Z, R, Zeta>>);
    static_assert(has_jacobian_v<CartesianToCylindrical<X, Y, Z, R, Zeta>>);
    static_assert(is_analytical_mapping_v<CartesianToCylindrical<X, Y, Z, R, Zeta>>);
}

TEST(MappingStaticAsserts, BarycentricToCart)
{
    struct Corner1;
    struct Corner2;
    struct Corner3;
    using Mapping = BarycentricToCartesian<Corner1, Corner2, Corner3, X, Y>;
    static_assert(is_mapping_v<Mapping>);
    static_assert(is_analytical_mapping_v<Mapping>);
}

TEST(MappingStaticAsserts, CartToBarycentric)
{
    struct Corner1;
    struct Corner2;
    struct Corner3;
    using Mapping = CartesianToBarycentric<X, Y, Corner1, Corner2, Corner3>;
    static_assert(is_mapping_v<Mapping>);
    static_assert(is_analytical_mapping_v<Mapping>);
}

TEST(MappingStaticAsserts, Identity)
{
    using Mapping = IdentityCoordinateChange<VectorIndexSet<X, Y, Z>, VectorIndexSet<X, Y, Z>>;
    static_assert(is_mapping_v<Mapping>);
    static_assert(has_jacobian_v<Mapping>);
    static_assert(has_inv_jacobian_v<Mapping>);
    static_assert(is_analytical_mapping_v<Mapping>);
}

TEST(MappingStaticAsserts, TorToCyl)
{
    using Mapping2D = CircularToCartesian<R, Theta, X, Y>;
    using Mapping = ToroidalToCylindrical<Mapping2D, Zeta, Phi>;
    static_assert(is_mapping_v<Mapping>);
    static_assert(has_jacobian_v<Mapping>);
    static_assert(has_inv_jacobian_v<Mapping>);
    static_assert(has_singular_o_point_inv_jacobian_v<Mapping>);
}
