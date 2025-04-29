#include <gtest/gtest.h>

#include "../mapping/geometry_mapping_tests.hpp"

#include "circular_to_cartesian.hpp"
#include "combined_mapping.hpp"
#include "cylindrical_to_cartesian.hpp"
#include "toroidal_to_cylindrical.hpp"


TEST(GyrokineticPoissonBracket, Anticommutativity)
{
    using Mapping2D = CircularToCartesian<R, Z, Rho, Theta>;
    using ToroidalMapping = ToroidalToCylindrical<Mapping2D, Zeta, Phi>;
    using CylindricalMapping = CylindricalToCartesian<R, Z, Zeta, X, Y>;
    Mapping2D polar_to_RZ;
    ToroidalMapping toroidal_to_cylindrical(polar_to_RZ);
    CylindricalMapping cylindrical_to_cartesian;
    CombinedMapping<CylindricalMapping, ToroidalMapping>
            mapping(toroidal_to_cylindrical, cylindrical_to_cartesian);
}

TEST(GyrokineticPoissonBracket, Bilinearity) {}

TEST(GyrokineticPoissonBracket, Leibniz) {}

TEST(GyrokineticPoissonBracket, Jacobi) {}
