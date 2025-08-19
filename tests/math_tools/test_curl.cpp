#include <gtest/gtest.h>

#include "../coord_transformations/geometry_coord_transformations_tests.hpp"

#include "circular_to_cartesian.hpp"
#include "combined_mapping.hpp"
#include "curl.hpp"
#include "cylindrical_to_cartesian.hpp"
#include "indexed_tensor.hpp"
#include "mesh_builder.hpp"
#include "metric_tensor_evaluator.hpp"
#include "tensor.hpp"
#include "toroidal_to_cylindrical.hpp"
#include "vector_field.hpp"

// Try to confirm that \curl (\grad psi) == 0
// psi = a * sinh(r) * cos(theta)
// B = \grad psi
// \curl B = \curl (\grad psi) = 0
TEST(CurlParamTests, irrotational_field_has_zero_curl)
{
    using Mapping2D = CircularToCartesian<Rho, Theta, R, Z>;
    using ToroidalMapping = ToroidalToCylindrical<Mapping2D, Zeta, Phi>;
    using CylindricalMapping = CylindricalToCartesian<R, Z, Zeta, X, Y>;

    double R0 = 6.2;
    Coord<R, Z> origin_point(R0, 0.0);
    Mapping2D polar_to_RZ(origin_point);
    ToroidalMapping toroidal_to_cylindrical(polar_to_RZ);
    CylindricalMapping cylindrical_to_cartesian;
    CombinedMapping<CylindricalMapping, ToroidalMapping, Coord<Rho, Theta, Phi>>
            mapping(cylindrical_to_cartesian, toroidal_to_cylindrical);
    Curl curl(mapping);

    using BasisSpatial = VectorIndexSet<Rho, Theta, Phi>;
    using CovBasisSpatial = get_covariant_dims_t<BasisSpatial>;
    IdxStepRho nrho(10);
    IdxStepTheta ntheta(10);
    IdxStepPhi nphi(2);
    ddc::init_discrete_space<GridRho>(
            GridRho::init<GridRho>(Coord<Rho>(0.0), Coord<Rho>(1.0), nrho));
    ddc::init_discrete_space<GridTheta>(GridTheta::init<GridTheta>(
            build_uniform_break_points(Coord<Theta>(0.0), Coord<Theta>(2 * M_PI), ntheta)));
    ddc::init_discrete_space<GridPhi>(
            GridPhi::init<GridPhi>(Coord<Phi>(0.0), Coord<Phi>(2 * M_PI), nphi));

    Coord<Rho, Theta, Phi> coord
            = ddc::coordinate(Idx<GridRho, GridTheta, GridPhi>(nrho - 1, ntheta - 5, 0));
    DTensor<CovBasisSpatial, CovBasisSpatial> grad_B;

    double const r = ddc::get<Rho>(coord);
    double const theta = ddc::get<Theta>(coord);
    // Compute dB_r/dr, dB_theta/dr, dB_r/dtheta, dB_theta/dtheta
    ddcHelper::get<R_cov, R_cov>(grad_B)
            = (std::sinh(r) * std::cosh(r) - std::sinh(r) * std::cos(theta)
               + std::sinh(r) * std::cos(theta) * std::cosh(r))
              * std::cos(theta);
    ddcHelper::get<R_cov, Theta_cov>(grad_B)
            = (-std::sinh(r) * std::sinh(r) - std::cosh(r) * std::cosh(r)
               + std::cosh(r) * std::cos(theta))
              * std::sin(theta);
    ddcHelper::get<R_cov, Phi_cov>(grad_B) = 0.0;
    ddcHelper::get<Theta_cov, R_cov>(grad_B)
            = -std::sin(theta) * std::cosh(r) * (std::cosh(r) - 2.0 * std::cos(theta));
    ddcHelper::get<Theta_cov, Theta_cov>(grad_B)
            = -std::sinh(r)
              * (std::sin(theta) * std::sin(theta)
                 + std::cos(theta) * (std::cos(theta) - std::cosh(r)));
    ddcHelper::get<Theta_cov, Phi_cov>(grad_B) = 0.0;
    ddcHelper::get<Phi_cov, R_cov>(grad_B) = 0.0;
    ddcHelper::get<Phi_cov, Theta_cov>(grad_B) = 0.0;
    ddcHelper::get<Phi_cov, Phi_cov>(grad_B) = 0.0;

    DVector<Rho, Theta, Phi> curl_B = curl(grad_B, coord);

    double eps = 1e-12;
    EXPECT_NEAR(ddcHelper::get<Rho>(curl_B), 0.0, eps);
    EXPECT_NEAR(ddcHelper::get<Theta>(curl_B), 0.0, eps);
    EXPECT_NEAR(ddcHelper::get<Phi>(curl_B), 0.0, eps);
}
