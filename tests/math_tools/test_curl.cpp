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

class CurlParamTests : public ::testing::TestWithParam<std::tuple<double, double>>
{
protected:
    void SetUp() override
    {
        m_rho = std::get<0>(GetParam());
        m_theta = std::get<1>(GetParam());
    }

protected:
    double const m_R0 = 6.2;
    double m_rho;
    double m_theta;
};

// Try to confirm that \curl (\grad psi) == 0
// psi = a * sinh(r) * cos(theta)
// B = \grad psi
// \curl B = \curl (\grad psi) = 0
TEST_P(CurlParamTests, irrotational_field_has_zero_curl)
{
    using Mapping2D = CircularToCartesian<Rho, Theta, R, Z>;
    using ToroidalMapping = ToroidalToCylindrical<Mapping2D, Zeta, Phi>;
    using CylindricalMapping = CylindricalToCartesian<R, Z, Zeta, X, Y>;

    Coord<R, Z> origin_point(m_R0, 0.0);
    Mapping2D polar_to_RZ(origin_point);
    ToroidalMapping toroidal_to_cylindrical(polar_to_RZ);
    CylindricalMapping cylindrical_to_cartesian;
    CombinedMapping<CylindricalMapping, ToroidalMapping, Coord<Rho, Theta, Phi>>
            mapping(cylindrical_to_cartesian, toroidal_to_cylindrical);
    Curl curl(mapping);

    using BasisSpatial = VectorIndexSet<Rho, Theta, Phi>;
    using CovBasisSpatial = get_covariant_dims_t<BasisSpatial>;

    // Consider a point with rho > 0
    Coord<Rho, Theta, Phi> coord(m_rho, m_theta, Kokkos::numbers::pi);
    DTensor<CovBasisSpatial, CovBasisSpatial> grad_B;

    double const r = ddc::get<Rho>(coord);
    double const theta = ddc::get<Theta>(coord);
    // dB_r/dr = R0 * std::cos(theta) * std::sinh(r)
    // dB_r/dtheta = -R0 * std::sin(theta) * std::cosh(r)
    // dB_theta/dr = -R0 * std::sin(theta) * std::cosh(r)
    // dB_theta/dtheta = -R0 * std::cos(theta) * std::sinh(r)
    ddcHelper::get<Rho_cov, Rho_cov>(grad_B) = m_R0 * std::cos(theta) * std::sinh(r);
    ddcHelper::get<Rho_cov, Theta_cov>(grad_B) = -m_R0 * std::sin(theta) * std::cosh(r);
    ddcHelper::get<Rho_cov, Phi_cov>(grad_B) = 0.0;
    ddcHelper::get<Theta_cov, Rho_cov>(grad_B) = -m_R0 * std::sin(theta) * std::cosh(r);
    ddcHelper::get<Theta_cov, Theta_cov>(grad_B) = -m_R0 * std::cos(theta) * std::sinh(r);
    ddcHelper::get<Theta_cov, Phi_cov>(grad_B) = 0.0;
    ddcHelper::get<Phi_cov, Rho_cov>(grad_B) = 0.0;
    ddcHelper::get<Phi_cov, Theta_cov>(grad_B) = 0.0;
    ddcHelper::get<Phi_cov, Phi_cov>(grad_B) = 0.0;

    DVector<Rho, Theta, Phi> curl_B = curl(grad_B, coord);

    double eps = 1e-12;
    EXPECT_NEAR(ddcHelper::get<Rho>(curl_B), 0.0, eps);
    EXPECT_NEAR(ddcHelper::get<Theta>(curl_B), 0.0, eps);
    EXPECT_NEAR(ddcHelper::get<Phi>(curl_B), 0.0, eps);
}

// Test the curl for a simple B-field where the
// result can be analytically computed.
// Let's use a contravariant vector
// B = (0, 0, 1/(R0+r))
TEST_P(CurlParamTests, specific_field_and_point)
{
    using Mapping2D = CircularToCartesian<Rho, Theta, R, Z>;
    using ToroidalMapping = ToroidalToCylindrical<Mapping2D, Zeta, Phi>;
    using CylindricalMapping = CylindricalToCartesian<R, Z, Zeta, X, Y>;

    Coord<R, Z> origin_point(m_R0, 0.0);
    Mapping2D polar_to_RZ(origin_point);
    ToroidalMapping toroidal_to_cylindrical(polar_to_RZ);
    CylindricalMapping cylindrical_to_cartesian;
    CombinedMapping<CylindricalMapping, ToroidalMapping, Coord<Rho, Theta, Phi>>
            mapping(cylindrical_to_cartesian, toroidal_to_cylindrical);
    Curl curl(mapping);

    using BasisSpatial = VectorIndexSet<Rho, Theta, Phi>;
    using CovBasisSpatial = get_covariant_dims_t<BasisSpatial>;

    // Consider a point with rho > 0
    Coord<Rho, Theta, Phi> coord(m_rho, m_theta, Kokkos::numbers::pi);
    DTensor<CovBasisSpatial, CovBasisSpatial> grad_B;

    double const r = ddc::get<Rho>(coord);
    double const theta = ddc::get<Theta>(coord);
    double const cos_t = std::cos(theta);
    double const sin_t = std::sin(theta);

    // Compute dB_phi/dr, dB_phi/dtheta
    // B_i needs to be computed by g_ij B^j
    // dB_phi/dr = (R0 + r * cos_t) * (-R0 - r * cos_t + 2 * (R0 + r) * cos_t) / ((R0 + r) * (R0 + r))
    // dB_phi/dtheta = -2.0 * r * (R0 + r * cos_t) * sin_t / (R0 + r)
    ddcHelper::get<Rho_cov, Rho_cov>(grad_B) = 0.0;
    ddcHelper::get<Rho_cov, Theta_cov>(grad_B) = 0.0;
    ddcHelper::get<Rho_cov, Phi_cov>(grad_B) = 0.0;
    ddcHelper::get<Theta_cov, Rho_cov>(grad_B) = 0.0;
    ddcHelper::get<Theta_cov, Theta_cov>(grad_B) = 0.0;
    ddcHelper::get<Theta_cov, Phi_cov>(grad_B) = 0.0;
    ddcHelper::get<Phi_cov, Rho_cov>(grad_B) = (m_R0 + r * cos_t)
                                               * (-m_R0 - r * cos_t + 2 * (m_R0 + r) * cos_t)
                                               / ((m_R0 + r) * (m_R0 + r));
    ddcHelper::get<Phi_cov, Theta_cov>(grad_B) = -2.0 * r * (m_R0 + r * cos_t) * sin_t / (m_R0 + r);
    ddcHelper::get<Phi_cov, Phi_cov>(grad_B) = 0.0;

    DVector<Rho, Theta, Phi> curl_B = curl(grad_B, coord);

    double eps = 1e-12;
    double ref_curl_B_r = -2 * sin_t / (r + m_R0);
    double ref_curl_B_theta = (m_R0 - (2 * m_R0 + r) * cos_t) / (r * (r + m_R0) * (r + m_R0));
    EXPECT_NEAR(ddcHelper::get<Rho>(curl_B), ref_curl_B_r, eps);
    EXPECT_NEAR(ddcHelper::get<Theta>(curl_B), ref_curl_B_theta, eps);
    EXPECT_NEAR(ddcHelper::get<Phi>(curl_B), 0.0, eps);
}

// Parameterisation over rho and theta
INSTANTIATE_TEST_SUITE_P(
        CurlTests,
        CurlParamTests,
        ::testing::Values(
                std::tuple<double, double>(0.1, Kokkos::numbers::pi * 0.5),
                std::tuple<double, double>(0.2, Kokkos::numbers::pi),
                std::tuple<double, double>(0.3, Kokkos::numbers::pi * 1.5)));
