#include <gtest/gtest.h>

#include "../mapping/geometry_mapping_tests.hpp"

#include "circular_to_cartesian.hpp"
#include "combined_mapping.hpp"
#include "cylindrical_to_cartesian.hpp"
#include "gyrokinetic_poisson_bracket.hpp"
#include "mesh_builder.hpp"
#include "toroidal_to_cylindrical.hpp"


TEST(GyrokineticPoissonBracket, Anticommutativity)
{
    using Mapping2D = CircularToCartesian<Rho, Theta, R, Z>;
    using ToroidalMapping = ToroidalToCylindrical<Mapping2D, Zeta, Phi>;
    using CylindricalMapping = CylindricalToCartesian<R, Z, Zeta, X, Y>;
    double major_radius(6.2);
    Mapping2D polar_to_RZ(major_radius);
    ToroidalMapping toroidal_to_cylindrical(polar_to_RZ);
    CylindricalMapping cylindrical_to_cartesian;
    CombinedMapping<CylindricalMapping, ToroidalMapping, Coord<Rho, Theta, Phi>>
            mapping(cylindrical_to_cartesian, toroidal_to_cylindrical);
    GyrokineticPoissonBracket calculate_poisson_bracket(mapping);

    using BasisSpatial = VectorIndexSet<Rho, Theta, Phi>;
    using CovBasisSpatial = get_covariant_dims_t<BasisSpatial>;

    IdxStepRho nrho(10);
    IdxStepTheta ntheta(10);
    IdxStepPhi nphi(10);
    ddc::init_discrete_space<GridRho>(
            GridRho::init<GridRho>(Coord<Rho>(0.0), Coord<Rho>(1.0), nrho));
    ddc::init_discrete_space<GridTheta>(GridTheta::init<GridTheta>(
            build_uniform_break_points(Coord<Theta>(0.0), Coord<Theta>(2 * M_PI), ntheta)));
    ddc::init_discrete_space<GridPhi>(
            GridPhi::init<GridPhi>(Coord<Phi>(0.0), Coord<Phi>(2 * M_PI), nphi));

    IdxRangeRho idx_range_rho(
            IdxRho(1), // Exclude central point where Jacobian/gradient is poorly defined
            nrho);
    IdxRangeTheta idx_range_theta(IdxTheta(0), ntheta);
    IdxRangePhi idx_range_phi(IdxPhi(0), nphi);
    IdxRangeRhoThetaPhi idx_range(idx_range_rho, idx_range_theta, idx_range_phi);

    DFieldMem<IdxRangeRhoThetaPhi> poisson_bracket_alloc(idx_range);
    DFieldMem<IdxRangeRhoThetaPhi> reverse_poisson_bracket_alloc(idx_range);
    DVectorFieldMem<IdxRangeRhoThetaPhi, CovBasisSpatial> df_alloc(idx_range);
    DVectorFieldMem<IdxRangeRhoThetaPhi, CovBasisSpatial> dg_alloc(idx_range);
    DVectorFieldMem<IdxRangeRhoThetaPhi, BasisSpatial> B_alloc(idx_range);

    DField<IdxRangeRhoThetaPhi> poisson_bracket = get_field(poisson_bracket_alloc);
    DField<IdxRangeRhoThetaPhi> reverse_poisson_bracket = get_field(reverse_poisson_bracket_alloc);
    DVectorField<IdxRangeRhoThetaPhi, CovBasisSpatial> df = get_field(df_alloc);
    DVectorField<IdxRangeRhoThetaPhi, CovBasisSpatial> dg = get_field(dg_alloc);
    DVectorField<IdxRangeRhoThetaPhi, BasisSpatial> B = get_field(B_alloc);

    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            idx_range,
            KOKKOS_LAMBDA(IdxRhoThetaPhi idx) {
                // f(rho, theta, phi) = rho
                ddcHelper::get<Rho_cov>(df)(idx) = 1;
                ddcHelper::get<Theta_cov>(df)(idx) = 0;
                ddcHelper::get<Phi_cov>(df)(idx) = 0;
                // g(rho, theta, phi) = phi
                ddcHelper::get<Rho_cov>(dg)(idx) = 0;
                ddcHelper::get<Theta_cov>(dg)(idx) = 0;
                ddcHelper::get<Phi_cov>(dg)(idx) = 1;
                // B(rho, theta, phi) = 0.1 \hat{theta} + 0.9 \hat{phi}
                ddcHelper::get<Rho>(B)(idx) = 0;
                ddcHelper::get<Theta>(B)(idx) = 0.1;
                ddcHelper::get<Phi>(B)(idx) = 0.9;
            });

    calculate_poisson_bracket(
            Kokkos::DefaultExecutionSpace(),
            poisson_bracket,
            get_const_field(df),
            get_const_field(dg),
            get_const_field(B));
    calculate_poisson_bracket(
            Kokkos::DefaultExecutionSpace(),
            reverse_poisson_bracket,
            get_const_field(dg),
            get_const_field(df),
            get_const_field(B));

    auto poisson_bracket_host = ddc::create_mirror_view_and_copy(poisson_bracket);
    auto reverse_poisson_bracket_host = ddc::create_mirror_view_and_copy(reverse_poisson_bracket);
    ddc::for_each(idx_range, [&](IdxRhoThetaPhi idx) {
        EXPECT_NEAR(poisson_bracket_host(idx), -reverse_poisson_bracket_host(idx), 1e-13);
    });
}

TEST(GyrokineticPoissonBracket, Bilinearity) {}

TEST(GyrokineticPoissonBracket, Leibniz) {}

TEST(GyrokineticPoissonBracket, Jacobi) {}
