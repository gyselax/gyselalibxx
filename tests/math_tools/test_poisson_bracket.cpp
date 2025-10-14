#include <cmath>
#include <vector>

#include <ddc/ddc.hpp>

#include <gtest/gtest.h>

#include <Kokkos_Core.hpp>

#include "../coord_transformations/geometry_coord_transformations_tests.hpp"

#include "circular_to_cartesian.hpp"
#include "combined_mapping.hpp"
#include "cylindrical_to_cartesian.hpp"
#include "ddc_alias_inline_functions.hpp"
#include "ddc_aliases.hpp"
#include "indexed_tensor.hpp"
#include "lie_poisson_bracket.hpp"
#include "math_tools.hpp"
#include "mesh_builder.hpp"
#include "metric_tensor_evaluator.hpp"
#include "static_tensors.hpp"
#include "tensor.hpp"
#include "tensor_common.hpp"
#include "toroidal_to_cylindrical.hpp"
#include "vector_field.hpp"
#include "vector_field_common.hpp"
#include "vector_field_mem.hpp"
#include "vector_index_tools.hpp"



/**
 * @brief Test the Lie Poisson Bracket vs analytical result in the case of an axisymmetric field,
 * i.e. d_\varhpi F = d_\varhpi G = 0. ASssuming on top of this that b_r = 0, which is always true
 * if nabla r is perpenduclar to the magnetic field (or said differently, r is a label of the flux surface)
 * Then in this case [F,G] =  b_{\varphi}/J_{x} \times (\partial_{r} F partial_{\theta} G)
 *                                                     -partial_{\theta} F (\partial_{r} G )
 */

void compute_and_test_Lie_Poisson_Bracket()
{
    using Mapping2D = CircularToCartesian<Rho, Theta, R, Z>;
    using ToroidalMapping = ToroidalToCylindrical<Mapping2D, Zeta, Phi>;
    using CylindricalMapping = CylindricalToCartesian<R, Z, Zeta, X, Y>;
    Coord<R, Z> origin_point(6.2, 0.0);
    Mapping2D polar_to_RZ(origin_point);
    ToroidalMapping toroidal_to_cylindrical(polar_to_RZ);
    CylindricalMapping cylindrical_to_cartesian;
    CombinedMapping<CylindricalMapping, ToroidalMapping, Coord<Rho, Theta, Phi>>
            mapping(cylindrical_to_cartesian, toroidal_to_cylindrical);
    LiePoissonBracket calculate_poisson_bracket(mapping);
    MetricTensorEvaluator get_metric_tensor(mapping);

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

    IdxRangeRho idx_range_rho(
            IdxRho(1), // Exclude central point where Jacobian/gradient is poorly defined
            nrho);
    IdxRangeTheta idx_range_theta(IdxTheta(0), ntheta);
    IdxRangePhi idx_range_phi(IdxPhi(0), nphi);
    IdxRangeRhoThetaPhi idx_range(idx_range_rho, idx_range_theta, idx_range_phi);

    DFieldMem<IdxRangeRhoThetaPhi> poisson_bracket_alloc(idx_range);
    DFieldMem<IdxRangeRhoThetaPhi> anal_alloc(idx_range);
    DVectorFieldMem<IdxRangeRhoThetaPhi, CovBasisSpatial> df_alloc(idx_range);
    DVectorFieldMem<IdxRangeRhoThetaPhi, CovBasisSpatial> dg_alloc(idx_range);
    DVectorFieldMem<IdxRangeRhoThetaPhi, BasisSpatial> B_alloc(idx_range);

    DField<IdxRangeRhoThetaPhi> poisson_bracket_a = get_field(poisson_bracket_alloc);
    DField<IdxRangeRhoThetaPhi> poisson_bracket_b = get_field(poisson_bracket_alloc);
    DField<IdxRangeRhoThetaPhi> analytical_matrix = get_field(anal_alloc);
    DVectorField<IdxRangeRhoThetaPhi, CovBasisSpatial> df_a = get_field(df_alloc);
    DVectorField<IdxRangeRhoThetaPhi, CovBasisSpatial> dg_a = get_field(dg_alloc);
    DVectorField<IdxRangeRhoThetaPhi, CovBasisSpatial> df_b = get_field(df_alloc);
    DVectorField<IdxRangeRhoThetaPhi, CovBasisSpatial> dg_b = get_field(dg_alloc);
    DVectorField<IdxRangeRhoThetaPhi, BasisSpatial> B = get_field(B_alloc);

    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            idx_range,
            KOKKOS_LAMBDA(IdxRhoThetaPhi idx) {
                const double rho = ddc::coordinate(ddc::select<GridRho>(idx));
                const double theta = ddc::coordinate(ddc::select<GridTheta>(idx));
                Coord<Rho, Theta, Phi> coord = ddc ::coordinate(idx);
                // f_a(rho, theta, phi) = 0.5 * rho^2
                ddcHelper::get<Rho_cov>(df_a)(idx) = rho;
                ddcHelper::get<Theta_cov>(df_a)(idx) = 0;
                ddcHelper::get<Phi_cov>(df_a)(idx) = 0;
                // g_a(rho, theta, phi) = 0.5 * theta^2
                ddcHelper::get<Rho_cov>(dg_a)(idx) = 0;
                ddcHelper::get<Theta_cov>(dg_a)(idx) = theta;
                ddcHelper::get<Phi_cov>(dg_a)(idx) = 0;
                // f_b(rho, theta, phi) = -0.5 * theta^2
                ddcHelper::get<Rho_cov>(df_b)(idx) = 0;
                ddcHelper::get<Theta_cov>(df_b)(idx) = -theta;
                ddcHelper::get<Phi_cov>(df_b)(idx) = 0;
                // g_b(rho, theta, phi) = 0.5 * rho^2
                ddcHelper::get<Rho_cov>(dg_b)(idx) = rho;
                ddcHelper::get<Theta_cov>(dg_b)(idx) = 0;
                ddcHelper::get<Phi_cov>(dg_b)(idx) = 0;
                // B(rho, theta, phi) = 0.1 \hat{theta} + 0.9 \hat{phi}
                const double B_theta = 0.1;
                const double B_phi = 0.9;
                DVector<Rho_cov, Theta_cov, Phi_cov> B_cov(0.0, B_theta, B_phi);
                DTensor<BasisSpatial, BasisSpatial> metric_tensor
                        = get_metric_tensor.inverse(coord);
                DVector<Rho, Theta, Phi> B_cont
                        = tensor_mul(index<'i', 'j'>(metric_tensor), index<'j'>(B_cov));
                ddcHelper::assign_vector_field_element(B, idx, B_cont);
                const double Jx = mapping.jacobian(coord);
                const double norm_B
                        = Kokkos::sqrt(tensor_mul(index<'i'>(B_cov), index<'i'>(B(idx))));
                analytical_matrix(idx) = (B_phi / norm_B) * rho * theta / Jx;
            });

    calculate_poisson_bracket(
            Kokkos::DefaultExecutionSpace(),
            poisson_bracket_a,
            get_const_field(df_a),
            get_const_field(dg_a),
            get_const_field(B));
    calculate_poisson_bracket(
            Kokkos::DefaultExecutionSpace(),
            poisson_bracket_b,
            get_const_field(df_b),
            get_const_field(dg_b),
            get_const_field(B));

    auto poisson_bracket_a_host = ddc::create_mirror_view_and_copy(poisson_bracket_a);
    auto poisson_bracket_b_host = ddc::create_mirror_view_and_copy(poisson_bracket_b);
    auto analytical_matrix_host = ddc::create_mirror_view_and_copy(analytical_matrix);

    ddc::for_each(idx_range, [&](IdxRhoThetaPhi idx) {
        EXPECT_NEAR(poisson_bracket_a_host(idx), poisson_bracket_b_host(idx), 1e-13);
    });
    ddc::for_each(idx_range, [&](IdxRhoThetaPhi idx) {
        EXPECT_NEAR(poisson_bracket_a_host(idx), analytical_matrix_host(idx), 1e-13);
    });
}

TEST(LiePoissonBracket, axisymmetric_tokamak)
{
    compute_and_test_Lie_Poisson_Bracket();
}

TEST(LiePoissonBracket, batched_elementwise)
{
    using Mapping2D = CircularToCartesian<Rho, Theta, R, Z>;
    using ToroidalMapping = ToroidalToCylindrical<Mapping2D, Zeta, Phi>;
    using CylindricalMapping = CylindricalToCartesian<R, Z, Zeta, X, Y>;
    Coord<R, Z> origin_point(6.2, 0.0);
    Mapping2D polar_to_RZ(origin_point);
    ToroidalMapping toroidal_to_cylindrical(polar_to_RZ);
    CylindricalMapping cylindrical_to_cartesian;
    CombinedMapping<CylindricalMapping, ToroidalMapping, Coord<Rho, Theta, Phi>>
            mapping(cylindrical_to_cartesian, toroidal_to_cylindrical);
    LiePoissonBracket calculate_poisson_bracket(mapping);
    MetricTensorEvaluator get_metric_tensor(mapping);

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
            = ddc::coordinate(Idx<GridRho, GridTheta, GridPhi>(nrho - 1, 0, 0));
    DVector<Rho_cov, Theta_cov, Phi_cov> df(0.0, 0.0, 1.0);
    DVector<Rho_cov, Theta_cov, Phi_cov> B_cov(1.0, 1.0, 1.0);
    DTensor<BasisSpatial, BasisSpatial> inv_G = get_metric_tensor.inverse(coord);
    double B_norm = norm(inv_G, B_cov);
    B_cov /= B_norm;
    DVector<Rho, Theta, Phi> B = tensor_mul(index<'i', 'j'>(inv_G), index<'j'>(B_cov));
    IdentityTensor<double, BasisSpatial, CovBasisSpatial> id;

    DVector<Rho, Theta, Phi> lp = calculate_poisson_bracket(df, id, B, coord);

    double J = mapping.jacobian(coord);

    double expected = 1.0 / (J * B_norm);
    EXPECT_DOUBLE_EQ(ddcHelper::get<Rho>(lp), expected);
    EXPECT_DOUBLE_EQ(ddcHelper::get<Theta>(lp), -expected);
    EXPECT_DOUBLE_EQ(ddcHelper::get<Phi>(lp), 0.0);
}
