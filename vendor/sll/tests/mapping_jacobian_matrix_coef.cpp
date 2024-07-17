#include <array>
#include <cassert>

#include <ddc/ddc.hpp>

#include "sll/mapping/analytical_invertible_curvilinear2d_to_cartesian.hpp"
#include "sll/mapping/circular_to_cartesian.hpp"
#include "sll/mapping/curvilinear2d_to_cartesian.hpp"
#include "sll/mapping/czarny_to_cartesian.hpp"
#include "sll/mapping/discrete_mapping_to_cartesian.hpp"

#include "test_utils.hpp"



namespace {
struct DimX
{
};
struct DimY
{
};
struct DimR
{
    static bool constexpr PERIODIC = false;
};

struct DimP
{
    static bool constexpr PERIODIC = true;
};

using CoordR = ddc::Coordinate<DimR>;
using CoordP = ddc::Coordinate<DimP>;
using CoordRP = ddc::Coordinate<DimR, DimP>;

int constexpr BSDegree = 3;

struct BSplinesR : ddc::NonUniformBSplines<DimR, BSDegree>
{
};
struct BSplinesP : ddc::NonUniformBSplines<DimP, BSDegree>
{
};

using InterpPointsR = ddc::
        GrevilleInterpolationPoints<BSplinesR, ddc::BoundCond::GREVILLE, ddc::BoundCond::GREVILLE>;
using InterpPointsP = ddc::
        GrevilleInterpolationPoints<BSplinesP, ddc::BoundCond::PERIODIC, ddc::BoundCond::PERIODIC>;

struct IDimR : InterpPointsR::interpolation_discrete_dimension_type
{
};
struct IDimP : InterpPointsP::interpolation_discrete_dimension_type
{
};

using SplineRPBuilder = ddc::SplineBuilder2D<
        Kokkos::DefaultHostExecutionSpace,
        Kokkos::DefaultHostExecutionSpace::memory_space,
        BSplinesR,
        BSplinesP,
        IDimR,
        IDimP,
        ddc::BoundCond::GREVILLE,
        ddc::BoundCond::GREVILLE,
        ddc::BoundCond::PERIODIC,
        ddc::BoundCond::PERIODIC,
        ddc::SplineSolver::GINKGO,
        IDimR,
        IDimP>;

using SplineRPEvaluator = ddc::SplineEvaluator2D<
        Kokkos::DefaultHostExecutionSpace,
        Kokkos::DefaultHostExecutionSpace::memory_space,
        BSplinesR,
        BSplinesP,
        IDimR,
        IDimP,
        ddc::NullExtrapolationRule,
        ddc::NullExtrapolationRule,
        ddc::PeriodicExtrapolationRule<DimP>,
        ddc::PeriodicExtrapolationRule<DimP>,
        IDimR,
        IDimP>;

using BSDomainR = ddc::DiscreteDomain<BSplinesR>;
using BSDomainP = ddc::DiscreteDomain<BSplinesP>;
using BSDomainRP = ddc::DiscreteDomain<BSplinesR, BSplinesP>;

using IDomainR = ddc::DiscreteDomain<IDimR>;
using IDomainP = ddc::DiscreteDomain<IDimP>;
using IDomainRP = ddc::DiscreteDomain<IDimR, IDimP>;

using IndexR = ddc::DiscreteElement<IDimR>;
using IndexP = ddc::DiscreteElement<IDimP>;
using IndexRP = ddc::DiscreteElement<IDimR, IDimP>;

using IVectR = ddc::DiscreteVector<IDimR>;
using IVectP = ddc::DiscreteVector<IDimP>;
using IVectRP = ddc::DiscreteVector<IDimR, IDimP>;

using IDomainRP = ddc::DiscreteDomain<IDimR, IDimP>;


template <class ElementType>
using FieldRP = ddc::Chunk<ElementType, IDomainRP>;

using Matrix_2x2 = std::array<std::array<double, 2>, 2>;

/**
 * @brief Check if the two matrix given as input are the same.
 *
 * The error tolerance is given at 1e-15.
 *
 * @param[in] matrix
 * 			The matrix defined with the matrix function.
 * @param[in] matrix_coeff
 * 			The matrix defined with the matrix coefficient functions.
 */
void check_matrix(Matrix_2x2 matrix, Matrix_2x2 matrix_coeff)
{
    const double TOL = 1e-15;
    constexpr std::size_t N = 2;

    for (std::size_t i(0); i < N; ++i) {
        for (std::size_t j(0); j < N; ++j) {
            const double id_val = fabs(matrix[i][j] - matrix_coeff[i][j]);
            EXPECT_NEAR(id_val, 0., TOL);
        }
    }
}

} // namespace


/**
 * @brief A class for the Google tests.
 */
class JacobianMatrixAndJacobianCoefficients
    : public testing::TestWithParam<std::tuple<std::size_t, std::size_t>>
{
};


// Circular mapping ------------------------------------------------
TEST_P(JacobianMatrixAndJacobianCoefficients, MatrixCircMap)
{
    auto const [Nr, Nt] = GetParam();
    const CircularToCartesian<DimX, DimY, DimR, DimP> mapping;

    CoordR const r_min(0.0);
    CoordR const r_max(1.0);
    IVectR const r_size(Nr);

    CoordP const p_min(0.0);
    CoordP const p_max(2.0 * M_PI);
    IVectP const p_size(Nt);

    IndexR const r_start(1); // avoid singular point at r = 0.
    IndexP const p_start(0);

    double const dr((r_max - r_min) / r_size);
    double const dp((p_max - p_min) / p_size);

    ddc::DiscreteDomain<IDimR> domain_r(r_start, r_size);
    ddc::DiscreteDomain<IDimP> domain_p(p_start, p_size);
    ddc::DiscreteDomain<IDimR, IDimP> grid(domain_r, domain_p);

    FieldRP<CoordRP> coords(grid);
    ddc::for_each(grid, [&](IndexRP const irp) {
        coords(irp) = CoordRP(
                r_min + dr * ddc::select<IDimR>(irp).uid(),
                p_min + dp * ddc::select<IDimR>(irp).uid());
    });

    // Test for each coordinates if the coefficients defined by the coefficients functions
    //are the same as the coefficients in the matrix function.
    // --- for the Jacobian matrix:
    ddc::for_each(grid, [&](IndexRP const irp) {
        Matrix_2x2 Jacobian_matrix;
        Matrix_2x2 Jacobian_matrix_coeff;

        mapping.jacobian_matrix(coords(irp), Jacobian_matrix);
        Jacobian_matrix_coeff[0][0] = mapping.jacobian_11(coords(irp));
        Jacobian_matrix_coeff[0][1] = mapping.jacobian_12(coords(irp));
        Jacobian_matrix_coeff[1][0] = mapping.jacobian_21(coords(irp));
        Jacobian_matrix_coeff[1][1] = mapping.jacobian_22(coords(irp));

        check_matrix(Jacobian_matrix, Jacobian_matrix_coeff);
    });

    // --- for the inverse Jacobian matrix:
    ddc::for_each(grid, [&](IndexRP const irp) {
        Matrix_2x2 inv_Jacobian_matrix;
        Matrix_2x2 inv_Jacobian_matrix_coeff;

        mapping.inv_jacobian_matrix(coords(irp), inv_Jacobian_matrix);
        inv_Jacobian_matrix_coeff[0][0] = mapping.inv_jacobian_11(coords(irp));
        inv_Jacobian_matrix_coeff[0][1] = mapping.inv_jacobian_12(coords(irp));
        inv_Jacobian_matrix_coeff[1][0] = mapping.inv_jacobian_21(coords(irp));
        inv_Jacobian_matrix_coeff[1][1] = mapping.inv_jacobian_22(coords(irp));

        check_matrix(inv_Jacobian_matrix, inv_Jacobian_matrix_coeff);
    });
}



// Czarny mapping --------------------------------------------------
TEST_P(JacobianMatrixAndJacobianCoefficients, MatrixCzarMap)
{
    auto const [Nr, Nt] = GetParam();
    const CzarnyToCartesian<DimX, DimY, DimR, DimP> mapping(0.3, 1.4);

    CoordR const r_min(0.0);
    CoordR const r_max(1.0);
    IVectR const r_size(Nr);

    CoordP const p_min(0.0);
    CoordP const p_max(2.0 * M_PI);
    IVectP const p_size(Nt);

    IndexR const r_start(1); // avoid singular point at r = 0.
    IndexP const p_start(0);

    double const dr((r_max - r_min) / r_size);
    double const dp((p_max - p_min) / p_size);

    ddc::DiscreteDomain<IDimR> domain_r(r_start, r_size);
    ddc::DiscreteDomain<IDimP> domain_p(p_start, p_size);
    ddc::DiscreteDomain<IDimR, IDimP> grid(domain_r, domain_p);

    FieldRP<CoordRP> coords(grid);
    ddc::for_each(grid, [&](IndexRP const irp) {
        coords(irp) = CoordRP(
                r_min + dr * ddc::select<IDimR>(irp).uid(),
                p_min + dp * ddc::select<IDimR>(irp).uid());
    });

    // Test for each coordinates if the coefficients defined by the coefficients functions
    //are the same as the coefficients in the matrix function.
    // --- for the Jacobian matrix:
    ddc::for_each(grid, [&](IndexRP const irp) {
        Matrix_2x2 Jacobian_matrix;
        Matrix_2x2 Jacobian_matrix_coeff;

        mapping.jacobian_matrix(coords(irp), Jacobian_matrix);
        Jacobian_matrix_coeff[0][0] = mapping.jacobian_11(coords(irp));
        Jacobian_matrix_coeff[0][1] = mapping.jacobian_12(coords(irp));
        Jacobian_matrix_coeff[1][0] = mapping.jacobian_21(coords(irp));
        Jacobian_matrix_coeff[1][1] = mapping.jacobian_22(coords(irp));

        check_matrix(Jacobian_matrix, Jacobian_matrix_coeff);
    });

    // --- for the inverseJacobian matrix:
    ddc::for_each(grid, [&](IndexRP const irp) {
        Matrix_2x2 inv_Jacobian_matrix;
        Matrix_2x2 inv_Jacobian_matrix_coeff;

        mapping.inv_jacobian_matrix(coords(irp), inv_Jacobian_matrix);
        inv_Jacobian_matrix_coeff[0][0] = mapping.inv_jacobian_11(coords(irp));
        inv_Jacobian_matrix_coeff[0][1] = mapping.inv_jacobian_12(coords(irp));
        inv_Jacobian_matrix_coeff[1][0] = mapping.inv_jacobian_21(coords(irp));
        inv_Jacobian_matrix_coeff[1][1] = mapping.inv_jacobian_22(coords(irp));

        check_matrix(inv_Jacobian_matrix, inv_Jacobian_matrix_coeff);
    });
}



// Discrete Czarny mapping -----------------------------------------
TEST_P(JacobianMatrixAndJacobianCoefficients, MatrixDiscCzarMap)
{
    auto const [Nr, Nt] = GetParam();
    const CzarnyToCartesian<DimX, DimY, DimR, DimP> analytical_mapping(0.3, 1.4);

    CoordR const r_min(0.0);
    CoordR const r_max(1.0);
    IVectR const r_size(Nr);

    CoordP const p_min(0.0);
    CoordP const p_max(2.0 * M_PI);
    IVectP const p_size(Nt);

    double const dr((r_max - r_min) / r_size);
    double const dp((p_max - p_min) / p_size);

    std::vector<CoordR> r_knots(r_size + 1);
    std::vector<CoordP> p_knots(p_size + 1);

    for (int i(0); i < r_size + 1; ++i) {
        r_knots[i] = CoordR(r_min + i * dr);
    }
    r_knots[r_size] = CoordR(r_max);
    for (int i(0); i < p_size + 1; ++i) {
        p_knots[i] = CoordP(p_min + i * dp);
    }

    ddc::init_discrete_space<BSplinesR>(r_knots);
    ddc::init_discrete_space<BSplinesP>(p_knots);

    ddc::init_discrete_space<IDimR>(InterpPointsR::get_sampling<IDimR>());
    ddc::init_discrete_space<IDimP>(InterpPointsP::get_sampling<IDimP>());

    IDomainR interpolation_domain_R(InterpPointsR::get_domain<IDimR>());
    IDomainP interpolation_domain_P(InterpPointsP::get_domain<IDimP>());
    IDomainRP grid(interpolation_domain_R, interpolation_domain_P);

    SplineRPBuilder builder(grid);
    ddc::NullExtrapolationRule r_extrapolation_rule;
    ddc::PeriodicExtrapolationRule<DimP> p_extrapolation_rule;
    SplineRPEvaluator evaluator(
            r_extrapolation_rule,
            r_extrapolation_rule,
            p_extrapolation_rule,
            p_extrapolation_rule);
    DiscreteToCartesian mapping
            = DiscreteToCartesian<DimX, DimY, SplineRPBuilder, SplineRPEvaluator>::
                    analytical_to_discrete(analytical_mapping, builder, evaluator);

    // Test for each coordinates if the coefficients defined by the coefficients functions
    //are the same as the coefficients in the matrix function.
    ddc::for_each(grid, [&](IndexRP const irp) {
        const CoordRP coord_rp(ddc::coordinate(irp));
        const double r = ddc::get<DimR>(coord_rp);
        if (fabs(r) > 1e-15) {
            // --- for the Jacobian matrix:
            Matrix_2x2 Jacobian_matrix;
            Matrix_2x2 Jacobian_matrix_coeff;

            mapping.jacobian_matrix(coord_rp, Jacobian_matrix);
            Jacobian_matrix_coeff[0][0] = mapping.jacobian_11(coord_rp);
            Jacobian_matrix_coeff[0][1] = mapping.jacobian_12(coord_rp);
            Jacobian_matrix_coeff[1][0] = mapping.jacobian_21(coord_rp);
            Jacobian_matrix_coeff[1][1] = mapping.jacobian_22(coord_rp);

            check_matrix(Jacobian_matrix, Jacobian_matrix_coeff);


            // --- for the inverse Jacobian matrix:
            Matrix_2x2 inv_Jacobian_matrix;
            Matrix_2x2 inv_Jacobian_matrix_coeff;

            mapping.inv_jacobian_matrix(coord_rp, inv_Jacobian_matrix);
            inv_Jacobian_matrix_coeff[0][0] = mapping.inv_jacobian_11(coord_rp);
            inv_Jacobian_matrix_coeff[0][1] = mapping.inv_jacobian_12(coord_rp);
            inv_Jacobian_matrix_coeff[1][0] = mapping.inv_jacobian_21(coord_rp);
            inv_Jacobian_matrix_coeff[1][1] = mapping.inv_jacobian_22(coord_rp);
        }
    });
}



INSTANTIATE_TEST_SUITE_P(
        MyGroup,
        JacobianMatrixAndJacobianCoefficients,
        testing::Combine(testing::Values<std::size_t>(40), testing::Values<std::size_t>(80)));
