#include <array>
#include <cassert>
#include <cstdlib>

#include <ddc/ddc.hpp>
#include <ddc/kernels/splines.hpp>

#include "sll/mapping/analytical_invertible_curvilinear2d_to_cartesian.hpp"
#include "sll/mapping/circular_to_cartesian.hpp"
#include "sll/mapping/curvilinear2d_to_cartesian.hpp"
#include "sll/mapping/czarny_to_cartesian.hpp"
#include "sll/mapping/refined_discrete_mapping_to_cartesian.hpp"

#include "test_utils.hpp"



namespace {
struct RDimX
{
};
struct RDimY
{
};
struct RDimR
{
    static bool constexpr PERIODIC = false;
};

struct RDimP
{
    static bool constexpr PERIODIC = true;
};

using CoordR = ddc::Coordinate<RDimR>;
using CoordP = ddc::Coordinate<RDimP>;
using CoordRP = ddc::Coordinate<RDimR, RDimP>;

using CoordXY = ddc::Coordinate<RDimX, RDimY>;

int constexpr BSDegree = 3;

using BSplinesR = ddc::NonUniformBSplines<RDimR, BSDegree>;
using BSplinesP = ddc::NonUniformBSplines<RDimP, BSDegree>;


using SplineInterpPointsR = ddc::
        GrevilleInterpolationPoints<BSplinesR, ddc::BoundCond::GREVILLE, ddc::BoundCond::GREVILLE>;
using SplineInterpPointsP = ddc::
        GrevilleInterpolationPoints<BSplinesP, ddc::BoundCond::PERIODIC, ddc::BoundCond::PERIODIC>;

using IDimR = typename SplineInterpPointsR::interpolation_mesh_type;
using IDimP = typename SplineInterpPointsP::interpolation_mesh_type;

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
        ddc::ConstantExtrapolationRule<RDimR, RDimP>,
        ddc::ConstantExtrapolationRule<RDimR, RDimP>,
        ddc::PeriodicExtrapolationRule<RDimP>,
        ddc::PeriodicExtrapolationRule<RDimP>,
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

template <class ElementType>
using FieldR = ddc::Chunk<ElementType, IDomainR>;

template <class ElementType>
using FieldP = ddc::Chunk<ElementType, IDomainP>;

template <class ElementType>
using FieldRP = ddc::Chunk<ElementType, IDomainRP>;



using IDomainRP = ddc::DiscreteDomain<IDimR, IDimP>;


using CzarnyMapping = CzarnyToCartesian<RDimX, RDimY, RDimR, RDimP>;
using CircularMapping = CircularToCartesian<RDimX, RDimY, RDimR, RDimP>;


template <class ElementType>
using FieldRP = ddc::Chunk<ElementType, IDomainRP>;

using Matrix_2x2 = std::array<std::array<double, 2>, 2>;

/**
 * @brief Compare the values given by the analytical mapping and
 * the other mapping on the grid points.
 *
 * The error tolerance is given at 5e-10.
 *
 * @param[in] mapping
 *          The mapping we are testing.
 * @param[in] analytical_mappping
 *          The mapping analytically defined.
 * @param[in] domain
 *          The domain on which we test the values.
 */
template <class Mapping, class DiscreteMapping>
double check_value_on_grid(
        DiscreteMapping const& mapping,
        Mapping const& analytical_mapping,
        IDomainRP const& domain)
{
    const double TOL = 1e-7;
    double max_err = 0.;
    ddc::for_each(domain, [&](IndexRP const irp) {
        const CoordRP coord(ddc::coordinate(irp));
        const CoordXY discrete_coord = mapping(coord);
        const CoordXY analytical_coord = analytical_mapping(coord);

        EXPECT_NEAR(ddc::get<RDimX>(discrete_coord), ddc::get<RDimX>(analytical_coord), TOL);
        EXPECT_NEAR(ddc::get<RDimY>(discrete_coord), ddc::get<RDimY>(analytical_coord), TOL);

        const double diff_x
                = double(ddc::get<RDimX>(discrete_coord) - ddc::get<RDimX>(analytical_coord));
        const double diff_y
                = double(ddc::get<RDimY>(discrete_coord) - ddc::get<RDimY>(analytical_coord));

        max_err = max_err > diff_x ? max_err : diff_x;
        max_err = max_err > diff_y ? max_err : diff_y;
    });

    return max_err;
}

/**
 * @brief Compare the values given by the analytical mapping and
 * the other mapping not on the grid points.
 *
 * The error tolerance is given at 5e-6.
 * The expected convergence order not on the grid points is
 * d + 1 where d is the degree of B-splines.
 *
 * @param[in] mapping
 *          The mapping we are testing.
 * @param[in] analytical_mappping
 *          The mapping analytically defined.
 * @param[in] domain
 *          The domain on which we test the values.
 */
template <class Mapping, class DiscreteMapping>
double check_value_not_on_grid(
        DiscreteMapping const& mapping,
        Mapping const& analytical_mapping,
        IDomainRP const& domain)
{
    std::srand(100);

    FieldRP<CoordRP> coords(domain);
    IndexR ir_max(ddc::select<IDimR>(domain).back());
    IndexP ip_max(ddc::select<IDimP>(domain).back());
    ddc::for_each(domain, [&](IndexRP const irp) {
        IndexR ir(ddc::select<IDimR>(irp));
        CoordR coordr_0 = ddc::coordinate(ir);
        CoordR coordr_1 = ddc::coordinate(ir + 1);
        double coord_r;
        if (ir.uid() < ir_max.uid()) {
            double factor = double(std::rand()) / RAND_MAX;
            coord_r = coordr_0 + (coordr_1 - coordr_0) * factor;
        } else {
            coord_r = coordr_0;
        }

        IndexP ip(ddc::select<IDimP>(irp));
        CoordP coordp_0 = ddc::coordinate(ip);
        CoordP coordp_1 = ddc::coordinate(ip + 1);
        double coord_p;
        if (ip.uid() < ip_max.uid()) {
            double factor = double(std::rand()) / RAND_MAX;
            coord_p = coordp_0 + (coordp_1 - coordp_0) * factor;
        } else {
            coord_p = coordp_0;
        }
        coords(irp) = CoordRP(coord_r, coord_p);
    });

    const double TOL = 5e-5;
    double max_err = 0.;
    ddc::for_each(domain, [&](IndexRP const irp) {
        const CoordRP coord(coords(irp));
        const CoordXY discrete_coord = mapping(coord);
        const CoordXY analytical_coord = analytical_mapping(coord);

        EXPECT_NEAR(ddc::get<RDimX>(discrete_coord), ddc::get<RDimX>(analytical_coord), TOL);
        EXPECT_NEAR(ddc::get<RDimY>(discrete_coord), ddc::get<RDimY>(analytical_coord), TOL);

        const double diff_x
                = double(ddc::get<RDimX>(discrete_coord) - ddc::get<RDimX>(analytical_coord));
        const double diff_y
                = double(ddc::get<RDimY>(discrete_coord) - ddc::get<RDimY>(analytical_coord));

        max_err = max_err > diff_x ? max_err : diff_x;
        max_err = max_err > diff_y ? max_err : diff_y;
    });

    return max_err;
}


template <class Mapping, int refined_Nr, int refined_Nt>
double test_on_grid_and_not_on_grid(
        int const Nr,
        int const Nt,
        Mapping const& analytical_mapping,
        RefinedDiscreteToCartesian<
                RDimX,
                RDimY,
                SplineRPBuilder,
                SplineRPEvaluator,
                refined_Nr,
                refined_Nt> const& refined_mapping,
        DiscreteToCartesian<RDimX, RDimY, SplineRPBuilder, SplineRPEvaluator> const&
                discrete_mapping,
        IDomainRP const& grid)
{
    std::cout << std::endl
              << "DOMAIN Nr x Nt = " << Nr << " x " << Nt
              << " - REFINED DOMAIN Nr x Nt = " << refined_Nr << " x " << refined_Nt << ":"
              << std::endl;

    std::cout << "           Errors discrete mapping   | Errors refined mapping" << std::endl;
    std::cout << "on grid:     " << check_value_on_grid(discrete_mapping, analytical_mapping, grid)
              << "             |  "
              << check_value_on_grid(refined_mapping, analytical_mapping, grid) << std::endl;

    double const not_on_grid_refined
            = check_value_not_on_grid(refined_mapping, analytical_mapping, grid);
    std::cout << "not on grid: "
              << check_value_not_on_grid(discrete_mapping, analytical_mapping, grid)
              << "             |  " << not_on_grid_refined << std::endl;

    return not_on_grid_refined;
}



/**
 * @brief Compare the Jacobian matrix given by the analytical mapping and
 * the other mapping on the grid points.
 *
 * The error tolerance is given at 5e-10.
 *
 * @param[in] mapping
 *          The mapping we are testing.
 * @param[in] analytical_mappping
 *          The mapping analytically defined.
 * @param[in] domain
 *          The domain on which we test the values.
 */
template <class Mapping, class DiscreteMapping>
double check_Jacobian_on_grid(
        DiscreteMapping const& mapping,
        Mapping const& analytical_mapping,
        IDomainRP const& domain)
{
    double max_err = 0.;
    ddc::for_each(domain, [&](IndexRP const irp) {
        const CoordRP coord(ddc::coordinate(irp));

        Matrix_2x2 discrete_Jacobian;
        Matrix_2x2 analytical_Jacobian;
        mapping.jacobian_matrix(coord, discrete_Jacobian);
        analytical_mapping.jacobian_matrix(coord, analytical_Jacobian);

        const double diff_11 = double(discrete_Jacobian[0][0] - analytical_Jacobian[0][0]);
        const double diff_12 = double(discrete_Jacobian[0][1] - analytical_Jacobian[0][1]);
        const double diff_21 = double(discrete_Jacobian[1][0] - analytical_Jacobian[1][0]);
        const double diff_22 = double(discrete_Jacobian[1][1] - analytical_Jacobian[1][1]);

        max_err = max_err > diff_11 ? max_err : diff_11;
        max_err = max_err > diff_12 ? max_err : diff_12;
        max_err = max_err > diff_21 ? max_err : diff_21;
        max_err = max_err > diff_22 ? max_err : diff_22;
    });

    return max_err;
}



/**
 * @brief Compare the Jacobian matrix given by the analytical mapping and
 * the other mapping not on the grid points.
 *
 * The error tolerance is given at 5e-6.
 * The expected convergence order not on the grid points is
 * d + 1 where d is the degree of B-splines.
 *
 * @param[in] mapping
 *          The mapping we are testing.
 * @param[in] analytical_mappping
 *          The mapping analytically defined.
 * @param[in] domain
 *          The domain on which we test the values.
 */
template <class Mapping, class DiscreteMapping>
double check_Jacobian_not_on_grid(
        DiscreteMapping const& mapping,
        Mapping const& analytical_mapping,
        IDomainRP const& domain)
{
    std::srand(100);

    FieldRP<CoordRP> coords(domain);
    IndexR ir_max(ddc::select<IDimR>(domain).back());
    IndexP ip_max(ddc::select<IDimP>(domain).back());
    ddc::for_each(domain, [&](IndexRP const irp) {
        IndexR ir(ddc::select<IDimR>(irp));
        CoordR coordr_0 = ddc::coordinate(ir);
        CoordR coordr_1 = ddc::coordinate(ir + 1);
        double coord_r;
        if (ir.uid() < ir_max.uid()) {
            double factor = double(std::rand()) / RAND_MAX;
            coord_r = coordr_0 + (coordr_1 - coordr_0) * factor;
        } else {
            coord_r = coordr_0;
        }

        IndexP ip(ddc::select<IDimP>(irp));
        CoordP coordp_0 = ddc::coordinate(ip);
        CoordP coordp_1 = ddc::coordinate(ip + 1);
        double coord_p;
        if (ip.uid() < ip_max.uid()) {
            double factor = double(std::rand()) / RAND_MAX;
            coord_p = coordp_0 + (coordp_1 - coordp_0) * factor;
        } else {
            coord_p = coordp_0;
        }
        coords(irp) = CoordRP(coord_r, coord_p);
    });


    double max_err = 0.;
    ddc::for_each(domain, [&](IndexRP const irp) {
        const CoordRP coord(coords(irp));

        Matrix_2x2 discrete_Jacobian;
        Matrix_2x2 analytical_Jacobian;
        mapping.jacobian_matrix(coord, discrete_Jacobian);
        analytical_mapping.jacobian_matrix(coord, analytical_Jacobian);

        const double diff_11 = double(discrete_Jacobian[0][0] - analytical_Jacobian[0][0]);
        const double diff_12 = double(discrete_Jacobian[0][1] - analytical_Jacobian[0][1]);
        const double diff_21 = double(discrete_Jacobian[1][0] - analytical_Jacobian[1][0]);
        const double diff_22 = double(discrete_Jacobian[1][1] - analytical_Jacobian[1][1]);

        max_err = max_err > diff_11 ? max_err : diff_11;
        max_err = max_err > diff_12 ? max_err : diff_12;
        max_err = max_err > diff_21 ? max_err : diff_21;
        max_err = max_err > diff_22 ? max_err : diff_22;
    });

    return max_err;
}



template <class Mapping, int refined_Nr, int refined_Nt>
double test_Jacobian(
        int const Nr,
        int const Nt,
        Mapping const& analytical_mapping,
        RefinedDiscreteToCartesian<
                RDimX,
                RDimY,
                SplineRPBuilder,
                SplineRPEvaluator,
                refined_Nr,
                refined_Nt> const& refined_mapping,
        DiscreteToCartesian<RDimX, RDimY, SplineRPBuilder, SplineRPEvaluator> const&
                discrete_mapping,
        IDomainRP const& grid)
{
    std::cout << std::endl
              << "DOMAIN Nr x Nt = " << Nr << " x " << Nt
              << " - REFINED DOMAIN Nr x Nt = " << refined_Nr << " x " << refined_Nt << ":"
              << std::endl;

    std::cout << "           Errors discrete mapping   | Errors refined mapping" << std::endl;
    std::cout << "on grid:     "
              << check_Jacobian_on_grid(discrete_mapping, analytical_mapping, grid)
              << "                 |  ";
    std::cout << check_Jacobian_on_grid(refined_mapping, analytical_mapping, grid) << std::endl;

    double const not_on_grid_refined
            = check_Jacobian_not_on_grid(refined_mapping, analytical_mapping, grid);
    std::cout << "not on grid: "
              << check_Jacobian_not_on_grid(discrete_mapping, analytical_mapping, grid);
    std::cout << "                 |  " << not_on_grid_refined << std::endl;

    return not_on_grid_refined;
}



/**
 * @brief Compare the pseudo Cartesian Jacobian matrix given
 * by the analytical mapping and the other mapping at the O-point.
 *
 * The error tolerance is given at 5e-10.
 *
 * @param[in] mapping
 *          The mapping we are testing.
 * @param[in] analytical_mappping
 *          The mapping analytically defined.
 * @param[in] domain
 *          The domain on which we test the values.
 */
template <class Mapping, class DiscreteMapping>
double check_pseudo_Cart(
        DiscreteMapping const& mapping,
        Mapping const& analytical_mapping,
        IDomainRP const& domain)
{
    double max_err = 0.;
    Matrix_2x2 discrete_pseudo_Cart;
    Matrix_2x2 analytical_pseudo_Cart;
    mapping.to_pseudo_cartesian_jacobian_center_matrix(domain, discrete_pseudo_Cart);
    analytical_mapping.to_pseudo_cartesian_jacobian_center_matrix(domain, analytical_pseudo_Cart);

    const double diff_11 = double(discrete_pseudo_Cart[0][0] - analytical_pseudo_Cart[0][0]);
    const double diff_12 = double(discrete_pseudo_Cart[0][1] - analytical_pseudo_Cart[0][1]);
    const double diff_21 = double(discrete_pseudo_Cart[1][0] - analytical_pseudo_Cart[1][0]);
    const double diff_22 = double(discrete_pseudo_Cart[1][1] - analytical_pseudo_Cart[1][1]);

    max_err = max_err > diff_11 ? max_err : diff_11;
    max_err = max_err > diff_12 ? max_err : diff_12;
    max_err = max_err > diff_21 ? max_err : diff_21;
    max_err = max_err > diff_22 ? max_err : diff_22;

    return max_err;
}



template <class Mapping, int refined_Nr, int refined_Nt>
double test_pseudo_Cart(
        int const Nr,
        int const Nt,
        Mapping const& analytical_mapping,
        RefinedDiscreteToCartesian<
                RDimX,
                RDimY,
                SplineRPBuilder,
                SplineRPEvaluator,
                refined_Nr,
                refined_Nt> const& refined_mapping,
        DiscreteToCartesian<RDimX, RDimY, SplineRPBuilder, SplineRPEvaluator> const&
                discrete_mapping,
        IDomainRP const& grid)
{
    std::cout << std::endl
              << "DOMAIN Nr x Nt = " << Nr << " x " << Nt
              << " - REFINED DOMAIN Nr x Nt = " << refined_Nr << " x " << refined_Nt << ":"
              << std::endl;

    double const pseudo_Cart_err = check_pseudo_Cart(refined_mapping, analytical_mapping, grid);
    std::cout << "           Errors discrete mapping   | Errors refined mapping" << std::endl;
    std::cout << "on grid:     " << check_pseudo_Cart(discrete_mapping, analytical_mapping, grid)
              << "                 |  ";
    std::cout << pseudo_Cart_err << std::endl;

    return pseudo_Cart_err;
}



} // namespace



TEST(RefinedDiscreteMapping, TestRefinedDiscreteMapping)
{
    const CzarnyMapping analytical_mapping(0.3, 1.4);
    //const CircularMapping analytical_mapping;

    // Discrete domain ---
    int const Nr = 16;
    int const Nt = 32;

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

    for (int i(1); i < r_size + 1; ++i) {
        r_knots[i] = CoordR(r_min + i * dr);
    }
    r_knots[0] = CoordR(r_min);
    r_knots[r_size] = CoordR(r_max);
    for (int i(1); i < p_size; ++i) {
        p_knots[i] = CoordP(p_min + i * dp);
    }
    p_knots[0] = CoordP(p_min);
    p_knots[p_size] = CoordP(p_max);

    ddc::init_discrete_space<BSplinesR>(r_knots);
    ddc::init_discrete_space<BSplinesP>(p_knots);

    ddc::init_discrete_space<IDimR>(SplineInterpPointsR::get_sampling());
    ddc::init_discrete_space<IDimP>(SplineInterpPointsP::get_sampling());

    IDomainR interpolation_domain_R(SplineInterpPointsR::get_domain());
    IDomainP interpolation_domain_P(SplineInterpPointsP::get_domain());
    IDomainRP grid(interpolation_domain_R, interpolation_domain_P);


    // Operators ---
    SplineRPBuilder builder(grid);
    ddc::ConstantExtrapolationRule<RDimR, RDimP> boundary_condition_r_left(r_min);
    ddc::ConstantExtrapolationRule<RDimR, RDimP> boundary_condition_r_right(r_max);
    SplineRPEvaluator spline_evaluator(
            boundary_condition_r_left,
            boundary_condition_r_right,
            ddc::PeriodicExtrapolationRule<RDimP>(),
            ddc::PeriodicExtrapolationRule<RDimP>());


    // Tests ---
    std::array<double, 3> results;

    DiscreteToCartesian discrete_mapping
            = DiscreteToCartesian<RDimX, RDimY, SplineRPBuilder, SplineRPEvaluator>::
                    analytical_to_discrete(analytical_mapping, builder, spline_evaluator);

    RefinedDiscreteToCartesian refined_mapping_16x32
            = RefinedDiscreteToCartesian<RDimX, RDimY, SplineRPBuilder, SplineRPEvaluator, 16, 32>::
                    analytical_to_refined(analytical_mapping, grid);

    RefinedDiscreteToCartesian refined_mapping_32x64
            = RefinedDiscreteToCartesian<RDimX, RDimY, SplineRPBuilder, SplineRPEvaluator, 32, 64>::
                    analytical_to_refined(analytical_mapping, grid);


    RefinedDiscreteToCartesian refined_mapping_64x128 = RefinedDiscreteToCartesian<
            RDimX,
            RDimY,
            SplineRPBuilder,
            SplineRPEvaluator,
            64,
            128>::analytical_to_refined(analytical_mapping, grid);


    std::cout << std::endl
              << "TESTS ON THE MAPPING VALUES: -------------------------------" << std::endl;
    results[0] = test_on_grid_and_not_on_grid(
            Nr,
            Nt,
            analytical_mapping,
            refined_mapping_16x32,
            discrete_mapping,
            grid);

    results[1] = test_on_grid_and_not_on_grid(
            Nr,
            Nt,
            analytical_mapping,
            refined_mapping_32x64,
            discrete_mapping,
            grid);


    results[2] = test_on_grid_and_not_on_grid(
            Nr,
            Nt,
            analytical_mapping,
            refined_mapping_64x128,
            discrete_mapping,
            grid);


    std::cout << std::endl << "Convergence order : " << std::endl << "  -" << std::endl;
    for (std::size_t i(0); i < results.size() - 1; i++) {
        double const order = std::log(results[i] / results[i + 1]) / std::log(2);
        std::cout << "  " << order << std::endl;

        int const BSDegree = 3;
        EXPECT_NEAR(order, BSDegree + 1, 0.25);
    }


    std::cout << std::endl
              << "TESTS ON THE JACOBIAN MATRIX: ------------------------------" << std::endl;
    results[0] = test_Jacobian(
            Nr,
            Nt,
            analytical_mapping,
            refined_mapping_16x32,
            discrete_mapping,
            grid);

    results[1] = test_Jacobian(
            Nr,
            Nt,
            analytical_mapping,
            refined_mapping_32x64,
            discrete_mapping,
            grid);


    results[2] = test_Jacobian(
            Nr,
            Nt,
            analytical_mapping,
            refined_mapping_64x128,
            discrete_mapping,
            grid);


    std::cout << std::endl << "Convergence order : " << std::endl << "  -" << std::endl;
    for (std::size_t i(0); i < results.size() - 1; i++) {
        double const order = std::log(results[i] / results[i + 1]) / std::log(2);
        std::cout << "  " << order << std::endl;

        int const BSDegree = 3;
        EXPECT_NEAR(order, BSDegree, 0.25);
    }


    std::cout << std::endl
              << "TESTS ON THE PSEUDO CARTESIAN MATRIX: ----------------------" << std::endl;
    results[0] = test_pseudo_Cart(
            Nr,
            Nt,
            analytical_mapping,
            refined_mapping_16x32,
            discrete_mapping,
            grid);

    results[1] = test_pseudo_Cart(
            Nr,
            Nt,
            analytical_mapping,
            refined_mapping_32x64,
            discrete_mapping,
            grid);


    results[2] = test_pseudo_Cart(
            Nr,
            Nt,
            analytical_mapping,
            refined_mapping_64x128,
            discrete_mapping,
            grid);


    std::cout << std::endl << "Convergence order : " << std::endl << "  -" << std::endl;
    for (std::size_t i(0); i < results.size() - 1; i++) {
        double const order = std::log(results[i] / results[i + 1]) / std::log(2);
        std::cout << "  " << order << std::endl;

        int const BSDegree = 3;
        EXPECT_NEAR(order, BSDegree + 1, 0.25);
    }
}
