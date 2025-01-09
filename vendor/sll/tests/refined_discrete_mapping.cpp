#include <array>
#include <cassert>
#include <cstdlib>

#include <ddc/ddc.hpp>
#include <ddc/kernels/splines.hpp>

#include "sll/mapping/circular_to_cartesian.hpp"
#include "sll/mapping/combined_mapping.hpp"
#include "sll/mapping/czarny_to_cartesian.hpp"
#include "sll/mapping/discrete_mapping_builder.hpp"
#include "sll/mapping/discrete_to_cartesian.hpp"
#include "sll/mapping/inv_jacobian_o_point.hpp"

#include "test_utils.hpp"



namespace {
struct X
{
};
struct Y
{
};
struct Xpc
{
};
struct Ypc
{
};
struct R
{
    static bool constexpr PERIODIC = false;
};

struct Theta
{
    static bool constexpr PERIODIC = true;
};

using CoordR = ddc::Coordinate<R>;
using CoordTheta = ddc::Coordinate<Theta>;
using CoordRTheta = ddc::Coordinate<R, Theta>;

using CoordXY = ddc::Coordinate<X, Y>;

int constexpr BSDegree = 3;

struct BSplinesR : ddc::NonUniformBSplines<R, BSDegree>
{
};
struct BSplinesTheta : ddc::NonUniformBSplines<Theta, BSDegree>
{
};


using SplineInterpPointsR = ddc::
        GrevilleInterpolationPoints<BSplinesR, ddc::BoundCond::GREVILLE, ddc::BoundCond::GREVILLE>;
using SplineInterpPointsTheta = ddc::GrevilleInterpolationPoints<
        BSplinesTheta,
        ddc::BoundCond::PERIODIC,
        ddc::BoundCond::PERIODIC>;

struct GridR : SplineInterpPointsR::interpolation_discrete_dimension_type
{
};
struct GridTheta : SplineInterpPointsTheta::interpolation_discrete_dimension_type
{
};

using SplineRThetaBuilder_host = ddc::SplineBuilder2D<
        Kokkos::DefaultHostExecutionSpace,
        Kokkos::DefaultHostExecutionSpace::memory_space,
        BSplinesR,
        BSplinesTheta,
        GridR,
        GridTheta,
        ddc::BoundCond::GREVILLE,
        ddc::BoundCond::GREVILLE,
        ddc::BoundCond::PERIODIC,
        ddc::BoundCond::PERIODIC,
        ddc::SplineSolver::LAPACK,
        GridR,
        GridTheta>;

using SplineRThetaEvaluator = ddc::SplineEvaluator2D<
        Kokkos::DefaultHostExecutionSpace,
        Kokkos::DefaultHostExecutionSpace::memory_space,
        BSplinesR,
        BSplinesTheta,
        GridR,
        GridTheta,
        ddc::ConstantExtrapolationRule<R, Theta>,
        ddc::ConstantExtrapolationRule<R, Theta>,
        ddc::PeriodicExtrapolationRule<Theta>,
        ddc::PeriodicExtrapolationRule<Theta>,
        GridR,
        GridTheta>;

using BSIdxRangeR = ddc::DiscreteDomain<BSplinesR>;
using BSIdxRangeTheta = ddc::DiscreteDomain<BSplinesTheta>;
using BSIdxRangeRTheta = ddc::DiscreteDomain<BSplinesR, BSplinesTheta>;

using IdxRangeR = ddc::DiscreteDomain<GridR>;
using IdxRangeTheta = ddc::DiscreteDomain<GridTheta>;
using IdxRangeRTheta = ddc::DiscreteDomain<GridR, GridTheta>;

using IdxR = ddc::DiscreteElement<GridR>;
using IdxTheta = ddc::DiscreteElement<GridTheta>;
using IdxRTheta = ddc::DiscreteElement<GridR, GridTheta>;

using IdxStepR = ddc::DiscreteVector<GridR>;
using IdxStepTheta = ddc::DiscreteVector<GridTheta>;
using IdxStepRTheta = ddc::DiscreteVector<GridR, GridTheta>;

template <class ElementType>
using FieldMemR = ddc::Chunk<ElementType, IdxRangeR>;

template <class ElementType>
using FieldMemTheta = ddc::Chunk<ElementType, IdxRangeTheta>;

template <class ElementType>
using FieldMemRTheta = ddc::Chunk<ElementType, IdxRangeRTheta>;



using IdxRangeRTheta = ddc::DiscreteDomain<GridR, GridTheta>;


using CzarnyMapping = CzarnyToCartesian<R, Theta, X, Y>;
using CircularMapping = CircularToCartesian<R, Theta, X, Y>;


template <class ElementType>
using FieldMemRTheta = ddc::Chunk<ElementType, IdxRangeRTheta>;

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
 * @param[in] idx_range
 *          The index range on which we test the values.
 */
template <class Mapping, class DiscreteMapping>
double check_value_on_grid(
        DiscreteMapping const& mapping,
        Mapping const& analytical_mapping,
        IdxRangeRTheta const& idx_range)
{
    const double TOL = 1e-7;
    double max_err = 0.;
    ddc::for_each(idx_range, [&](IdxRTheta const irp) {
        const CoordRTheta coord(ddc::coordinate(irp));
        const CoordXY discrete_coord = mapping(coord);
        const CoordXY analytical_coord = analytical_mapping(coord);

        EXPECT_NEAR(ddc::get<X>(discrete_coord), ddc::get<X>(analytical_coord), TOL);
        EXPECT_NEAR(ddc::get<Y>(discrete_coord), ddc::get<Y>(analytical_coord), TOL);

        const double diff_x = double(ddc::get<X>(discrete_coord) - ddc::get<X>(analytical_coord));
        const double diff_y = double(ddc::get<Y>(discrete_coord) - ddc::get<Y>(analytical_coord));

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
 * @param[in] idx_range
 *          The index range on which we test the values.
 */
template <class Mapping, class DiscreteMapping>
double check_value_not_on_grid(
        DiscreteMapping const& mapping,
        Mapping const& analytical_mapping,
        IdxRangeRTheta const& idx_range)
{
    std::srand(100);

    FieldMemRTheta<CoordRTheta> coords(idx_range);
    IdxR ir_max(ddc::select<GridR>(idx_range).back());
    IdxTheta itheta_max(ddc::select<GridTheta>(idx_range).back());
    ddc::for_each(idx_range, [&](IdxRTheta const irp) {
        IdxR ir(ddc::select<GridR>(irp));
        CoordR coord_r_0 = ddc::coordinate(ir);
        CoordR coord_r_1 = ddc::coordinate(ir + 1);
        double coord_r;
        if (ir.uid() < ir_max.uid()) {
            double factor = double(std::rand()) / RAND_MAX;
            coord_r = coord_r_0 + (coord_r_1 - coord_r_0) * factor;
        } else {
            coord_r = coord_r_0;
        }

        IdxTheta ip(ddc::select<GridTheta>(irp));
        CoordTheta coord_theta_0 = ddc::coordinate(ip);
        CoordTheta coord_theta_1 = ddc::coordinate(ip + 1);
        double coord_theta;
        if (ip.uid() < itheta_max.uid()) {
            double factor = double(std::rand()) / RAND_MAX;
            coord_theta = coord_theta_0 + (coord_theta_1 - coord_theta_0) * factor;
        } else {
            coord_theta = coord_theta_0;
        }
        coords(irp) = CoordRTheta(coord_r, coord_theta);
    });

    const double TOL = 5e-5;
    double max_err = 0.;
    ddc::for_each(idx_range, [&](IdxRTheta const irp) {
        const CoordRTheta coord(coords(irp));
        const CoordXY discrete_coord = mapping(coord);
        const CoordXY analytical_coord = analytical_mapping(coord);

        EXPECT_NEAR(ddc::get<X>(discrete_coord), ddc::get<X>(analytical_coord), TOL);
        EXPECT_NEAR(ddc::get<Y>(discrete_coord), ddc::get<Y>(analytical_coord), TOL);

        const double diff_x = double(ddc::get<X>(discrete_coord) - ddc::get<X>(analytical_coord));
        const double diff_y = double(ddc::get<Y>(discrete_coord) - ddc::get<Y>(analytical_coord));

        max_err = max_err > diff_x ? max_err : diff_x;
        max_err = max_err > diff_y ? max_err : diff_y;
    });

    return max_err;
}


template <class Mapping, class DiscreteMapping, class RefinedDiscreteMapping>
double test_on_grid_and_not_on_grid(
        int const Nr,
        int const Nt,
        int const refined_Nr,
        int const refined_Nt,
        Mapping const& analytical_mapping,
        RefinedDiscreteMapping const& refined_mapping,
        DiscreteMapping const& discrete_mapping,
        IdxRangeRTheta const& grid)
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
 * @param[in] idx_range
 *          The index range on which we test the values.
 */
template <class Mapping, class DiscreteMapping>
double check_Jacobian_on_grid(
        DiscreteMapping const& mapping,
        Mapping const& analytical_mapping,
        IdxRangeRTheta const& idx_range)
{
    static_assert(has_2d_jacobian_v<Mapping, CoordRTheta>);
    double max_err = 0.;
    ddc::for_each(idx_range, [&](IdxRTheta const irp) {
        const CoordRTheta coord(ddc::coordinate(irp));

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
 * @param[in] idx_range
 *          The index range on which we test the values.
 */
template <class Mapping, class DiscreteMapping>
double check_Jacobian_not_on_grid(
        DiscreteMapping const& mapping,
        Mapping const& analytical_mapping,
        IdxRangeRTheta const& idx_range)
{
    std::srand(100);

    FieldMemRTheta<CoordRTheta> coords(idx_range);
    IdxR ir_max(ddc::select<GridR>(idx_range).back());
    IdxTheta itheta_max(ddc::select<GridTheta>(idx_range).back());
    ddc::for_each(idx_range, [&](IdxRTheta const irp) {
        IdxR ir(ddc::select<GridR>(irp));
        CoordR coord_r_0 = ddc::coordinate(ir);
        CoordR coord_r_1 = ddc::coordinate(ir + 1);
        double coord_r;
        if (ir.uid() < ir_max.uid()) {
            double factor = double(std::rand()) / RAND_MAX;
            coord_r = coord_r_0 + (coord_r_1 - coord_r_0) * factor;
        } else {
            coord_r = coord_r_0;
        }

        IdxTheta ip(ddc::select<GridTheta>(irp));
        CoordTheta coord_theta_0 = ddc::coordinate(ip);
        CoordTheta coord_theta_1 = ddc::coordinate(ip + 1);
        double coord_theta;
        if (ip.uid() < itheta_max.uid()) {
            double factor = double(std::rand()) / RAND_MAX;
            coord_theta = coord_theta_0 + (coord_theta_1 - coord_theta_0) * factor;
        } else {
            coord_theta = coord_theta_0;
        }
        coords(irp) = CoordRTheta(coord_r, coord_theta);
    });


    double max_err = 0.;
    ddc::for_each(idx_range, [&](IdxRTheta const irp) {
        const CoordRTheta coord(coords(irp));

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



template <class Mapping, class RefinedDiscreteMapping, class DiscreteMapping>
double test_Jacobian(
        int const Nr,
        int const Nt,
        int const refined_Nr,
        int const refined_Nt,
        Mapping const& analytical_mapping,
        RefinedDiscreteMapping const& refined_mapping,
        DiscreteMapping const& discrete_mapping,
        IdxRangeRTheta const& grid)
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
 * @param[in] idx_range
 *          The index range on which we test the values.
 */
template <class Mapping, class DiscreteMapping>
double check_pseudo_Cart(
        DiscreteMapping const& mapping,
        Mapping const& analytical_mapping,
        IdxRangeRTheta const& idx_range)
{
    const double epsilon(1.e-16);
    using PseudoCartToCirc = CartesianToCircular<Xpc, Ypc, R, Theta>;
    PseudoCartToCirc pseudo_cart_to_circ;
    CombinedMapping<DiscreteMapping, PseudoCartToCirc>
            pseudo_cart_to_cart(mapping, pseudo_cart_to_circ, epsilon);
    CombinedMapping<Mapping, PseudoCartToCirc>
            analytical_pseudo_cart_to_cart(analytical_mapping, pseudo_cart_to_circ, epsilon);
    InvJacobianOPoint<CombinedMapping<DiscreteMapping, PseudoCartToCirc>, CoordRTheta>
            inv_jacob_mapping(pseudo_cart_to_cart);
    InvJacobianOPoint<CombinedMapping<Mapping, PseudoCartToCirc>, CoordRTheta>
            inv_jacob_analytical_mapping(analytical_pseudo_cart_to_cart);

    Matrix_2x2 discrete_pseudo_Cart = inv_jacob_mapping();
    Matrix_2x2 analytical_pseudo_Cart = inv_jacob_analytical_mapping();

    const double diff_11 = double(discrete_pseudo_Cart[0][0] - analytical_pseudo_Cart[0][0]);
    const double diff_12 = double(discrete_pseudo_Cart[0][1] - analytical_pseudo_Cart[0][1]);
    const double diff_21 = double(discrete_pseudo_Cart[1][0] - analytical_pseudo_Cart[1][0]);
    const double diff_22 = double(discrete_pseudo_Cart[1][1] - analytical_pseudo_Cart[1][1]);

    double max_err = 0.;
    max_err = max_err > diff_11 ? max_err : diff_11;
    max_err = max_err > diff_12 ? max_err : diff_12;
    max_err = max_err > diff_21 ? max_err : diff_21;
    max_err = max_err > diff_22 ? max_err : diff_22;

    return max_err;
}



template <class Mapping, class RefinedDiscreteMapping, class DiscreteMapping>
double test_pseudo_Cart(
        int const Nr,
        int const Nt,
        int const refined_Nr,
        int const refined_Nt,
        Mapping const& analytical_mapping,
        RefinedDiscreteMapping const& refined_mapping,
        DiscreteMapping const& discrete_mapping,
        IdxRangeRTheta const& grid)
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

    // Discrete index range ---
    int const Nr = 16;
    int const Nt = 32;

    CoordR const r_min(0.0);
    CoordR const r_max(1.0);
    IdxStepR const r_size(Nr);

    CoordTheta const theta_min(0.0);
    CoordTheta const theta_max(2.0 * M_PI);
    IdxStepTheta const theta_size(Nt);

    double const dr((r_max - r_min) / r_size);
    double const dp((theta_max - theta_min) / theta_size);

    std::vector<CoordR> r_knots(r_size + 1);
    std::vector<CoordTheta> theta_knots(theta_size + 1);

    for (int i(1); i < r_size + 1; ++i) {
        r_knots[i] = CoordR(r_min + i * dr);
    }
    r_knots[0] = CoordR(r_min);
    r_knots[r_size] = CoordR(r_max);
    for (int i(1); i < theta_size; ++i) {
        theta_knots[i] = CoordTheta(theta_min + i * dp);
    }
    theta_knots[0] = CoordTheta(theta_min);
    theta_knots[theta_size] = CoordTheta(theta_max);

    ddc::init_discrete_space<BSplinesR>(r_knots);
    ddc::init_discrete_space<BSplinesTheta>(theta_knots);

    ddc::init_discrete_space<GridR>(SplineInterpPointsR::get_sampling<GridR>());
    ddc::init_discrete_space<GridTheta>(SplineInterpPointsTheta::get_sampling<GridTheta>());

    IdxRangeR interpolation_idx_range_R(SplineInterpPointsR::get_domain<GridR>());
    IdxRangeTheta interpolation_idx_range_Theta(SplineInterpPointsTheta::get_domain<GridTheta>());
    IdxRangeRTheta grid(interpolation_idx_range_R, interpolation_idx_range_Theta);


    // Operators ---
    SplineRThetaBuilder_host builder(grid);
    ddc::ConstantExtrapolationRule<R, Theta> boundary_condition_r_left(r_min);
    ddc::ConstantExtrapolationRule<R, Theta> boundary_condition_r_right(r_max);
    SplineRThetaEvaluator spline_evaluator(
            boundary_condition_r_left,
            boundary_condition_r_right,
            ddc::PeriodicExtrapolationRule<Theta>(),
            ddc::PeriodicExtrapolationRule<Theta>());


    // Tests ---
    std::array<double, 3> results;

    DiscreteToCartesianBuilder<X, Y, SplineRThetaBuilder_host, SplineRThetaEvaluator>
            mapping_builder(
                    Kokkos::DefaultHostExecutionSpace(),
                    analytical_mapping,
                    builder,
                    spline_evaluator);
    DiscreteToCartesian discrete_mapping = mapping_builder();

    RefinedDiscreteToCartesianBuilder<X, Y, SplineRThetaBuilder_host, SplineRThetaEvaluator, 16, 32>
            mapping_builder_16x32(
                    Kokkos::DefaultHostExecutionSpace(),
                    analytical_mapping,
                    builder,
                    spline_evaluator);
    DiscreteToCartesian refined_mapping_16x32 = mapping_builder_16x32();

    RefinedDiscreteToCartesianBuilder<X, Y, SplineRThetaBuilder_host, SplineRThetaEvaluator, 32, 64>
            mapping_builder_32x64(
                    Kokkos::DefaultHostExecutionSpace(),
                    analytical_mapping,
                    builder,
                    spline_evaluator);
    DiscreteToCartesian refined_mapping_32x64 = mapping_builder_32x64();

    RefinedDiscreteToCartesianBuilder<
            X,
            Y,
            SplineRThetaBuilder_host,
            SplineRThetaEvaluator,
            64,
            128>
            mapping_builder_64x128(
                    Kokkos::DefaultHostExecutionSpace(),
                    analytical_mapping,
                    builder,
                    spline_evaluator);
    DiscreteToCartesian refined_mapping_64x128 = mapping_builder_64x128();


    std::cout << std::endl
              << "TESTS ON THE MAPPING VALUES: -------------------------------" << std::endl;
    results[0] = test_on_grid_and_not_on_grid(
            Nr,
            Nt,
            16,
            32,
            analytical_mapping,
            refined_mapping_16x32,
            discrete_mapping,
            grid);

    results[1] = test_on_grid_and_not_on_grid(
            Nr,
            Nt,
            32,
            64,
            analytical_mapping,
            refined_mapping_32x64,
            discrete_mapping,
            grid);


    results[2] = test_on_grid_and_not_on_grid(
            Nr,
            Nt,
            64,
            128,
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
            16,
            32,
            analytical_mapping,
            refined_mapping_16x32,
            discrete_mapping,
            grid);

    results[1] = test_Jacobian(
            Nr,
            Nt,
            32,
            64,
            analytical_mapping,
            refined_mapping_32x64,
            discrete_mapping,
            grid);


    results[2] = test_Jacobian(
            Nr,
            Nt,
            64,
            128,
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
            16,
            32,
            analytical_mapping,
            refined_mapping_16x32,
            discrete_mapping,
            grid);

    results[1] = test_pseudo_Cart(
            Nr,
            Nt,
            32,
            64,
            analytical_mapping,
            refined_mapping_32x64,
            discrete_mapping,
            grid);


    results[2] = test_pseudo_Cart(
            Nr,
            Nt,
            64,
            128,
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
