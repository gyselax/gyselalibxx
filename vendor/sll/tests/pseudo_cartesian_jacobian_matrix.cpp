#include <ddc/ddc.hpp>
#include <ddc/kernels/splines.hpp>

#include <sll/mapping/circular_to_cartesian.hpp>
#include <sll/mapping/combined_mapping.hpp>
#include <sll/mapping/czarny_to_cartesian.hpp>
#include <sll/mapping/discrete_mapping_builder.hpp>
#include <sll/mapping/discrete_to_cartesian.hpp>
#include <sll/math_tools.hpp>

#include <gtest/gtest.h>


namespace {
/**
 * @brief A class for the Google test.
 */
template <int N>
class PseudoCartesianJacobianMatrixTest
{
public:
    struct X
    {
    };
    struct Y
    {
    };
    struct X_pc
    {
    };
    struct Y_pc
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

    static int constexpr BSDegree = 3;

    struct BSplinesR : ddc::NonUniformBSplines<R, BSDegree>
    {
    };
    struct BSplinesP : ddc::NonUniformBSplines<Theta, BSDegree>
    {
    };


    using InterpPointsR = ddc::GrevilleInterpolationPoints<
            BSplinesR,
            ddc::BoundCond::GREVILLE,
            ddc::BoundCond::GREVILLE>;
    using InterpPointsP = ddc::GrevilleInterpolationPoints<
            BSplinesP,
            ddc::BoundCond::PERIODIC,
            ddc::BoundCond::PERIODIC>;


    struct GridR : InterpPointsR::interpolation_discrete_dimension_type
    {
    };
    struct GridP : InterpPointsP::interpolation_discrete_dimension_type
    {
    };


    using BSIdxRangeR = ddc::DiscreteDomain<BSplinesR>;
    using BSIdxRangeP = ddc::DiscreteDomain<BSplinesP>;
    using BSIdxRangeRP = ddc::DiscreteDomain<BSplinesR, BSplinesP>;


    using IdxR = ddc::DiscreteElement<GridR>;
    using IdxP = ddc::DiscreteElement<GridP>;
    using IdxRP = ddc::DiscreteElement<GridR, GridP>;


    using IdxStepR = ddc::DiscreteVector<GridR>;
    using IdxStepP = ddc::DiscreteVector<GridP>;
    using IdxStepRP = ddc::DiscreteVector<GridR, GridP>;


    using IdxRangeR = ddc::DiscreteDomain<GridR>;
    using IdxRangeP = ddc::DiscreteDomain<GridP>;
    using IdxRangeRP = ddc::DiscreteDomain<GridR, GridP>;


    using SplineRThetaBuilder_host = ddc::SplineBuilder2D<
            Kokkos::DefaultHostExecutionSpace,
            Kokkos::DefaultHostExecutionSpace::memory_space,
            BSplinesR,
            BSplinesP,
            GridR,
            GridP,
            ddc::BoundCond::GREVILLE,
            ddc::BoundCond::GREVILLE,
            ddc::BoundCond::PERIODIC,
            ddc::BoundCond::PERIODIC,
            ddc::SplineSolver::LAPACK,
            GridR,
            GridP>;

    using SplineRThetaEvaluator = ddc::SplineEvaluator2D<
            Kokkos::DefaultHostExecutionSpace,
            Kokkos::DefaultHostExecutionSpace::memory_space,
            BSplinesR,
            BSplinesP,
            GridR,
            GridP,
            ddc::NullExtrapolationRule,
            ddc::NullExtrapolationRule,
            ddc::PeriodicExtrapolationRule<Theta>,
            ddc::PeriodicExtrapolationRule<Theta>,
            GridR,
            GridP>;


    using spline_idx_range = ddc::DiscreteDomain<BSplinesR, BSplinesP>;

    template <class ElementType>
    using FieldMemRP = ddc::Chunk<ElementType, IdxRangeRP>;

    using Matrix_2x2 = std::array<std::array<double, 2>, 2>;

    int const Nr;
    int const Nt;

public:
    /**
     * Instantiate a PseudoCartesianJacobianMatrixTest object.
     */
    PseudoCartesianJacobianMatrixTest() : Nr(N), Nt(2 * N) {};


    /**
     * @brief Test the pseudo Cartesian Jacobian matrix for a discrete
     * mapping of a circular mapping and a discrete mapping of a Czarny
     * mapping.
     *
     * Create a Nr x Nt discrete index range with Nr = N, Nt = 2*N and N a
     * templated parameter.
     * Then compute the infinity norm of the difference between the
     * the pseudo Cartesian Jacobian matrix of the analytical mapping
     * and the pseudo Cartesian Jacobian matrix of the discrete mapping
     * for the circular and Czarny mapping.
     */
    std::array<double, 2> test_circ_and_czar()
    {
        // --- Define the grid. ---------------------------------------------------------------------------
        CoordR const r_min(0.0);
        CoordR const r_max(1.0);
        IdxStepR const r_size(Nr);

        CoordTheta const p_min(0.0);
        CoordTheta const p_max(2.0 * M_PI);
        IdxStepP const p_size(Nt);

        double const dr((r_max - r_min) / r_size);
        double const dp((p_max - p_min) / p_size);


        std::vector<CoordR> r_knots(r_size + 1);
        std::vector<CoordTheta> p_knots(p_size + 1);


        r_knots[0] = r_min;
        for (int i(1); i < r_size; ++i) {
            r_knots[i] = CoordR(r_min + i * dr);
        }
        r_knots[r_size] = r_max;
        for (int i(0); i < p_size + 1; ++i) {
            p_knots[i] = CoordTheta(p_min + i * dp);
        }

        // Creating mesh & supports
        ddc::init_discrete_space<BSplinesR>(r_knots);
        ddc::init_discrete_space<BSplinesP>(p_knots);

        ddc::init_discrete_space<GridR>(InterpPointsR::template get_sampling<GridR>());
        ddc::init_discrete_space<GridP>(InterpPointsP::template get_sampling<GridP>());

        IdxRangeR interpolation_idx_range_R(InterpPointsR::template get_domain<GridR>());
        IdxRangeP interpolation_idx_range_P(InterpPointsP::template get_domain<GridP>());
        IdxRangeRP grid(interpolation_idx_range_R, interpolation_idx_range_P);

        // --- Define the operators. ----------------------------------------------------------------------
        SplineRThetaBuilder_host const builder(grid);
        ddc::NullExtrapolationRule r_extrapolation_rule;
        ddc::PeriodicExtrapolationRule<Theta> p_extrapolation_rule;
        SplineRThetaEvaluator spline_evaluator(
                r_extrapolation_rule,
                r_extrapolation_rule,
                p_extrapolation_rule,
                p_extrapolation_rule);

        Matrix_2x2 analytical_matrix;
        Matrix_2x2 discrete_matrix;


        // --- CIRCULAR MAPPING ---------------------------------------------------------------------------
        std::cout << " - Nr x Nt  = " << Nr << " x " << Nt << std::endl
                  << "   - Circular mapping: ";
        const CircularToCartesian<R, Theta, X, Y> circ_to_cart;
        const CartesianToCircular<X_pc, Y_pc, R, Theta> pseudo_cart_to_circ;
        using PseudoCartToCircToCart = CombinedMapping<
                CircularToCartesian<R, Theta, X, Y>,
                CartesianToCircular<X_pc, Y_pc, R, Theta>>;
        const PseudoCartToCircToCart
                pseudo_cart_to_circ_to_cart(circ_to_cart, pseudo_cart_to_circ, 1e-12);
        DiscreteToCartesianBuilder<X, Y, SplineRThetaBuilder_host, SplineRThetaEvaluator>
                mapping_builder_circ(
                        Kokkos::DefaultHostExecutionSpace(),
                        circ_to_cart,
                        builder,
                        spline_evaluator);
        DiscreteToCartesian discrete_mapping_circ_to_cart = mapping_builder_circ();
        using DiscreteMappingCirc = CombinedMapping<
                decltype(discrete_mapping_circ_to_cart),
                CartesianToCircular<X_pc, Y_pc, R, Theta>>;
        DiscreteMappingCirc discrete_pseudo_cart_to_circ_to_cart(
                discrete_mapping_circ_to_cart,
                pseudo_cart_to_circ,
                1e-12);

        InvJacobianOPoint<PseudoCartToCircToCart, CoordRTheta> inv_o_point_analytical_circ(
                pseudo_cart_to_circ_to_cart);
        InvJacobianOPoint<DiscreteMappingCirc, CoordRTheta> inv_o_point_discrete_circ(
                discrete_pseudo_cart_to_circ_to_cart);
        analytical_matrix = inv_o_point_analytical_circ();
        discrete_matrix = inv_o_point_discrete_circ();
        double max_diff_circ
                = check_same(analytical_matrix, discrete_matrix, 1e-5 * ipow(16. / double(N), 4));
        std::cout << max_diff_circ << std::endl;



        // --- CZARNY MAPPING -----------------------------------------------------------------------------
        //std::cout << "   - Czarny mapping:   ";
        const CzarnyToCartesian<R, Theta, X, Y> czarny_to_cart(0.3, 1.4);
        using PseudoCartToCzarnyToCart = CombinedMapping<
                CzarnyToCartesian<R, Theta, X, Y>,
                CartesianToCircular<X_pc, Y_pc, R, Theta>>;
        const PseudoCartToCzarnyToCart
                pseudo_cart_to_czarny_to_cart(czarny_to_cart, pseudo_cart_to_circ, 1e-12);
        DiscreteToCartesianBuilder<X, Y, SplineRThetaBuilder_host, SplineRThetaEvaluator>
                mapping_builder_czarny(
                        Kokkos::DefaultHostExecutionSpace(),
                        czarny_to_cart,
                        builder,
                        spline_evaluator);
        DiscreteToCartesian discrete_mapping_czarny_to_cart = mapping_builder_czarny();
        using DiscreteMappingCzarny = CombinedMapping<
                decltype(discrete_mapping_czarny_to_cart),
                CartesianToCircular<X_pc, Y_pc, R, Theta>>;
        DiscreteMappingCzarny discrete_pseudo_cart_to_czarny_to_cart(
                discrete_mapping_czarny_to_cart,
                pseudo_cart_to_circ,
                1e-12);

        InvJacobianOPoint<PseudoCartToCzarnyToCart, CoordRTheta> inv_o_point_analytical_czarny(
                pseudo_cart_to_czarny_to_cart);
        InvJacobianOPoint<DiscreteMappingCzarny, CoordRTheta> inv_o_point_discrete_czarny(
                discrete_pseudo_cart_to_czarny_to_cart);
        analytical_matrix = inv_o_point_analytical_czarny();
        discrete_matrix = inv_o_point_discrete_czarny();
        double max_diff_czarny
                = check_same(analytical_matrix, discrete_matrix, 1e-5 * ipow(16. / double(N), 4));
        std::cout << max_diff_czarny << std::endl;


        std::array<double, 2> results;
        results[0] = max_diff_circ;
        results[1] = max_diff_czarny;
        return results;
    };

private:
    /**
     * @brief Check if two matrix have the same values.
     *
     * @param[in] matrix_1
     *      One matrix.
     * @param[in] matrix_2
     *      A second matrix.
     * @param[in] TOL
     *      Tolerance error.
     *
     * @return A double with the infinity norm of the difference.
     */
    double check_same(Matrix_2x2 const& matrix_1, Matrix_2x2 const& matrix_2, double TOL)
    {
        std::size_t size = 2;

        double max_diff = 0.;
        for (std::size_t i(0); i < size; ++i) {
            for (std::size_t j(0); j < size; ++j) {
                double const diff = fabs(matrix_1[i][j] - matrix_2[i][j]);
                max_diff = diff > max_diff ? diff : max_diff;

                EXPECT_NEAR(matrix_1[i][j], matrix_2[i][j], TOL);
            }
        }
        return max_diff;
    }
};


} // namespace



TEST(PseudoCartesianJacobianMatrix, TestDiscreteMapping)
{
    std::cout << "Comparison of the pseudo cartesian Jacobian matrix between analytical and "
              << std::endl
              << "discrete mappings at the center point: -----------------------------------"
              << std::endl
              << ">>> L_inf norm of | M_map - M_dis_map |" << std::endl;

    std::array<std::array<double, 2>, 4> results;
    results[0] = (PseudoCartesianJacobianMatrixTest<16>()).test_circ_and_czar();
    results[1] = (PseudoCartesianJacobianMatrixTest<32>()).test_circ_and_czar();
    results[2] = (PseudoCartesianJacobianMatrixTest<64>()).test_circ_and_czar();
    results[3] = (PseudoCartesianJacobianMatrixTest<128>()).test_circ_and_czar();



    std::cout << "Convergence order : ------------------------------------------------------"
              << std::endl
              << "Circular |  Czarny" << std::endl
              << "         |         " << std::endl;

    for (int i(0); i < results.size() - 1; i++) {
        double const order_circ = std::log(results[i][0] / results[i + 1][0]) / std::log(2);
        double const order_czar = std::log(results[i][1] / results[i + 1][1]) / std::log(2);
        std::cout << order_circ << "  |  " << order_czar << std::endl;
        int const BSDegree = 3;
        EXPECT_NEAR(order_circ, BSDegree + 1, 0.25);
        EXPECT_NEAR(order_czar, BSDegree + 1, 0.25);
    }
}
