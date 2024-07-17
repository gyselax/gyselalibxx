#pragma once

#include <ddc/ddc.hpp>
#include <ddc/kernels/splines.hpp>

#include <sll/mapping/curvilinear2d_to_cartesian.hpp>
#include <sll/mapping/discrete_mapping_to_cartesian.hpp>


/**
 * @brief A class for describing refined discrete 2D mappings from the logical domain to the physical domain.
 *
 * The mapping is a DiscreteToCartesian mapping but built on a different grid.
 * To be independent of errors due to the use of a discrete mapping, it could
 * be interesting to work with a discrete mapping built on a refined grid.
 *
 * Here, we wrap a DiscreteToCartesian mapping defined on a refined grid in a
 * RefinedDiscreteToCartesian mapping. The creation of the refined grid of size
 * Nr x Nt (@f$= N_r \times N_\theta @f$) is made in the class.
 *
 *
 * @see DiscreteToCartesian
 * @see Curvilinear2DToCartesian
 */
template <class RDimX, class RDimY, class SplineRPBuilder, class SplineRPEvaluator, int Nr, int Nt>
class RefinedDiscreteToCartesian
    : public Curvilinear2DToCartesian<
              RDimX,
              RDimY,
              typename SplineRPBuilder::continuous_dimension_type1,
              typename SplineRPBuilder::continuous_dimension_type2>
{
private:
    /**
     * @brief Indicate the bspline type of the first logical dimension.
     */
    using BSplineR = typename SplineRPBuilder::bsplines_type1;
    /**
     * @brief Indicate the bspline type of the second logical dimension.
     */
    using BSplineP = typename SplineRPBuilder::bsplines_type2;
    /**
     * @brief Indicate the first logical coordinate.
     */
    using RDimR = typename BSplineR::continuous_dimension_type;
    /**
     * @brief Indicate the second logical coordinate.
     */
    using RDimP = typename BSplineP::continuous_dimension_type;

    /**
     * @brief Indicate the first logical coordinate in the discrete space.
     */
    using IDimR = typename SplineRPBuilder::interpolation_discrete_dimension_type1;

    /**
     * @brief Indicate the second logical coordinate in the discrete space.
     */
    using IDimP = typename SplineRPBuilder::interpolation_discrete_dimension_type2;



    /**
     * @brief Boolean: true is BsplineR is defined on uniform sampling;
     * false if it is on non-uniform sampling.
     */
    static constexpr bool BSplineR_uniform = BSplineR::is_uniform();
    /**
     * @brief Boolean: true is BsplineP is defined on uniform sampling;
     * false if it is on non-uniform sampling.
     */
    static constexpr bool BSplineP_uniform = BSplineP::is_uniform();

    static int constexpr BSDegreeRRefined = BSplineR::degree();
    static int constexpr BSDegreePRefined = BSplineP::degree();


public:
    /**
     * @brief Define non periodic real refined R dimension.
     */
    struct RDimRRefined
    {
        /**
         * @brief Define periodicity of the dimension.
         * Here, not periodic.
         */
        static bool constexpr PERIODIC = false;
    };
    /**
     * @brief Define periodic real refined P dimension.
     */
    struct RDimPRefined
    {
        /**
         * @brief Define periodicity of the dimension.
         * Here, periodic.
         */
        static bool constexpr PERIODIC = true;
    };
    /**
     * @brief Define non periodic real refined X dimension.
     */
    struct RDimXRefined
    {
        /**
         * @brief Define periodicity of the dimension.
         * Here, not periodic.
         */
        static bool constexpr PERIODIC = false;
    };
    /**
     * @brief Define non periodic real refined X dimension.
     */
    struct RDimYRefined
    {
        /**
         * @brief Define periodicity of the dimension.
         * Here, not periodic.
         */
        static bool constexpr PERIODIC = false;
    };

public:
    struct BSplineRRefined
        : std::conditional_t<
                  BSplineR_uniform,
                  ddc::UniformBSplines<RDimRRefined, BSDegreeRRefined>,
                  ddc::NonUniformBSplines<RDimRRefined, BSDegreeRRefined>>
    {
    };
    struct BSplinePRefined
        : std::conditional_t<
                  BSplineP_uniform,
                  ddc::UniformBSplines<RDimPRefined, BSDegreePRefined>,
                  ddc::NonUniformBSplines<RDimPRefined, BSDegreePRefined>>
    {
    };

    static auto constexpr SplineRBoundaryRefined_min = SplineRPBuilder::builder_type1::s_bc_xmin;
    static auto constexpr SplineRBoundaryRefined_max = SplineRPBuilder::builder_type1::s_bc_xmax;
    static auto constexpr SplinePBoundaryRefined_min = SplineRPBuilder::builder_type2::s_bc_xmin;
    static auto constexpr SplinePBoundaryRefined_max = SplineRPBuilder::builder_type2::s_bc_xmax;


    static bool constexpr UniformMeshR = ddc::is_uniform_point_sampling_v<
            typename SplineRPBuilder::interpolation_discrete_dimension_type1>;
    static bool constexpr UniformMeshP = ddc::is_uniform_point_sampling_v<
            typename SplineRPBuilder::interpolation_discrete_dimension_type2>;


    struct IDimRRefined
        : std::conditional_t<
                  UniformMeshR,
                  ddc::UniformPointSampling<RDimRRefined>,
                  ddc::NonUniformPointSampling<RDimRRefined>>
    {
    };
    struct IDimPRefined
        : std::conditional_t<
                  UniformMeshP,
                  ddc::UniformPointSampling<RDimPRefined>,
                  ddc::NonUniformPointSampling<RDimPRefined>>
    {
    };


    using REvalBoundary = ddc::ConstantExtrapolationRule<RDimRRefined, RDimPRefined>;

    using SplineRPBuilderRefined = ddc::SplineBuilder2D<
            Kokkos::DefaultHostExecutionSpace,
            Kokkos::DefaultHostExecutionSpace::memory_space,
            BSplineRRefined,
            BSplinePRefined,
            IDimRRefined,
            IDimPRefined,
            SplineRBoundaryRefined_min,
            SplineRBoundaryRefined_max,
            SplinePBoundaryRefined_min,
            SplinePBoundaryRefined_max,
            ddc::SplineSolver::GINKGO,
            IDimRRefined,
            IDimPRefined>;

    using SplineRPEvaluatorRefined = ddc::SplineEvaluator2D<
            Kokkos::DefaultHostExecutionSpace,
            Kokkos::DefaultHostExecutionSpace::memory_space,
            BSplineRRefined,
            BSplinePRefined,
            IDimRRefined,
            IDimPRefined,
            REvalBoundary,
            REvalBoundary,
            ddc::PeriodicExtrapolationRule<RDimPRefined>,
            ddc::PeriodicExtrapolationRule<RDimPRefined>,
            IDimRRefined,
            IDimPRefined>;


    using CoordRRefined = ddc::Coordinate<RDimRRefined>;
    using CoordPRefined = ddc::Coordinate<RDimPRefined>;
    using CoordRPRefined = ddc::Coordinate<RDimRRefined, RDimPRefined>;

    using SplineInterpPointsRRefined = ddc::GrevilleInterpolationPoints<
            BSplineRRefined,
            SplineRBoundaryRefined_min,
            SplineRBoundaryRefined_max>;
    using SplineInterpPointsPRefined = ddc::GrevilleInterpolationPoints<
            BSplinePRefined,
            SplinePBoundaryRefined_min,
            SplinePBoundaryRefined_max>;

    using BSDomainRRefined = ddc::DiscreteDomain<BSplineRRefined>;
    using BSDomainPRefined = ddc::DiscreteDomain<BSplinePRefined>;
    using BSDomainRPRefined = ddc::DiscreteDomain<BSplineRRefined, BSplinePRefined>;

    using IndexRRefined = ddc::DiscreteElement<IDimRRefined>;
    using IndexPRefined = ddc::DiscreteElement<IDimPRefined>;
    using IndexRPRefined = ddc::DiscreteElement<IDimRRefined, IDimPRefined>;

    using IVectRRefined = ddc::DiscreteVector<IDimRRefined>;
    using IVectPRefined = ddc::DiscreteVector<IDimPRefined>;
    using IVectRPRefined = ddc::DiscreteVector<IDimRRefined, IDimPRefined>;

    using IDomainRRefined = ddc::DiscreteDomain<IDimRRefined>;
    using IDomainPRefined = ddc::DiscreteDomain<IDimPRefined>;
    using IDomainRPRefined = ddc::DiscreteDomain<IDimRRefined, IDimPRefined>;


    using spline_domain = ddc::DiscreteDomain<BSplineR, BSplineP>;

    /**
     * @brief Define a 2x2 matrix with an 2D array of an 2D array.
     */
    using Matrix_2x2 = std::array<std::array<double, 2>, 2>;

    IDomainRPRefined const m_refined_domain;
    REvalBoundary const boundary_condition_r_left;
    REvalBoundary const boundary_condition_r_right;
    SplineRPEvaluatorRefined const refined_evaluator;
    DiscreteToCartesian<
            RDimXRefined,
            RDimYRefined,
            SplineRPBuilderRefined,
            SplineRPEvaluatorRefined> const m_mapping;


    static inline ddc::Coordinate<RDimXRefined> to_refined(ddc::Coordinate<RDimX> const& coord)
    {
        return ddc::Coordinate<RDimXRefined>(ddc::get<RDimX>(coord));
    }

    static inline ddc::Coordinate<RDimYRefined> to_refined(ddc::Coordinate<RDimY> const& coord)
    {
        return ddc::Coordinate<RDimYRefined>(ddc::get<RDimY>(coord));
    }

    static inline ddc::Coordinate<RDimXRefined, RDimYRefined> to_refined(
            ddc::Coordinate<RDimX, RDimY> const& coord)
    {
        const double coord1 = ddc::get<RDimX>(coord);
        const double coord2 = ddc::get<RDimY>(coord);
        return ddc::Coordinate<RDimXRefined, RDimYRefined>(coord1, coord2);
    }


    static inline CoordRRefined to_refined(ddc::Coordinate<RDimR> const& coord)
    {
        return CoordRRefined(ddc::get<RDimR>(coord));
    }

    static inline CoordPRefined to_refined(ddc::Coordinate<RDimP> const& coord)
    {
        return CoordPRefined(ddc::get<RDimP>(coord));
    }

    static inline CoordRPRefined to_refined(ddc::Coordinate<RDimR, RDimP> const& coord)
    {
        const double coord1 = ddc::get<RDimR>(coord);
        const double coord2 = ddc::get<RDimP>(coord);
        return CoordRPRefined(coord1, coord2);
    }


    static inline ddc::Coordinate<RDimX> from_refined(ddc::Coordinate<RDimXRefined> const& coord)
    {
        return ddc::Coordinate<RDimX>(ddc::get<RDimXRefined>(coord));
    }

    static inline ddc::Coordinate<RDimY> from_refined(ddc::Coordinate<RDimYRefined> const& coord)
    {
        return ddc::Coordinate<RDimY>(ddc::get<RDimYRefined>(coord));
    }

    static inline ddc::Coordinate<RDimX, RDimY> from_refined(
            ddc::Coordinate<RDimXRefined, RDimYRefined> const& coord)
    {
        const double coord1 = ddc::get<RDimXRefined>(coord);
        const double coord2 = ddc::get<RDimYRefined>(coord);
        return ddc::Coordinate<RDimX, RDimY>(coord1, coord2);
    }


    static inline ddc::Coordinate<RDimR> from_refined(CoordRRefined const& coord)
    {
        return ddc::Coordinate<RDimR>(ddc::get<RDimRRefined>(coord));
    }

    static inline ddc::Coordinate<RDimP> from_refined(CoordPRefined const& coord)
    {
        return ddc::Coordinate<RDimP>(ddc::get<RDimPRefined>(coord));
    }

    static inline ddc::Coordinate<RDimR, RDimP> from_refined(CoordRPRefined const& coord)
    {
        const double coord1 = ddc::get<RDimRRefined>(coord);
        const double coord2 = ddc::get<RDimPRefined>(coord);
        return ddc::Coordinate<RDimR, RDimP>(coord1, coord2);
    }


    /**
     * @brief Instantiate a RefinedDiscreteToCartesian from the coefficients of 2D splines
     * approximating the mapping in the refined domain.
     *
     * This function is private. If the user has enough information to call this constructor
     * then they should create a DiscreteToCartesian directly.
     *
     * @param[in] curvilinear_to_x
     *      Bsplines coefficients of the first physical dimension in the refined logical domain.
     *
     * @param[in] curvilinear_to_y
     *      Bsplines coefficients of the second physical dimension in the refined logical domain.
     */
    RefinedDiscreteToCartesian(
            IDomainRPRefined const& refined_domain,
            ddc::Chunk<double, BSDomainRPRefined>&& curvilinear_to_x,
            ddc::Chunk<double, BSDomainRPRefined>&& curvilinear_to_y,
            CoordRRefined r_min,
            CoordRRefined r_max)
        : m_refined_domain(refined_domain)
        , boundary_condition_r_left(r_min)
        , boundary_condition_r_right(r_max)
        , refined_evaluator(
                  boundary_condition_r_left,
                  boundary_condition_r_right,
                  ddc::PeriodicExtrapolationRule<RDimPRefined>(),
                  ddc::PeriodicExtrapolationRule<RDimPRefined>())
        , m_mapping(std::move(curvilinear_to_x), std::move(curvilinear_to_y), refined_evaluator)
    {
    }


public:
    /**
     * @brief Instantiate a RefinedDiscreteToCartesian from another
     * RefinedDiscreteToCartesian.
     *
     * @param[in] x
     *      Another RefinedDiscreteToCartesian mapping used to
     *      instantiate the new one. (lvalue)
     */
    RefinedDiscreteToCartesian(RefinedDiscreteToCartesian&& x) = default;
    /**
     * @brief Instantiate a RefinedDiscreteToCartesian from another
     * temporary RefinedDiscreteToCartesian.
     *
     * @param[in] x
     *      Another temporary RefinedDiscreteToCartesian mapping used to
     *      instantiate the new one. (rvalue)
     */
    RefinedDiscreteToCartesian(RefinedDiscreteToCartesian const& x) = default;

    /**
     * @brief Compute the physical coordinates from the logical coordinates.
     *
     * It evaluates the decomposed mapping on B-splines at the coordinate point
     * with a SplineEvaluator2D.
     *
     * @param[in] coord
     *          The coordinates in the initial logical domain.
     *
     * @return The coordinates of the mapping in the initial physical domain.
     *
     * @see DiscreteToCartesian
     * @see SplineEvaluator2D
     */
    ddc::Coordinate<RDimX, RDimY> operator()(ddc::Coordinate<RDimR, RDimP> const& coord) const final
    {
        return from_refined(m_mapping(to_refined(coord)));
    }


    /**
     * @brief Compute full Jacobian matrix.
     *
     * For some computations, we need the complete Jacobian matrix or just the
     * coefficients.
     * The coefficients can be given indendently with the functions
     * jacobian_11, jacobian_12, jacobian_21 and jacobian_22.
     *
     * @param[in] coord
     *              The coordinate in the initial domain where we evaluate the Jacobian matrix.
     * @param[out] matrix
     *              The Jacobian matrix returned.
     *
     * @see DiscreteToCartesian
     * @see Curvilinear2DToCartesian::jacobian_11
     * @see Curvilinear2DToCartesian::jacobian_12
     * @see Curvilinear2DToCartesian::jacobian_21
     * @see Curvilinear2DToCartesian::jacobian_22
     */
    void jacobian_matrix(ddc::Coordinate<RDimR, RDimP> const& coord, Matrix_2x2& matrix) const final
    {
        m_mapping.jacobian_matrix(to_refined(coord), matrix);
    }

    /**
     * @brief Compute the (1,1) coefficient of the Jacobian matrix.
     *
     * @param[in] coord
     *              The coordinate in the initial domain where we evaluate the Jacobian matrix.
     *
     * @return A double with the value of the (1,1) coefficient of the Jacobian matrix.
     *
     * @see DiscreteToCartesian
     * @see SplineEvaluator2D
     */
    double jacobian_11(ddc::Coordinate<RDimR, RDimP> const& coord) const final
    {
        return m_mapping.jacobian_11(to_refined(coord));
    }

    /**
     * @brief Compute the (1,2) coefficient of the Jacobian matrix.
     *
     * @param[in] coord
     *              The coordinate in the initial domain where we evaluate the Jacobian matrix.
     *
     * @return A double with the value of the (1,2) coefficient of the Jacobian matrix.
     *
     * @see DiscreteToCartesian
     * @see SplineEvaluator2D
     */
    double jacobian_12(ddc::Coordinate<RDimR, RDimP> const& coord) const final
    {
        return m_mapping.jacobian_12(to_refined(coord));
    }

    /**
     * @brief Compute the (2,1) coefficient of the Jacobian matrix.
     *
     * @param[in] coord
     *              The coordinate in the initial domain where we evaluate the Jacobian matrix.
     *
     * @return A double with the value of the (2,1) coefficient of the Jacobian matrix.
     *
     * @see DiscreteToCartesian
     * @see SplineEvaluator2D
     */
    double jacobian_21(ddc::Coordinate<RDimR, RDimP> const& coord) const final
    {
        return m_mapping.jacobian_21(to_refined(coord));
    }

    /**
     * @brief Compute the (2,2) coefficient of the Jacobian matrix.
     *
     * @param[in] coord
     *              The coordinate in the initial domain where we evaluate the Jacobian matrix.
     *
     * @return A double with the value of the (2,2) coefficient of the Jacobian matrix.
     *
     * @see DiscreteToCartesian
     * @see SplineEvaluator2D
     */
    double jacobian_22(ddc::Coordinate<RDimR, RDimP> const& coord) const final
    {
        return m_mapping.jacobian_22(to_refined(coord));
    }


    /**
     * @brief  Compute the full Jacobian matrix from the mapping to the pseudo-Cartesian mapping at the central point.
     *
     * Call the DiscreteToCartesian::to_pseudo_cartesian_jacobian_center_matrix() function.
     *
     * @param[in] grid
     *      The domain where the mapping is defined.
     * @param[out] matrix
     *      The pseudo-Cartesian matrix evaluated at the central point.
     *
     *
     * @see Discrete2DToCartesian
     * @see BslAdvection
     * @see AdvectionDomain
     */
    void to_pseudo_cartesian_jacobian_center_matrix(
            ddc::DiscreteDomain<IDimR, IDimP> const& grid,
            Matrix_2x2& matrix) const
    {
        m_mapping.to_pseudo_cartesian_jacobian_center_matrix(m_refined_domain, matrix);
    }


    /**
     * @brief Compute the (1,1) coefficient of the pseudo-Cartesian Jacobian matrix at the central point.
     *
     * @param[in] grid
     *      The domain where the mapping is defined.
     *
     * @return A double with the (1,1) coefficient of the pseudo-Cartesian Jacobian matrix at the central point.
     *
     * @see to_pseudo_cartesian_jacobian_center_matrix
     */
    double to_pseudo_cartesian_jacobian_11_center(
            ddc::DiscreteDomain<IDimR, IDimP> const& grid) const
    {
        Matrix_2x2 jacobian;
        to_pseudo_cartesian_jacobian_center_matrix(grid, jacobian);
        return jacobian[0][0];
    }


    /**
     * @brief Compute the (1,2) coefficient of the pseudo-Cartesian Jacobian matrix at the central point.
     *
     * @param[in] grid
     *      The domain where the mapping is defined.
     *
     * @return A double with the (1,2) coefficient of the pseudo-Cartesian Jacobian matrix at the central point.
     *
     * @see Curvilinear2DToCartesian::to_pseudo_cartesian_jacobian_center_matrix
     */
    double to_pseudo_cartesian_jacobian_12_center(
            ddc::DiscreteDomain<IDimR, IDimP> const& grid) const
    {
        Matrix_2x2 jacobian;
        to_pseudo_cartesian_jacobian_center_matrix(grid, jacobian);
        return jacobian[0][1];
    }


    /**
     * @brief Compute the (2,1) coefficient of the pseudo-Cartesian Jacobian matrix at the central point.
     *
     * @param[in] grid
     *      The domain where the mapping is defined.
     *
     * @return A double with the (2,1) coefficient of the pseudo-Cartesian Jacobian matrix at the central point.
     *
     * @see to_pseudo_cartesian_jacobian_center_matrix
     */
    double to_pseudo_cartesian_jacobian_21_center(
            ddc::DiscreteDomain<IDimR, IDimP> const& grid) const
    {
        Matrix_2x2 jacobian;
        to_pseudo_cartesian_jacobian_center_matrix(grid, jacobian);
        return jacobian[1][0];
    }


    /**
     * @brief Compute the (2,2) coefficient of the pseudo-Cartesian Jacobian matrix at the central point.
     *
     * @param[in] grid
     *      The domain where the mapping is defined.
     *
     * @return A double with the (2,2) coefficient of the pseudo-Cartesian Jacobian matrix at the central point.
     *
     * @see to_pseudo_cartesian_jacobian_center_matrix
     */
    double to_pseudo_cartesian_jacobian_22_center(
            ddc::DiscreteDomain<IDimR, IDimP> const& grid) const
    {
        ddc::DiscreteDomain<IDimP> const theta_domain = ddc::select<IDimP>(grid);

        Matrix_2x2 jacobian;
        to_pseudo_cartesian_jacobian_center_matrix(grid, jacobian);
        return jacobian[1][1];
    }



    /**
     * @brief Define a RefinedDiscreteToCartesian mapping from an analytical mapping.
     *
     * Create a refined grid of size @f$ N_r \times N_\theta @f$ (with  @f$N_r @f$ and  @f$N_\theta @f$
     * the templated parameters of the class). Define new operators on the refined grid and decompose
     * the analytical mapping on the B-splines basis.
     *
     * @param[in] analytical_mapping
     *          The mapping defined analytically.
     * @param[in] domain
     *          The initial domain.
     *
     * @return A RefinedDiscreteToCartesian version of the analytical mapping.
     *
     * @see DiscreteToCartesian
     */
    template <class Mapping>
    static RefinedDiscreteToCartesian analytical_to_refined(
            Mapping const& analytical_mapping,
            ddc::DiscreteDomain<IDimR, IDimP> const& domain)
    {
        const double rmin = ddc::discrete_space<BSplineR>().rmin();
        const double rmax = ddc::discrete_space<BSplineR>().rmax();

        const double pmin = ddc::discrete_space<BSplineP>().rmin();
        const double pmax = ddc::discrete_space<BSplineP>().rmax();

        // Create refined grid
        CoordRRefined const r_min(rmin);
        CoordRRefined const r_max(rmax);
        IVectRRefined const r_size(Nr);

        CoordPRefined const p_min(pmin);
        CoordPRefined const p_max(pmax);
        IVectPRefined const p_size(Nt);

        if constexpr (BSplineRRefined::is_uniform()) {
            ddc::init_discrete_space<BSplineRRefined>(r_min, r_max, r_size);
            ddc::init_discrete_space<BSplinePRefined>(p_min, p_max, p_size);
        } else {
            double const dr((r_max - r_min) / r_size);
            double const dp((p_max - p_min) / p_size);

            std::vector<CoordRRefined> r_knots(r_size + 1);
            std::vector<CoordPRefined> p_knots(p_size + 1);

            for (int i(0); i < r_size + 1; ++i) {
                r_knots[i] = CoordRRefined(r_min + i * dr);
            }
            for (int i(0); i < p_size + 1; ++i) {
                p_knots[i] = CoordPRefined(p_min + i * dp);
            }

            ddc::init_discrete_space<BSplineRRefined>(r_knots);
            ddc::init_discrete_space<BSplinePRefined>(p_knots);
        }

        ddc::init_discrete_space<IDimRRefined>(
                SplineInterpPointsRRefined::template get_sampling<IDimRRefined>());
        ddc::init_discrete_space<IDimPRefined>(
                SplineInterpPointsPRefined::template get_sampling<IDimPRefined>());

        IDomainRRefined const interpolation_domain_R(
                SplineInterpPointsRRefined::template get_domain<IDimRRefined>());
        IDomainPRefined const interpolation_domain_P(
                SplineInterpPointsPRefined::template get_domain<IDimPRefined>());
        IDomainRPRefined const refined_domain(interpolation_domain_R, interpolation_domain_P);

        // Operators on the refined grid
        SplineRPBuilderRefined const refined_builder(refined_domain);

        BSDomainRPRefined const spline_domain = refined_builder.spline_domain();

        // Compute the B-splines coefficients of the analytical mapping
        ddc::Chunk<double, BSDomainRPRefined> curvilinear_to_x_spline(spline_domain);
        ddc::Chunk<double, BSDomainRPRefined> curvilinear_to_y_spline(spline_domain);
        ddc::Chunk<double, IDomainRPRefined> curvilinear_to_x_vals(refined_domain);
        ddc::Chunk<double, IDomainRPRefined> curvilinear_to_y_vals(refined_domain);
        ddc::for_each(
                refined_builder.interpolation_domain(),
                [&](typename ddc::DiscreteDomain<IDimRRefined, IDimPRefined>::
                            discrete_element_type const& el) {
                    ddc::Coordinate<RDimR, RDimP> polar_coord(from_refined(ddc::coordinate(el)));
                    ddc::Coordinate<RDimXRefined, RDimYRefined> cart_coord
                            = to_refined(analytical_mapping(polar_coord));

                    curvilinear_to_x_vals(el) = ddc::select<RDimXRefined>(cart_coord);
                    curvilinear_to_y_vals(el) = ddc::select<RDimYRefined>(cart_coord);
                });
        refined_builder(curvilinear_to_x_spline.span_view(), curvilinear_to_x_vals.span_cview());
        refined_builder(curvilinear_to_y_spline.span_view(), curvilinear_to_y_vals.span_cview());

        return RefinedDiscreteToCartesian(
                refined_domain,
                std::move(curvilinear_to_x_spline),
                std::move(curvilinear_to_y_spline),
                r_min,
                r_max);
    }
};
