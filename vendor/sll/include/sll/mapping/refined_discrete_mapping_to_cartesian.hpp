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
template <class X, class Y, class SplineRThetaBuilder, class SplineRThetaEvaluator, int Nr, int Nt>
class RefinedDiscreteToCartesian
    : public Curvilinear2DToCartesian<
              X,
              Y,
              typename SplineRThetaBuilder::continuous_dimension_type1,
              typename SplineRThetaBuilder::continuous_dimension_type2>
{
private:
    /**
     * @brief Indicate the bspline type of the first logical dimension.
     */
    using BSplineR = typename SplineRThetaBuilder::bsplines_type1;
    /**
     * @brief Indicate the bspline type of the second logical dimension.
     */
    using BSplineTheta = typename SplineRThetaBuilder::bsplines_type2;
    /**
     * @brief Indicate the first logical coordinate.
     */
    using R = typename BSplineR::continuous_dimension_type;
    /**
     * @brief Indicate the second logical coordinate.
     */
    using Theta = typename BSplineTheta::continuous_dimension_type;

    /**
     * @brief Indicate the first logical coordinate in the discrete space.
     */
    using GridR = typename SplineRThetaBuilder::interpolation_discrete_dimension_type1;

    /**
     * @brief Indicate the second logical coordinate in the discrete space.
     */
    using GridTheta = typename SplineRThetaBuilder::interpolation_discrete_dimension_type2;



    /**
     * @brief Boolean: true is BsplineR is defined on uniform sampling;
     * false if it is on non-uniform sampling.
     */
    static constexpr bool BSplineR_uniform = BSplineR::is_uniform();
    /**
     * @brief Boolean: true is BsplineP is defined on uniform sampling;
     * false if it is on non-uniform sampling.
     */
    static constexpr bool BSplineP_uniform = BSplineTheta::is_uniform();

    static int constexpr BSDegreeRRefined = BSplineR::degree();
    static int constexpr BSDegreePRefined = BSplineTheta::degree();


public:
    /**
     * @brief Define non periodic real refined R dimension.
     */
    struct RRefined
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
    struct ThetaRefined
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
    struct XRefined
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
    struct YRefined
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
                  ddc::UniformBSplines<RRefined, BSDegreeRRefined>,
                  ddc::NonUniformBSplines<RRefined, BSDegreeRRefined>>
    {
    };
    struct BSplineThetaRefined
        : std::conditional_t<
                  BSplineP_uniform,
                  ddc::UniformBSplines<ThetaRefined, BSDegreePRefined>,
                  ddc::NonUniformBSplines<ThetaRefined, BSDegreePRefined>>
    {
    };

    static auto constexpr SplineRBoundaryRefined_min
            = SplineRThetaBuilder::builder_type1::s_bc_xmin;
    static auto constexpr SplineRBoundaryRefined_max
            = SplineRThetaBuilder::builder_type1::s_bc_xmax;
    static auto constexpr SplinePBoundaryRefined_min
            = SplineRThetaBuilder::builder_type2::s_bc_xmin;
    static auto constexpr SplinePBoundaryRefined_max
            = SplineRThetaBuilder::builder_type2::s_bc_xmax;


    static bool constexpr UniformMeshR = ddc::is_uniform_point_sampling_v<
            typename SplineRThetaBuilder::interpolation_discrete_dimension_type1>;
    static bool constexpr UniformMeshP = ddc::is_uniform_point_sampling_v<
            typename SplineRThetaBuilder::interpolation_discrete_dimension_type2>;


    struct GridRRefined
        : std::conditional_t<
                  UniformMeshR,
                  ddc::UniformPointSampling<RRefined>,
                  ddc::NonUniformPointSampling<RRefined>>
    {
    };
    struct GridThetaRefined
        : std::conditional_t<
                  UniformMeshP,
                  ddc::UniformPointSampling<ThetaRefined>,
                  ddc::NonUniformPointSampling<ThetaRefined>>
    {
    };


    using REvalBoundary = ddc::ConstantExtrapolationRule<RRefined, ThetaRefined>;

    using SplineRThetaBuilderRefined = ddc::SplineBuilder2D<
            Kokkos::DefaultHostExecutionSpace,
            Kokkos::DefaultHostExecutionSpace::memory_space,
            BSplineRRefined,
            BSplineThetaRefined,
            GridRRefined,
            GridThetaRefined,
            SplineRBoundaryRefined_min,
            SplineRBoundaryRefined_max,
            SplinePBoundaryRefined_min,
            SplinePBoundaryRefined_max,
            ddc::SplineSolver::LAPACK,
            GridRRefined,
            GridThetaRefined>;

    using SplineRThetaEvaluatorRefined = ddc::SplineEvaluator2D<
            Kokkos::DefaultHostExecutionSpace,
            Kokkos::DefaultHostExecutionSpace::memory_space,
            BSplineRRefined,
            BSplineThetaRefined,
            GridRRefined,
            GridThetaRefined,
            REvalBoundary,
            REvalBoundary,
            ddc::PeriodicExtrapolationRule<ThetaRefined>,
            ddc::PeriodicExtrapolationRule<ThetaRefined>,
            GridRRefined,
            GridThetaRefined>;


    using CoordRRefined = ddc::Coordinate<RRefined>;
    using CoordThetaRefined = ddc::Coordinate<ThetaRefined>;
    using CoordRThetaRefined = ddc::Coordinate<RRefined, ThetaRefined>;

    using SplineInterpThetaointsRRefined = ddc::GrevilleInterpolationPoints<
            BSplineRRefined,
            SplineRBoundaryRefined_min,
            SplineRBoundaryRefined_max>;
    using SplineInterpThetaointsThetaRefined = ddc::GrevilleInterpolationPoints<
            BSplineThetaRefined,
            SplinePBoundaryRefined_min,
            SplinePBoundaryRefined_max>;

    using BSIdxRangeRRefined = ddc::DiscreteDomain<BSplineRRefined>;
    using BSIdxRangeThetaRefined = ddc::DiscreteDomain<BSplineThetaRefined>;
    using BSIdxRangeRThetaRefined = ddc::DiscreteDomain<BSplineRRefined, BSplineThetaRefined>;

    using IdxRRefined = ddc::DiscreteElement<GridRRefined>;
    using IdxThetaRefined = ddc::DiscreteElement<GridThetaRefined>;
    using IdxRThetaRefined = ddc::DiscreteElement<GridRRefined, GridThetaRefined>;

    using IdxStepRRefined = ddc::DiscreteVector<GridRRefined>;
    using IdxStepThetaRefined = ddc::DiscreteVector<GridThetaRefined>;
    using IdxStepRThetaRefined = ddc::DiscreteVector<GridRRefined, GridThetaRefined>;

    using IdxRangeRRefined = ddc::DiscreteDomain<GridRRefined>;
    using IdxRangeThetaRefined = ddc::DiscreteDomain<GridThetaRefined>;
    using IdxRangeRThetaRefined = ddc::DiscreteDomain<GridRRefined, GridThetaRefined>;


    using spline_domain = ddc::DiscreteDomain<BSplineR, BSplineTheta>;

    /**
     * @brief Define a 2x2 matrix with an 2D array of an 2D array.
     */
    using Matrix_2x2 = std::array<std::array<double, 2>, 2>;

    IdxRangeRThetaRefined const m_refined_domain;
    REvalBoundary const boundary_condition_r_left;
    REvalBoundary const boundary_condition_r_right;
    SplineRThetaEvaluatorRefined const refined_evaluator;
    DiscreteToCartesian<
            XRefined,
            YRefined,
            SplineRThetaBuilderRefined,
            SplineRThetaEvaluatorRefined> const m_mapping;


    static inline ddc::Coordinate<XRefined> to_refined(ddc::Coordinate<X> const& coord)
    {
        return ddc::Coordinate<XRefined>(ddc::get<X>(coord));
    }

    static inline ddc::Coordinate<YRefined> to_refined(ddc::Coordinate<Y> const& coord)
    {
        return ddc::Coordinate<YRefined>(ddc::get<Y>(coord));
    }

    static inline ddc::Coordinate<XRefined, YRefined> to_refined(ddc::Coordinate<X, Y> const& coord)
    {
        const double coord1 = ddc::get<X>(coord);
        const double coord2 = ddc::get<Y>(coord);
        return ddc::Coordinate<XRefined, YRefined>(coord1, coord2);
    }


    static inline CoordRRefined to_refined(ddc::Coordinate<R> const& coord)
    {
        return CoordRRefined(ddc::get<R>(coord));
    }

    static inline CoordThetaRefined to_refined(ddc::Coordinate<Theta> const& coord)
    {
        return CoordThetaRefined(ddc::get<Theta>(coord));
    }

    static inline CoordRThetaRefined to_refined(ddc::Coordinate<R, Theta> const& coord)
    {
        const double coord1 = ddc::get<R>(coord);
        const double coord2 = ddc::get<Theta>(coord);
        return CoordRThetaRefined(coord1, coord2);
    }


    static inline ddc::Coordinate<X> from_refined(ddc::Coordinate<XRefined> const& coord)
    {
        return ddc::Coordinate<X>(ddc::get<XRefined>(coord));
    }

    static inline ddc::Coordinate<Y> from_refined(ddc::Coordinate<YRefined> const& coord)
    {
        return ddc::Coordinate<Y>(ddc::get<YRefined>(coord));
    }

    static inline ddc::Coordinate<X, Y> from_refined(
            ddc::Coordinate<XRefined, YRefined> const& coord)
    {
        const double coord1 = ddc::get<XRefined>(coord);
        const double coord2 = ddc::get<YRefined>(coord);
        return ddc::Coordinate<X, Y>(coord1, coord2);
    }


    static inline ddc::Coordinate<R> from_refined(CoordRRefined const& coord)
    {
        return ddc::Coordinate<R>(ddc::get<RRefined>(coord));
    }

    static inline ddc::Coordinate<Theta> from_refined(CoordThetaRefined const& coord)
    {
        return ddc::Coordinate<Theta>(ddc::get<ThetaRefined>(coord));
    }

    static inline ddc::Coordinate<R, Theta> from_refined(CoordRThetaRefined const& coord)
    {
        const double coord1 = ddc::get<RRefined>(coord);
        const double coord2 = ddc::get<ThetaRefined>(coord);
        return ddc::Coordinate<R, Theta>(coord1, coord2);
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
            IdxRangeRThetaRefined const& refined_domain,
            ddc::Chunk<double, BSIdxRangeRThetaRefined>&& curvilinear_to_x,
            ddc::Chunk<double, BSIdxRangeRThetaRefined>&& curvilinear_to_y,
            CoordRRefined r_min,
            CoordRRefined r_max)
        : m_refined_domain(refined_domain)
        , boundary_condition_r_left(r_min)
        , boundary_condition_r_right(r_max)
        , refined_evaluator(
                  boundary_condition_r_left,
                  boundary_condition_r_right,
                  ddc::PeriodicExtrapolationRule<ThetaRefined>(),
                  ddc::PeriodicExtrapolationRule<ThetaRefined>())
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
    ddc::Coordinate<X, Y> operator()(ddc::Coordinate<R, Theta> const& coord) const final
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
    void jacobian_matrix(ddc::Coordinate<R, Theta> const& coord, Matrix_2x2& matrix) const final
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
    double jacobian_11(ddc::Coordinate<R, Theta> const& coord) const final
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
    double jacobian_12(ddc::Coordinate<R, Theta> const& coord) const final
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
    double jacobian_21(ddc::Coordinate<R, Theta> const& coord) const final
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
    double jacobian_22(ddc::Coordinate<R, Theta> const& coord) const final
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
            ddc::DiscreteDomain<GridR, GridTheta> const& grid,
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
            ddc::DiscreteDomain<GridR, GridTheta> const& grid) const
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
            ddc::DiscreteDomain<GridR, GridTheta> const& grid) const
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
            ddc::DiscreteDomain<GridR, GridTheta> const& grid) const
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
            ddc::DiscreteDomain<GridR, GridTheta> const& grid) const
    {
        ddc::DiscreteDomain<GridTheta> const theta_domain = ddc::select<GridTheta>(grid);

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
            ddc::DiscreteDomain<GridR, GridTheta> const& domain)
    {
        const double rmin = ddc::discrete_space<BSplineR>().rmin();
        const double rmax = ddc::discrete_space<BSplineR>().rmax();

        const double pmin = ddc::discrete_space<BSplineTheta>().rmin();
        const double pmax = ddc::discrete_space<BSplineTheta>().rmax();

        // Create refined grid
        CoordRRefined const r_min(rmin);
        CoordRRefined const r_max(rmax);
        IdxStepRRefined const r_size(Nr);

        CoordThetaRefined const theta_min(pmin);
        CoordThetaRefined const theta_max(pmax);
        IdxStepThetaRefined const theta_size(Nt);

        if constexpr (BSplineRRefined::is_uniform()) {
            ddc::init_discrete_space<BSplineRRefined>(r_min, r_max, r_size);
            ddc::init_discrete_space<BSplineThetaRefined>(theta_min, theta_max, theta_size);
        } else {
            double const dr((r_max - r_min) / r_size);
            double const dp((theta_max - theta_min) / theta_size);

            std::vector<CoordRRefined> r_knots(r_size + 1);
            std::vector<CoordThetaRefined> theta_knots(theta_size + 1);

            for (int i(0); i < r_size + 1; ++i) {
                r_knots[i] = CoordRRefined(r_min + i * dr);
            }
            for (int i(0); i < theta_size + 1; ++i) {
                theta_knots[i] = CoordThetaRefined(theta_min + i * dp);
            }

            ddc::init_discrete_space<BSplineRRefined>(r_knots);
            ddc::init_discrete_space<BSplineThetaRefined>(theta_knots);
        }

        ddc::init_discrete_space<GridRRefined>(
                SplineInterpThetaointsRRefined::template get_sampling<GridRRefined>());
        ddc::init_discrete_space<GridThetaRefined>(
                SplineInterpThetaointsThetaRefined::template get_sampling<GridThetaRefined>());

        IdxRangeRRefined const interpolation_domain_R(
                SplineInterpThetaointsRRefined::template get_domain<GridRRefined>());
        IdxRangeThetaRefined const interpolation_domain_P(
                SplineInterpThetaointsThetaRefined::template get_domain<GridThetaRefined>());
        IdxRangeRThetaRefined const refined_domain(interpolation_domain_R, interpolation_domain_P);

        // Operators on the refined grid
        SplineRThetaBuilderRefined const refined_builder(refined_domain);

        BSIdxRangeRThetaRefined const spline_domain = refined_builder.spline_domain();

        // Compute the B-splines coefficients of the analytical mapping
        ddc::Chunk<double, BSIdxRangeRThetaRefined> curvilinear_to_x_spline(spline_domain);
        ddc::Chunk<double, BSIdxRangeRThetaRefined> curvilinear_to_y_spline(spline_domain);
        ddc::Chunk<double, IdxRangeRThetaRefined> curvilinear_to_x_vals(refined_domain);
        ddc::Chunk<double, IdxRangeRThetaRefined> curvilinear_to_y_vals(refined_domain);
        ddc::for_each(
                refined_builder.interpolation_domain(),
                [&](typename ddc::DiscreteDomain<GridRRefined, GridThetaRefined>::
                            discrete_element_type const& el) {
                    ddc::Coordinate<R, Theta> polar_coord(from_refined(ddc::coordinate(el)));
                    ddc::Coordinate<XRefined, YRefined> cart_coord
                            = to_refined(analytical_mapping(polar_coord));

                    curvilinear_to_x_vals(el) = ddc::select<XRefined>(cart_coord);
                    curvilinear_to_y_vals(el) = ddc::select<YRefined>(cart_coord);
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
