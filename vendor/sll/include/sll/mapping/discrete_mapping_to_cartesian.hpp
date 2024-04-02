#pragma once

#include <ddc/ddc.hpp>
#include <ddc/kernels/splines.hpp>

#include <sll/mapping/curvilinear2d_to_cartesian.hpp>


/**
 * @brief A class for describing discrete 2D mappings from the logical domain to the physical domain.
 *
 * The mapping describe here is only defined on a grid. The DiscreteToCartesian class decomposes the mapping
 * on B-splines to evaluate it on the physical domain.
 *
 * @f$ x(r,\theta) = \sum_k c_{x,k} B_k(r,\theta),@f$
 *
 * @f$ y(r,\theta) = \sum_k c_{y,k} B_k(r,\theta).@f$
 *
 * This mapping could be costly to inverse.
 *
 * @see Curvilinear2DToCartesian
 */
template <class DimX, class DimY, class SplineBuilder, class SplineEvaluator>
class DiscreteToCartesian
    : public Curvilinear2DToCartesian<
              DimX,
              DimY,
              typename SplineBuilder::bsplines_type1::tag_type,
              typename SplineBuilder::bsplines_type2::tag_type>
{
public:
    /**
     * @brief Indicate the bspline type of the first logical dimension.
     */
    using BSplineR = typename SplineBuilder::bsplines_type1;
    /**
     * @brief Indicate the bspline type of the second logical dimension.
     */
    using BSplineP = typename SplineBuilder::bsplines_type2;

    /**
     * @brief Indicate the first physical coordinate.
     */
    using cartesian_tag_x = DimX;
    /**
     * @brief Indicate the second physical coordinate.
     */
    using cartesian_tag_y = DimY;
    /**
     * @brief Indicate the first logical coordinate.
     */
    using circular_tag_r = typename BSplineR::tag_type;
    /**
     * @brief Indicate the second logical coordinate.
     */
    using circular_tag_p = typename BSplineP::tag_type;

    /**
     * @brief Indicate the first logical coordinate in the discrete space.
     */
    using IDimR = typename SplineBuilder::interpolation_mesh_type1;

    /**
     * @brief Indicate the second logical coordinate in the discrete space.
     */
    using IDimP = typename SplineBuilder::interpolation_mesh_type2;

    /**
     * @brief Define a 2x2 matrix with an 2D array of an 2D array.
     */
    using Matrix_2x2 = std::array<std::array<double, 2>, 2>;

private:
    using interpolation_domain = typename SplineBuilder::interpolation_domain_type;
    using spline_domain = ddc::DiscreteDomain<BSplineR, BSplineP>;

    using SplineType = ddc::Chunk<
            double,
            spline_domain,
            ddc::KokkosAllocator<double, typename SplineBuilder::memory_space>>;

private:
    SplineType x_spline_representation;
    SplineType y_spline_representation;
    SplineEvaluator const& spline_evaluator;

public:
    /**
     * @brief Instantiate a DiscreteToCartesian from the coefficients of 2D splines approximating the mapping.
     *
     * A discrete mapping is a mapping whose values are known only at the mesh points of the grid.
     * To interpolate the mapping, we use B-splines. The DiscreteToCartesian mapping is initialized
     * from the coefficients in front of the basis splines which arise when we approximate the
     * functions @f$ x(r,\theta) @f$, and @f$ y(r,\theta) @f$ (with @f$ x @f$ and @f$ y @f$ the physical dimensions in
     * the logical domain) with Splines (using SplineBuilder2D). Then to interpolate the mapping,
     * we will evaluate the decomposed functions on B-splines (see DiscreteToCartesian::operator()).

     *
     * Here, the evaluator is given as input.
     *
     * @param[in] curvilinear_to_x
     * 		Bsplines coefficients of the first physical dimension in the logical domain.
     *
     * @param[in] curvilinear_to_y
     * 		Bsplines coefficients of the second physical dimension in the logical domain.
     *
     * @param[in] evaluator
     * 		The evaluator used to evaluate the mapping.
     *
     *
     * @see SplineBuilder2D
     * @see DiscreteToCartesian::operator()
     * @see SplineBoundaryValue
     */
    DiscreteToCartesian(
            SplineType&& curvilinear_to_x,
            SplineType&& curvilinear_to_y,
            SplineEvaluator const& evaluator)
        : x_spline_representation(std::move(curvilinear_to_x))
        , y_spline_representation(std::move(curvilinear_to_y))
        , spline_evaluator(evaluator)
    {
    }

    /**
     * @brief Compute the physical coordinates from the logical coordinates.
     *
     * It evaluates the decomposed mapping on B-splines at the coordinate point
     * with a SplineEvaluator2D.
     *
     * @param[in] coord
     * 			The coordinates in the logical domain.
     *
     * @return The coordinates of the mapping in the physical domain.
     *
     * @see SplineEvaluator2D
     */
    ddc::Coordinate<DimX, DimY> operator()(
            ddc::Coordinate<circular_tag_r, circular_tag_p> const& coord) const final
    {
        const double x = spline_evaluator(coord, x_spline_representation.span_cview());
        const double y = spline_evaluator(coord, y_spline_representation.span_cview());
        return ddc::Coordinate<DimX, DimY>(x, y);
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
     * 				The coordinate where we evaluate the Jacobian matrix.
     * @param[out] matrix
     * 				The Jacobian matrix returned.
     *
     * @see Curvilinear2DToCartesian::jacobian_11
     * @see Curvilinear2DToCartesian::jacobian_12
     * @see Curvilinear2DToCartesian::jacobian_21
     * @see Curvilinear2DToCartesian::jacobian_22
     */
    void jacobian_matrix(
            ddc::Coordinate<circular_tag_r, circular_tag_p> const& coord,
            Matrix_2x2& matrix) const final
    {
        matrix[0][0] = spline_evaluator.deriv_dim_1(coord, x_spline_representation.span_cview());
        matrix[0][1] = spline_evaluator.deriv_dim_2(coord, x_spline_representation.span_cview());
        matrix[1][0] = spline_evaluator.deriv_dim_1(coord, y_spline_representation.span_cview());
        matrix[1][1] = spline_evaluator.deriv_dim_2(coord, y_spline_representation.span_cview());
    }

    /**
     * @brief Compute the (1,1) coefficient of the Jacobian matrix.
     *
     * For a mapping given by @f$ \mathcal{F} : (r,\theta)\mapsto (x,y) @f$, the
     * (1,1) coefficient of the Jacobian matrix is given by @f$ \frac{\partial x}{\partial r} @f$.
     * As the mapping is decomposed on B-splines, it means it computes the derivatives of B-splines
     * @f$ \frac{\partial x}{\partial r} (r,\theta)= \sum_k c_{x,k} \frac{\partial B_k}{\partial r}(r,\theta)@f$
     * (the derivatives are implemented in SplineEvaluator2D).
     *
     * @param[in] coord
     * 				The coordinate where we evaluate the Jacobian matrix.
     *
     * @return A double with the value of the (1,1) coefficient of the Jacobian matrix.
     *
     * @see SplineEvaluator2D
     */
    double jacobian_11(ddc::Coordinate<circular_tag_r, circular_tag_p> const& coord) const final
    {
        return spline_evaluator.deriv_dim_1(coord, x_spline_representation.span_cview());
    }

    /**
     * @brief Compute the (1,2) coefficient of the Jacobian matrix.
     *
     * For a mapping given by @f$ \mathcal{F} : (r,\theta)\mapsto (x,y) @f$, the
     * (1,2) coefficient of the Jacobian matrix is given by @f$ \frac{\partial x}{\partial \theta} @f$.
     * As the mapping is decomposed on B-splines, it means it computes
     * @f$ \frac{\partial x}{\partial \theta}(r,\theta) = \sum_k c_{x,k} \frac{\partial B_k}{\partial \theta}(r,\theta) @f$
     * (the derivatives of B-splines are implemented in SplineEvaluator2D).
     *
     * @param[in] coord
     * 				The coordinate where we evaluate the Jacobian matrix.
     *
     * @return A double with the value of the (1,2) coefficient of the Jacobian matrix.
     *
     * @see SplineEvaluator2D
     */
    double jacobian_12(ddc::Coordinate<circular_tag_r, circular_tag_p> const& coord) const final
    {
        return spline_evaluator.deriv_dim_2(coord, x_spline_representation.span_cview());
    }

    /**
     * @brief Compute the (2,1) coefficient of the Jacobian matrix.
     *
     *For a mapping given by @f$ \mathcal{F} : (r,\theta)\mapsto (x,y) @f$, the
     * (2,1) coefficient of the Jacobian matrix is given by @f$ \frac{\partial y}{\partial r} @f$.
     * As the mapping is decomposed on B-splines, it means it computes
     * @f$ \frac{\partial y}{\partial r}(r,\theta) = \sum_k c_{y,k} \frac{\partial B_k}{\partial r}(r,\theta)@f$
     * (the derivatives of B-splines are implemented in SplineEvaluator2D).
     *
     * @param[in] coord
     * 				The coordinate where we evaluate the Jacobian matrix. .
     *
     * @return A double with the value of the (2,1) coefficient of the Jacobian matrix.
     *
     * @see SplineEvaluator2D
     */
    double jacobian_21(ddc::Coordinate<circular_tag_r, circular_tag_p> const& coord) const final
    {
        return spline_evaluator.deriv_dim_1(coord, y_spline_representation.span_cview());
    }

    /**
     * @brief Compute the (2,2) coefficient of the Jacobian matrix.
     *
     *For a mapping given by @f$ \mathcal{F} : (r,\theta)\mapsto (x,y) @f$, the
     * (2,2) coefficient of the Jacobian matrix is given by @f$ \frac{\partial y}{\partial \theta} @f$.
     * As the mapping is decomposed on B-splines, it means it computes
     * @f$ \frac{\partial y}{\partial \theta} (r,\theta) = \sum_k c_{y,k} \frac{\partial B_k}{\partial \theta}(r,\theta) @f$
     * (the derivatives of B-splines are implemented in SplineEvaluator2D).
     *
     * @param[in] coord
     * 				The coordinate where we evaluate the Jacobian matrix.
     *
     * @return A double with the value of the (2,2) coefficient of the Jacobian matrix.
     *
     * @see SplineEvaluator2D
     */
    double jacobian_22(ddc::Coordinate<circular_tag_r, circular_tag_p> const& coord) const final
    {
        return spline_evaluator.deriv_dim_2(coord, y_spline_representation.span_cview());
    }



    /**
     * @brief  Compute the full Jacobian matrix from the mapping to the pseudo-Cartesian mapping at the central point.
     *
     *
     * The discrete mappings can be difficult to inverse especially at the central point.
     * In case of non analytical invertible mapping, we can work in another domain called pseudo-Cartesian domain.
     * In this domain, it is easier to inverse the Jacobian matrix. The idea is detailed in Edoardo Zoni's article
     * (*Solving hyperbolic-elliptic problems on singular mapped disk-like domains with the method of characteristics and
     * spline finite elements*, https://doi.org/10.1016/j.jcp.2019.108889)
     *
     * The current mapping maps from the logical domain to the physical domain @f$ \mathcal{F}: (r,\theta) \mapsto (x, y)@f$.
     * The pseudo-Cartesian domain is built with @f$ \mathcal{F} @f$ and using the circular mapping @f$ \mathcal{G} @f$:
     * - @f$ \mathcal{G}_{11}(r, \theta) = \cos(\theta),
     *      \qquad\quad \mathcal{G}_{12}(r, \theta) = \sin(\theta) @f$,
     * - @f$ \mathcal{G}_{21}(r, \theta) = -\frac{1}{r}\sin(\theta),
     *      \qquad\quad \mathcal{G}_{22}(r, \theta) = \frac{1}{r}\cos(\theta) @f$.
     *
     * The pseudo-Cartesian domain is obtained by the composition of both mappings:
     * @f$ (\mathcal{F} \circ \mathcal{G}^{-1})^{-1} @f$.
     * This new mapping is invertible and its inverse at the central point is given by
     * - @f$ (J_{\mathcal{F}}J_{\mathcal{G}}^{-1})_{11}(0, \theta) = \partial_r x (0, \theta) \cos(\theta)
     *              - \partial_{r \theta} x (0, \theta) \sin(\theta), @f$
     * - @f$ (J_{\mathcal{F}}J_{\mathcal{G}}^{-1})_{12}(0, \theta) = \partial_r x (0, \theta) \sin(\theta)
     *              + \partial_{r \theta} x (0, \theta) \cos(\theta), @f$
     * - @f$ (J_{\mathcal{F}}J_{\mathcal{G}}^{-1})_{21}(0, \theta) = \partial_r y (0, \theta) \cos(\theta)
     *              - \partial_{r \theta} y (0, \theta) \sin(\theta), @f$
     * - @f$ (J_{\mathcal{F}}J_{\mathcal{G}}^{-1})_{22}(0, \theta) = \partial_r y (0, \theta) \sin(\theta)
     *              + \partial_{r \theta} y (0, \theta) \cos(\theta). @f$
     *
     *
     * So the pseudo-Cartesian Jacobian matrix at the central point,
     * @f$ (J_{\mathcal{F}}J_{\mathcal{G}}^{-1})^{-1}(0, \theta) @f$, is obtained by inversing this matrix.     *
     *
     *
     *
     * @param[in] grid
     *      The domain where the mapping is defined.
     * @param[out] matrix
     *      The pseudo-Cartesian matrix evaluated at the central point.
     *
     *
     * @see Curvilinear2DToCartesian
     * @see BslAdvection
     * @see AdvectionDomain
     */
    void to_pseudo_cartesian_jacobian_center_matrix(
            ddc::DiscreteDomain<IDimR, IDimP> const& grid,
            Matrix_2x2& matrix) const
    {
        ddc::DiscreteDomain<IDimP> const theta_domain = ddc::select<IDimP>(grid);

        matrix[0][0] = 0;
        matrix[0][1] = 0;
        matrix[1][0] = 0;
        matrix[1][1] = 0;

        // Average the values at (r = 0, theta):
        ddc::for_each(theta_domain, [&](auto const ip) {
            const double th = ddc::coordinate(ip);
            ddc::Coordinate<circular_tag_r, circular_tag_p> const coord(0, th);
            double const deriv_1_x
                    = spline_evaluator.deriv_dim_1(coord, x_spline_representation.span_cview());
            double const deriv_1_2_x
                    = spline_evaluator.deriv_1_and_2(coord, x_spline_representation.span_cview());
            double const deriv_1_y
                    = spline_evaluator.deriv_dim_1(coord, y_spline_representation.span_cview());
            double const deriv_1_2_y
                    = spline_evaluator.deriv_1_and_2(coord, y_spline_representation.span_cview());

            // Matrix from pseudo-Cart domain to physical domain by logical domain
            double const j11 = deriv_1_x * std::cos(th) - deriv_1_2_x * std::sin(th);
            double const j12 = deriv_1_x * std::sin(th) + deriv_1_2_x * std::cos(th);
            double const j21 = deriv_1_y * std::cos(th) - deriv_1_2_y * std::sin(th);
            double const j22 = deriv_1_y * std::sin(th) + deriv_1_2_y * std::cos(th);

            double const jacobian = j11 * j22 - j12 * j21;
            // Matrix from physical domain to pseudo_cart domain by logical domain
            if (fabs(jacobian) <= 1e-16) {
                std::cout << "WARNING! - Non invertible Jacobian matrix. ((r, theta) = (0, " << th
                          << "))" << std::endl;
            }
            assert(fabs(jacobian) >= 1e-16);

            matrix[0][0] += j22 / jacobian;
            matrix[0][1] += -j12 / jacobian;
            matrix[1][0] += -j21 / jacobian;
            matrix[1][1] += j11 / jacobian;
        });

        int const theta_size = theta_domain.size();
        matrix[0][0] /= theta_size;
        matrix[0][1] /= theta_size;
        matrix[1][0] /= theta_size;
        matrix[1][1] /= theta_size;
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
     * @brief Get a control point of the mapping on B-splines.
     *
     * The mapping @f$ (r,\theta) \mapsto (x,y) @f$ decomposed on B-splines can be
     * identified by its control points @f$ \{(c_{x,k}, c_{y,k})\}_{k}  @f$ where
     * @f$ c_{x,k} @f$ and @f$ c_{y,k} @f$ are the B-splines coefficients:
     *
     * @f$ x(r,\theta) = \sum_{k=0}^{N_r\times N_{\theta}-1} c_{x, k} B_k{r,\theta} @f$,
     *
     * @f$ y(r,\theta) = \sum_{k=0}^{N_r\times N_{\theta}-1} c_{y, k} B_k{r,\theta} @f$,
     *
     * where @f$ N_r\times N_{\theta} @f$ is the number of B-splines.
     *
     * The control points can be obtained by interpolating the mapping on interpolation
     * points (see GrevilleInterpolationPoints or KnotsAsInterpolationPoints).
     * We can also note that the first control points @f$ \{(c_{x,k}, c_{y,k})\}_{k=0}^{N_{\theta}} @f$
     * are equal to the pole @f$ (c_{x,k}, c_{y,k}) = (x_0,y_0) @f$, @f$ \forall k = 0, ..., N_{\theta}-1 @f$
     * where @f$ x(0,\theta), y(0,\theta) = (x_0,y_0) @f$ @f$ \forall \theta @f$.
     *
     *
     * @param[in] el
     * 			The number of the control point.
     *
     * @return The el-th control point.
     *
     * @see GrevilleInterpolationPoints
     * @see KnotsAsInterpolationPoints
     */
    inline const ddc::Coordinate<DimX, DimY> control_point(
            ddc::DiscreteElement<BSplineR, BSplineP> const& el) const
    {
        return ddc::
                Coordinate<DimX, DimY>(x_spline_representation(el), y_spline_representation(el));
    }


    /**
     * @brief Define a DiscreteToCartesian mapping from an analytical mapping.
     *
     * @param[in] analytical_mapping
     * 			The mapping defined analytically.
     * @param[in] builder
     * 			The spline builder on the B-splines on which we want to decompose the mapping.
     * @param[in] evaluator
     * 			The spline evaluator with which we want to evaluate the mapping.
     * @tparam Mapping
     * 			The analytical mapping described by this discrete mapping.
     *
     * @return A DiscreteToCartesian version of the analytical mapping.
     *
     * @see ddc::SplineBuilder2D
     * @see ddc::SplineEvaluator2D
     */
    template <class Mapping>
    static DiscreteToCartesian analytical_to_discrete(
            Mapping const& analytical_mapping,
            SplineBuilder const& builder,
            SplineEvaluator const& evaluator)
    {
        using Domain = typename SplineBuilder::interpolation_domain_type;
        SplineType curvilinear_to_x_spline(builder.spline_domain());
        SplineType curvilinear_to_y_spline(builder.spline_domain());
        ddc::Chunk<double, Domain> curvilinear_to_x_vals(builder.interpolation_domain());
        ddc::Chunk<double, Domain> curvilinear_to_y_vals(builder.interpolation_domain());
        ddc::for_each(
                builder.interpolation_domain(),
                [&](typename Domain::discrete_element_type const& el) {
                    ddc::Coordinate<circular_tag_r, circular_tag_p> polar_coord(
                            ddc::coordinate(el));
                    ddc::Coordinate<DimX, DimY> cart_coord = analytical_mapping(polar_coord);
                    curvilinear_to_x_vals(el) = ddc::select<DimX>(cart_coord);
                    curvilinear_to_y_vals(el) = ddc::select<DimY>(cart_coord);
                });
        builder(curvilinear_to_x_spline.span_view(), curvilinear_to_x_vals.span_cview());
        builder(curvilinear_to_y_spline.span_view(), curvilinear_to_y_vals.span_cview());
        return DiscreteToCartesian(
                std::move(curvilinear_to_x_spline),
                std::move(curvilinear_to_y_spline),
                evaluator);
    }
};
