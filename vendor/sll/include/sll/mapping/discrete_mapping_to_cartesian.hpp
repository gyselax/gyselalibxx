#pragma once

#include <ddc/ddc.hpp>

#include <sll/mapping/curvilinear2d_to_cartesian.hpp>
#include <sll/null_boundary_value.hpp>
#include <sll/spline_builder_2d.hpp>
#include <sll/spline_evaluator_2d.hpp>


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
template <class DimX, class DimY, class SplineBuilder>
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
     * @brief Define a 2x2 matrix with an 2D array of an 2D array.
     */
    using Matrix_2x2 = std::array<std::array<double, 2>, 2>;

private:
    using interpolation_domain = typename SplineBuilder::interpolation_domain_type;
    using spline_domain = ddc::DiscreteDomain<BSplineR, BSplineP>;

private:
    ddc::Chunk<double, spline_domain> x_spline_representation;
    ddc::Chunk<double, spline_domain> y_spline_representation;
    SplineEvaluator2D<BSplineR, BSplineP> const& spline_evaluator;

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
            ddc::Chunk<double, spline_domain>&& curvilinear_to_x,
            ddc::Chunk<double, spline_domain>&& curvilinear_to_y,
            SplineEvaluator2D<BSplineR, BSplineP> const& evaluator)
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
        const double x = spline_evaluator(coord, x_spline_representation);
        const double y = spline_evaluator(coord, y_spline_representation);
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
        matrix[0][0] = spline_evaluator.deriv_dim_1(coord, x_spline_representation);
        matrix[0][1] = spline_evaluator.deriv_dim_2(coord, x_spline_representation);
        matrix[1][0] = spline_evaluator.deriv_dim_1(coord, y_spline_representation);
        matrix[1][1] = spline_evaluator.deriv_dim_2(coord, y_spline_representation);
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
        return spline_evaluator.deriv_dim_1(coord, x_spline_representation);
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
        return spline_evaluator.deriv_dim_2(coord, x_spline_representation);
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
        return spline_evaluator.deriv_dim_1(coord, y_spline_representation);
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
        return spline_evaluator.deriv_dim_2(coord, y_spline_representation);
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
     *
     * @return A DiscreteToCartesian version of the analytical mapping.
     *
     * @see SplineBuilder2D
     * @see SplineEvaluator2D
     */
    template <class Mapping, class Builder2D, class Evaluator2D>
    static DiscreteToCartesian analytical_to_discrete(
            Mapping const& analytical_mapping,
            Builder2D const& builder,
            Evaluator2D const& evaluator)
    {
        using Domain = typename Builder2D::interpolation_domain_type;
        ddc::Chunk<double, spline_domain> curvilinear_to_x_spline(builder.spline_domain());
        ddc::Chunk<double, spline_domain> curvilinear_to_y_spline(builder.spline_domain());
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
        builder(curvilinear_to_x_spline, curvilinear_to_x_vals);
        builder(curvilinear_to_y_spline, curvilinear_to_y_vals);
        return DiscreteToCartesian(
                std::move(curvilinear_to_x_spline),
                std::move(curvilinear_to_y_spline),
                evaluator);
    }
};
