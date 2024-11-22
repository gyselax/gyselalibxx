// SPDX-License-Identifier: MIT
#pragma once
#include <iostream>

#include <ddc/ddc.hpp>

#include <sll/view.hpp>

#include "coordinate_converter.hpp"
#include "curvilinear2d_to_cartesian.hpp"
#include "jacobian.hpp"
#include "mapping_tools.hpp"
#include "pseudo_cartesian_compatible_mapping.hpp"

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
template <
        class X,
        class Y,
        class SplineEvaluator,
        class R = typename SplineEvaluator::continuous_dimension_type1,
        class Theta = typename SplineEvaluator::continuous_dimension_type2,
        class MemorySpace = typename SplineEvaluator::memory_space>
class DiscreteToCartesian
    : public CoordinateConverter<ddc::Coordinate<R, Theta>, ddc::Coordinate<X, Y>>
    , public NonAnalyticalJacobian<ddc::Coordinate<R, Theta>>
    , public PseudoCartesianCompatibleMapping
    , public Curvilinear2DToCartesian<X, Y, R, Theta>
{
    static_assert(std::is_same_v<MemorySpace, typename SplineEvaluator::memory_space>);

public:
    /**
     * @brief Indicate the bspline type of the first logical dimension.
     */
    using BSplineR = typename SplineEvaluator::bsplines_type1;
    /**
     * @brief Indicate the bspline type of the second logical dimension.
     */
    using BSplineTheta = typename SplineEvaluator::bsplines_type2;

    /// @brief Indicate the first physical coordinate.
    using cartesian_tag_x = typename Curvilinear2DToCartesian<X, Y, R, Theta>::cartesian_tag_x;
    /// @brief Indicate the second physical coordinate.
    using cartesian_tag_y = typename Curvilinear2DToCartesian<X, Y, R, Theta>::cartesian_tag_y;
    /// @brief Indicate the first logical coordinate.
    using curvilinear_tag_r = typename Curvilinear2DToCartesian<X, Y, R, Theta>::curvilinear_tag_r;
    /// @brief Indicate the second logical coordinate.
    using curvilinear_tag_theta =
            typename Curvilinear2DToCartesian<X, Y, R, Theta>::curvilinear_tag_theta;

    /// The type of the argument of the function described by this mapping
    using CoordArg = ddc::Coordinate<R, Theta>;
    /// The type of the result of the function described by this mapping
    using CoordResult = ddc::Coordinate<X, Y>;

private:
    using spline_idx_range = ddc::DiscreteDomain<BSplineR, BSplineTheta>;

    using SplineType = ddc::
            ChunkView<double, spline_idx_range, std::experimental::layout_right, MemorySpace>;

    using IdxRangeTheta = typename SplineEvaluator::evaluation_domain_type2;
    using IdxTheta = typename IdxRangeTheta::discrete_element_type;

private:
    SplineType m_x_spline_representation;
    SplineType m_y_spline_representation;
    SplineEvaluator m_spline_evaluator;
    IdxRangeTheta m_idx_range_theta;

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
     * @param[in] idx_range_theta
     *      The index range describing the poloidal points which should be used to average the derivatives of
     *      the pseudo-Cartesian matrix at the central point.
     *
     *
     * @see SplineBuilder2D
     * @see DiscreteToCartesian::operator()
     * @see SplineBoundaryValue
     */
    KOKKOS_FUNCTION DiscreteToCartesian(
            SplineType curvilinear_to_x,
            SplineType curvilinear_to_y,
            SplineEvaluator const& evaluator,
            IdxRangeTheta idx_range_theta)
        : m_x_spline_representation(curvilinear_to_x)
        , m_y_spline_representation(curvilinear_to_y)
        , m_spline_evaluator(evaluator)
        , m_idx_range_theta(idx_range_theta)
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
    KOKKOS_FUNCTION ddc::Coordinate<X, Y> operator()(
            ddc::Coordinate<curvilinear_tag_r, curvilinear_tag_theta> const& coord) const final
    {
        const double x = m_spline_evaluator(coord, m_x_spline_representation.span_cview());
        const double y = m_spline_evaluator(coord, m_y_spline_representation.span_cview());
        return ddc::Coordinate<X, Y>(x, y);
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
    KOKKOS_FUNCTION void jacobian_matrix(
            ddc::Coordinate<curvilinear_tag_r, curvilinear_tag_theta> const& coord,
            Matrix_2x2& matrix) const final
    {
        matrix[0][0]
                = m_spline_evaluator.deriv_dim_1(coord, m_x_spline_representation.span_cview());
        matrix[0][1]
                = m_spline_evaluator.deriv_dim_2(coord, m_x_spline_representation.span_cview());
        matrix[1][0]
                = m_spline_evaluator.deriv_dim_1(coord, m_y_spline_representation.span_cview());
        matrix[1][1]
                = m_spline_evaluator.deriv_dim_2(coord, m_y_spline_representation.span_cview());
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
    KOKKOS_FUNCTION double jacobian_11(
            ddc::Coordinate<curvilinear_tag_r, curvilinear_tag_theta> const& coord) const final
    {
        return m_spline_evaluator.deriv_dim_1(coord, m_x_spline_representation.span_cview());
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
    KOKKOS_FUNCTION double jacobian_12(
            ddc::Coordinate<curvilinear_tag_r, curvilinear_tag_theta> const& coord) const final
    {
        return m_spline_evaluator.deriv_dim_2(coord, m_x_spline_representation.span_cview());
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
    KOKKOS_FUNCTION double jacobian_21(
            ddc::Coordinate<curvilinear_tag_r, curvilinear_tag_theta> const& coord) const final
    {
        return m_spline_evaluator.deriv_dim_1(coord, m_y_spline_representation.span_cview());
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
    KOKKOS_FUNCTION double jacobian_22(
            ddc::Coordinate<curvilinear_tag_r, curvilinear_tag_theta> const& coord) const final
    {
        return m_spline_evaluator.deriv_dim_2(coord, m_y_spline_representation.span_cview());
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
     * @f$ (J_{\mathcal{F}}J_{\mathcal{G}}^{-1})^{-1}(0, \theta) @f$, is obtained by inversing this matrix. 
     *
     *
     *
     * @param[out] matrix
     *      The pseudo-Cartesian matrix evaluated at the central point.
     *
     *
     * @see Curvilinear2DToCartesian
     * @see BslAdvection
     * @see AdvectionDomain
     */
    KOKKOS_FUNCTION void to_pseudo_cartesian_jacobian_center_matrix(Matrix_2x2& matrix) const final
    {
        matrix[0][0] = 0;
        matrix[0][1] = 0;
        matrix[1][0] = 0;
        matrix[1][1] = 0;

        // Average the values at (r = 0, theta):
        for (IdxTheta ip : m_idx_range_theta) {
            const double th = ddc::coordinate(ip);
            ddc::Coordinate<curvilinear_tag_r, curvilinear_tag_theta> const coord(0, th);
            double const deriv_1_x
                    = m_spline_evaluator.deriv_dim_1(coord, m_x_spline_representation.span_cview());
            double const deriv_1_2_x
                    = m_spline_evaluator
                              .deriv_1_and_2(coord, m_x_spline_representation.span_cview());
            double const deriv_1_y
                    = m_spline_evaluator.deriv_dim_1(coord, m_y_spline_representation.span_cview());
            double const deriv_1_2_y
                    = m_spline_evaluator
                              .deriv_1_and_2(coord, m_y_spline_representation.span_cview());

            // Matrix from pseudo-Cart domain to physical domain by logical domain
            double const j11 = deriv_1_x * Kokkos::cos(th) - deriv_1_2_x * Kokkos::sin(th);
            double const j12 = deriv_1_x * Kokkos::sin(th) + deriv_1_2_x * Kokkos::cos(th);
            double const j21 = deriv_1_y * Kokkos::cos(th) - deriv_1_2_y * Kokkos::sin(th);
            double const j22 = deriv_1_y * Kokkos::sin(th) + deriv_1_2_y * Kokkos::cos(th);

            double const jacobian = j11 * j22 - j12 * j21;
            // Matrix from physical domain to pseudo_cart domain by logical domain
            assert(fabs(jacobian) >= 1e-16);

            matrix[0][0] += j22 / jacobian;
            matrix[0][1] += -j12 / jacobian;
            matrix[1][0] += -j21 / jacobian;
            matrix[1][1] += j11 / jacobian;
        }

        int const theta_size = m_idx_range_theta.size();
        matrix[0][0] /= theta_size;
        matrix[0][1] /= theta_size;
        matrix[1][0] /= theta_size;
        matrix[1][1] /= theta_size;
    }


    /**
     * @brief Compute the (1,1) coefficient of the pseudo-Cartesian Jacobian matrix at the central point.
     *
     * @return A double with the (1,1) coefficient of the pseudo-Cartesian Jacobian matrix at the central point.
     *
     * @see to_pseudo_cartesian_jacobian_center_matrix
     */
    KOKKOS_FUNCTION double to_pseudo_cartesian_jacobian_11_center() const final
    {
        Matrix_2x2 jacobian;
        to_pseudo_cartesian_jacobian_center_matrix(jacobian);
        return jacobian[0][0];
    }


    /**
     * @brief Compute the (1,2) coefficient of the pseudo-Cartesian Jacobian matrix at the central point.
     *
     * @return A double with the (1,2) coefficient of the pseudo-Cartesian Jacobian matrix at the central point.
     *
     * @see to_pseudo_cartesian_jacobian_center_matrix
     */
    KOKKOS_FUNCTION double to_pseudo_cartesian_jacobian_12_center() const final
    {
        Matrix_2x2 jacobian;
        to_pseudo_cartesian_jacobian_center_matrix(jacobian);
        return jacobian[0][1];
    }


    /**
     * @brief Compute the (2,1) coefficient of the pseudo-Cartesian Jacobian matrix at the central point.
     *
     * @return A double with the (2,1) coefficient of the pseudo-Cartesian Jacobian matrix at the central point.
     *
     * @see to_pseudo_cartesian_jacobian_center_matrix
     */
    KOKKOS_FUNCTION double to_pseudo_cartesian_jacobian_21_center() const final
    {
        Matrix_2x2 jacobian;
        to_pseudo_cartesian_jacobian_center_matrix(jacobian);
        return jacobian[1][0];
    }


    /**
     * @brief Compute the (2,2) coefficient of the pseudo-Cartesian Jacobian matrix at the central point.
     *
     * @return A double with the (2,2) coefficient of the pseudo-Cartesian Jacobian matrix at the central point.
     *
     * @see to_pseudo_cartesian_jacobian_center_matrix
     */
    KOKKOS_FUNCTION double to_pseudo_cartesian_jacobian_22_center() const final
    {
        Matrix_2x2 jacobian;
        to_pseudo_cartesian_jacobian_center_matrix(jacobian);
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
    KOKKOS_INLINE_FUNCTION const ddc::Coordinate<X, Y> control_point(
            ddc::DiscreteElement<BSplineR, BSplineTheta> const& el) const
    {
        return ddc::Coordinate<X, Y>(m_x_spline_representation(el), m_y_spline_representation(el));
    }
};


namespace mapping_detail {
template <
        class X,
        class Y,
        class SplineEvaluator,
        class R,
        class Theta,
        class MemorySpace,
        class ExecSpace>
struct MappingAccessibility<
        ExecSpace,
        DiscreteToCartesian<X, Y, SplineEvaluator, R, Theta, MemorySpace>>
{
    static constexpr bool value = Kokkos::SpaceAccessibility<ExecSpace, MemorySpace>::accessible;
};
} // namespace mapping_detail
