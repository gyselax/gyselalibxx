// SPDX-License-Identifier: MIT
#pragma once
#include <iostream>

#include <ddc/ddc.hpp>

#include "ddc_alias_inline_functions.hpp"
#include "ddc_aliases.hpp"
#include "mapping_tools.hpp"
#include "math_tools.hpp"
#include "tensor.hpp"
#include "view.hpp"

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
 */
template <
        class X,
        class Y,
        class SplineEvaluator,
        class R = typename SplineEvaluator::continuous_dimension_type1,
        class Theta = typename SplineEvaluator::continuous_dimension_type2,
        class MemorySpace = typename SplineEvaluator::memory_space>
class DiscreteToCartesian
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
    using cartesian_tag_x = X;
    /// @brief Indicate the second physical coordinate.
    using cartesian_tag_y = Y;
    /// @brief Indicate the first logical coordinate.
    using curvilinear_tag_r = R;
    /// @brief Indicate the second logical coordinate.
    using curvilinear_tag_theta = Theta;

    /// The type of the argument of the function described by this mapping
    using CoordArg = Coord<R, Theta>;
    /// The type of the result of the function described by this mapping
    using CoordResult = Coord<X, Y>;

    /// @brief The covariant form of the first physical coordinate.
    using X_cov = typename X::Dual;
    /// @brief The covariant form of the second physical coordinate.
    using Y_cov = typename Y::Dual;
    /// @brief The covariant form of the first logical coordinate.
    using R_cov = typename R::Dual;
    /// @brief The covariant form of the second logical coordinate.
    using Theta_cov = typename Theta::Dual;

private:
    using spline_idx_range = IdxRange<BSplineR, BSplineTheta>;

    using SplineType = DConstField<spline_idx_range, MemorySpace>;

    using IdxRangeRTheta = typename SplineEvaluator::evaluation_domain_type;
    using IdxRangeTheta = typename SplineEvaluator::evaluation_domain_type2;
    using IdxTheta = typename IdxRangeTheta::discrete_element_type;

private:
    SplineType m_x_spline_representation;
    SplineType m_y_spline_representation;
    SplineEvaluator m_spline_evaluator;
    IdxRangeRTheta m_idx_range_singular_point;

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
     * @param[in] idx_range_singular_point
     *      The index range describing the points which should be used to evaluate functions at the central point.
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
            IdxRangeRTheta idx_range_singular_point)
        : m_x_spline_representation(curvilinear_to_x)
        , m_y_spline_representation(curvilinear_to_y)
        , m_spline_evaluator(evaluator)
        , m_idx_range_singular_point(idx_range_singular_point)
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
    KOKKOS_FUNCTION Coord<X, Y> operator()(
            Coord<curvilinear_tag_r, curvilinear_tag_theta> const& coord) const
    {
        const double x = m_spline_evaluator(coord, get_const_field(m_x_spline_representation));
        const double y = m_spline_evaluator(coord, get_const_field(m_y_spline_representation));
        return Coord<X, Y>(x, y);
    }

    /**
     * @brief Compute full Jacobian matrix.
     *
     * For some computations, we need the complete Jacobian matrix or just the
     * coefficients.
     * The coefficients can be given independently with the functions
     * jacobian_11, jacobian_12, jacobian_21 and jacobian_22.
     *
     * @param[in] coord
     * 				The coordinate where we evaluate the Jacobian matrix.
     * @param[out] matrix
     * 				The Jacobian matrix returned.
     */
    KOKKOS_FUNCTION DTensor<VectorIndexSet<X, Y>, VectorIndexSet<R_cov, Theta_cov>> jacobian_matrix(
            Coord<R, Theta> const& coord) const
    {
        DTensor<VectorIndexSet<X, Y>, VectorIndexSet<R_cov, Theta_cov>> jacobian_matrix;
        ddcHelper::get<X, R_cov>(jacobian_matrix)
                = m_spline_evaluator.deriv_dim_1(coord, get_const_field(m_x_spline_representation));
        ddcHelper::get<X, Theta_cov>(jacobian_matrix)
                = m_spline_evaluator.deriv_dim_2(coord, get_const_field(m_x_spline_representation));
        ddcHelper::get<Y, R_cov>(jacobian_matrix)
                = m_spline_evaluator.deriv_dim_1(coord, get_const_field(m_y_spline_representation));
        ddcHelper::get<Y, Theta_cov>(jacobian_matrix)
                = m_spline_evaluator.deriv_dim_2(coord, get_const_field(m_y_spline_representation));
        return jacobian_matrix;
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
            Coord<curvilinear_tag_r, curvilinear_tag_theta> const& coord) const
    {
        return m_spline_evaluator.deriv_dim_1(coord, get_const_field(m_x_spline_representation));
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
            Coord<curvilinear_tag_r, curvilinear_tag_theta> const& coord) const
    {
        return m_spline_evaluator.deriv_dim_2(coord, get_const_field(m_x_spline_representation));
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
            Coord<curvilinear_tag_r, curvilinear_tag_theta> const& coord) const
    {
        return m_spline_evaluator.deriv_dim_1(coord, get_const_field(m_y_spline_representation));
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
            Coord<curvilinear_tag_r, curvilinear_tag_theta> const& coord) const
    {
        return m_spline_evaluator.deriv_dim_2(coord, get_const_field(m_y_spline_representation));
    }

    /**
     * @brief Compute the Jacobian, the determinant of the Jacobian matrix of the mapping.
     *
     * @param[in] coord
     *          The coordinate where we evaluate the Jacobian.
     *
     * @return A double with the value of the determinant of the Jacobian matrix.
     */
    KOKKOS_FUNCTION double jacobian(
            Coord<curvilinear_tag_r, curvilinear_tag_theta> const& coord) const
    {
        Tensor J = jacobian_matrix(coord);
        return ddcHelper::get<X, R_cov>(J) * ddcHelper::get<Y, Theta_cov>(J)
               - ddcHelper::get<Y, R_cov>(J) * ddcHelper::get<X, Theta_cov>(J);
    }

    /**
     * @brief Get the first order expansion of the Jacobian matrix with the theta component divided by r.
     * The expansion is carried out around @f$ r=0 @f$.
     * The returned matrix @f$ J@f$ is defined as:
     * @f$ J_{00} = \frac{\partial x}{\partial r}(r, \theta) @f$
     * @f$ J_{01} = \frac{1}{r} \frac{\partial x}{\partial \theta}(r, \theta) + O(r^2) = \frac{\partial^2 x}{\partial r \partial \theta}(0, \theta) @f$
     * @f$ J_{10} = \frac{\partial y}{\partial r}(r, \theta) @f$
     * @f$ J_{11} = \frac{1}{r} \frac{\partial y}{\partial \theta}(r, \theta) + O(r^2) = \frac{\partial^2 y}{\partial r \partial \theta}(0, \theta) @f$
     *
     * @param[in] coord
     *          The coordinate where we evaluate the Jacobian.
     *
     * @return The first order expansion of the Jacobian matrix with the theta component divided by r.
     */
    KOKKOS_INLINE_FUNCTION DTensor<VectorIndexSet<X, Y>, VectorIndexSet<R_cov, Theta_cov>>
    first_order_jacobian_matrix_r_rtheta(
            Coord<curvilinear_tag_r, curvilinear_tag_theta> const& coord) const
    {
        DTensor<VectorIndexSet<X, Y>, VectorIndexSet<R_cov, Theta_cov>> jacobian_matrix;
        ddcHelper::get<X, R_cov>(jacobian_matrix)
                = m_spline_evaluator.deriv_dim_1(coord, get_const_field(m_x_spline_representation));
        ddcHelper::get<X, Theta_cov>(jacobian_matrix)
                = m_spline_evaluator
                          .deriv_1_and_2(coord, get_const_field(m_x_spline_representation));
        ddcHelper::get<Y, R_cov>(jacobian_matrix)
                = m_spline_evaluator.deriv_dim_1(coord, get_const_field(m_y_spline_representation));
        ddcHelper::get<X, Theta_cov>(jacobian_matrix)
                = m_spline_evaluator
                          .deriv_1_and_2(coord, get_const_field(m_y_spline_representation));
        return jacobian_matrix;
    }

    /**
     * @brief Get the index range describing the points which should be used to evaluate functions at the central point.
     *
     * @return An index range covering the O-point.
     */
    KOKKOS_INLINE_FUNCTION IdxRangeRTheta idx_range_singular_point() const
    {
        return m_idx_range_singular_point;
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
    KOKKOS_INLINE_FUNCTION const Coord<X, Y> control_point(
            Idx<BSplineR, BSplineTheta> const& el) const
    {
        return Coord<X, Y>(m_x_spline_representation(el), m_y_spline_representation(el));
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

template <class X, class Y, class SplineEvaluator, class R, class Theta, class MemorySpace>
struct IsCurvilinear2DMapping<DiscreteToCartesian<X, Y, SplineEvaluator, R, Theta, MemorySpace>>
    : std::true_type
{
};

template <class X, class Y, class SplineEvaluator, class R, class Theta, class MemorySpace>
struct SingularOPointInvJacobian<DiscreteToCartesian<X, Y, SplineEvaluator, R, Theta, MemorySpace>>
    : std::true_type
{
};

} // namespace mapping_detail
