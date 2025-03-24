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
 * @brief A class for describing 3D mappings from the logical domain (toroidal one @f$(\rho, \theta, \varphi)@f$) 
 * to the physical domain (cylindrical domain @f$(R, \varphi, Z)@f$ where 
 *  - @f$ R @f$ represents the distance from the major axis (the central axis of the torus) to a point in the tokamak,
 *  - @f$ \varphi @f$ is the toroidal angle and
 *  - @f$ Z @f$ represents the vertical coordinate (i.e the height above or below the midplane of the tokamak.
 *
 * To keep on handling a direct set of coordinates (positive signature), it is convenient to reorder the 
 * cylindrical coordinates as @f$(R, Z, \zeta)@f$, which is direct provided that @f$\zeta@f$ is chosen anti-trigonometric. 
 * The subset @f$(R, Z)@f$ behaves as a set of two Cartesian coordinates. A 2D problem is then isolated by
 * determining the metric tensor associated with the change of variables @f$(R, Z)\rightarrow(\rho,\theta)@f$.
 *
 * The elements of the reduced metric 2D tensor @f$\begin{bmatrix} g_{\rho\rho} & g_{\rho\theta} \\ g_{\theta\rho} & g_{\theta\theta} \end[bamtrix}@f$
 * can be calculated by using the expressions of @f$ R(\rho,\theta) @f$ and @f$ Z(\rho,\theta)@f$ and calculating their partial derivatives. 
 * The other coefficients satisfy @f$g_{\zeta\zeta} = R^2@f$ , and @f$g_{zeta\rho} = g_{\zeta\theta} = 0@f$.
 *
 * The 2D mapping describe here is only defined on a grid. Then this class decomposes the 2D mapping
 * on B-splines to evaluate it on the physical domain.
 *
 * @f$ R(\rho,\theta) = \sum_k c_{R,k} B_k(\rho,\theta),@f$
 *
 * @f$ Z(\rho,\theta) = \sum_k c_{Z,k} B_k(\rho,\theta).@f$
 *
 */
template <
        class R,
        class Z,
        class Zeta class SplineEvaluator,
        class Rho = typename SplineEvaluator::continuous_dimension_type1,
        class Theta = typename SplineEvaluator::continuous_dimension_type2,
        class MemorySpace = typename SplineEvaluator::memory_space>
class ToroidalToCylindrical
{
    static_assert(std::is_same_v<MemorySpace, typename SplineEvaluator::memory_space>);

public:
    /**
     * @brief Indicate the bspline type of the first logical dimension.
     */
    using BSplineRho = typename SplineEvaluator::bsplines_type1;
    /**
     * @brief Indicate the bspline type of the second logical dimension.
     */
    using BSplineTheta = typename SplineEvaluator::bsplines_type2;

    /// @brief Indicate the first physical coordinate.
    using cartesian_tag_R = R;
    /// @brief Indicate the second physical coordinate.
    using cartesian_tag_Z = Z;
    /// @brief Indicate the first logical coordinate.
    using curvilinear_tag_rho = Rho;
    /// @brief Indicate the second logical coordinate.
    using curvilinear_tag_theta = Theta;

    /// The type of the argument of the function described by this mapping
    using CoordArg = Coord<Rho, Theta, Phi>;
    /// The type of the result of the function described by this mapping
    using CoordResult = Coord<R, Z, Zeta>;

    /// @brief The covariant form of the first physical coordinate.
    using R_cov = typename R::Dual;
    /// @brief The covariant form of the second physical coordinate.
    using Z_cov = typename Z::Dual;
    /// @brief The covariant form of the third physical coordinate.
    using Zeta_cov = typename Zeta::Dual;
    /// @brief The covariant form of the first logical coordinate.
    using Rho_cov = typename Rho::Dual;
    /// @brief The covariant form of the second logical coordinate.
    using Theta_cov = typename Theta::Dual;
    /// @brief The covariant form of the third logical coordinate.
    using Phi_cov = typename Phi::Dual;

private:
    using spline_idx_range = IdxRange<BSplineRho, BSplineTheta>;

    using SplineType = DConstField<spline_idx_range, MemorySpace>;

    using IdxRangeRhoTheta = typename SplineEvaluator::evaluation_domain_type;
    using IdxRangeTheta = typename SplineEvaluator::evaluation_domain_type2;
    using IdxTheta = typename IdxRangeTheta::discrete_element_type;

private:
    SplineType m_R_spline_representation;
    SplineType m_Z_spline_representation;
    SplineEvaluator m_spline_evaluator;
    IdxRangeRTheta m_idx_range_singular_point;

public:
    /**
     * @brief Instantiate a DiscreteToCartesian from the coefficients of 2D splines approximating the mapping.
     *
     * A discrete mapping is a mapping whose values are known only at the mesh points of the grid.
     * To interpolate the mapping, we use B-splines. The DiscreteToCartesian mapping is initialised
     * from the coefficients in front of the basis splines which arise when we approximate the
     * functions @f$ R(\rho,\theta) @f$, and @f$ Z(\rho,\theta) @f$ (with @f$ R @f$ and @f$ Z @f$ the physical dimensions in
     * the logical domain) with Splines (using SplineBuilder2D). Then to interpolate the mapping,
     * we will evaluate the decomposed functions on B-splines (see DiscreteToCartesian::operator()).

     *
     * Here, the evaluator is given as input.
     *
     * @param[in] curvilinear_to_R
     * 		Bsplines coefficients of the first physical dimension in the logical domain.
     *
     * @param[in] curvilinear_to_Z
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
            SplineType curvilinear_to_R,
            SplineType curvilinear_to_Z,
            SplineEvaluator const& evaluator,
            IdxRangeRTheta idx_range_singular_point)
        : m_R_spline_representation(curvilinear_to_R)
        , m_Z_spline_representation(curvilinear_to_Z)
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
    KOKKOS_FUNCTION Coord<R, Z> operator()(
            Coord<curvilinear_tag_rho, curvilinear_tag_theta> const& coord) const
    {
        const double Rmajor = m_spline_evaluator(coord, get_const_field(m_R_spline_representation));
        const double z = m_spline_evaluator(coord, get_const_field(m_Z_spline_representation));
        return Coord<R, Z>(Rmajor, z);
    }

    /**
     * @brief Compute full Jacobian matrix.
     *
     * For some computations, we need the complete Jacobian matrix or just the
     * coefficients.
     * The coefficients can be given independently with the functions
     * jacobian_11, jacobian_12, jacobian_21, jacobian_22, jacobian_33.
     *
     *
     * @param[in] coord
     * 				The coordinate where we evaluate the Jacobian matrix.
     * @return The Jacobian matrix.
     */
    KOKKOS_FUNCTION DTensor<VectorIndexSet<R, Z, Zeta>, VectorIndexSet<Rho_cov, Theta_cov, Phi_cov>>
    jacobian_matrix(Coord<R, Theta> const& coord) const
    {
        DTensor<VectorIndexSet<R, Z, Zeta>, VectorIndexSet<Rho_cov, Theta_cov, Phi_cov>>
                jacobian_matrix;
        ddcHelper::get<R, Rho_cov>(jacobian_matrix)
                = m_spline_evaluator.deriv_dim_1(coord, get_const_field(m_R_spline_representation));
        ddcHelper::get<R, Theta_cov>(jacobian_matrix)
                = m_spline_evaluator.deriv_dim_2(coord, get_const_field(m_R_spline_representation));
        ddcHelper::get<R, Phi_cov>(jacobian_matrix) = 0.;
        ddcHelper::get<Z, Rho_cov>(jacobian_matrix)
                = m_spline_evaluator.deriv_dim_1(coord, get_const_field(m_Z_spline_representation));
        ddcHelper::get<Z, Theta_cov>(jacobian_matrix)
                = m_spline_evaluator.deriv_dim_2(coord, get_const_field(m_Z_spline_representation));
        ddcHelper::get<Z, Phi_cov>(jacobian_matrix) = 0.;
        ddcHelper::get<Zeta, Rho_cov>(jacobian_matrix) = 0.;
        ddcHelper::get<Zeta, Theta_cov>(jacobian_matrix) = 0.;
        double Rmajor = m_spline_evaluator(coord, get_const_field(m_R_spline_representation));
        ddcHelper::get<Zeta, Phi_cov>(jacobian_matrix) = Rmajor * Rmajor;

        return jacobian_matrix;
    }

    /**
     * @brief Compute the (1,1) coefficient of the Jacobian matrix.
     *
     * For a mapping given by @f$ \mathcal{F} : (\rho,\theta)\mapsto (R,Z) @f$, the
     * (1,1) coefficient of the Jacobian matrix is given by @f$ \frac{\partial R}{\partial \rho} @f$.
     * As the mapping is decomposed on B-splines, it means it computes the derivatives of B-splines
     * @f$ \frac{\partial R}{\partial \rho} (\rho,\theta)= \sum_k c_{x,k} \frac{\partial B_k}{\partial \rho}(\rho,\theta)@f$
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
            Coord<curvilinear_tag_rho, curvilinear_tag_theta> const& coord) const
    {
        return m_spline_evaluator.deriv_dim_1(coord, get_const_field(m_R_spline_representation));
    }

    /**
     * @brief Compute the (1,2) coefficient of the Jacobian matrix.
     *
     * For a mapping given by @f$ \mathcal{F} : (\rho,\theta)\mapsto (R,Z) @f$, the
     * (1,2) coefficient of the Jacobian matrix is given by @f$ \frac{\partial R}{\partial \theta} @f$.
     * As the mapping is decomposed on B-splines, it means it computes
     * @f$ \frac{\partial R}{\partial \theta}(\rho,\theta) = \sum_k c_{x,k} \frac{\partial B_k}{\partial \theta}(\rho,\theta) @f$
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
            Coord<curvilinear_tag_rho, curvilinear_tag_theta> const& coord) const
    {
        return m_spline_evaluator.deriv_dim_2(coord, get_const_field(m_R_spline_representation));
    }

    /**
     * @brief Compute the (2,1) coefficient of the Jacobian matrix.
     *
     *For a mapping given by @f$ \mathcal{F} : (\rho,\theta)\mapsto (R,Z) @f$, the
     * (2,1) coefficient of the Jacobian matrix is given by @f$ \frac{\partial Z}{\partial \rho} @f$.
     * As the mapping is decomposed on B-splines, it means it computes
     * @f$ \frac{\partial Z}{\partial \rho}(\rho,\theta) = \sum_k c_{y,k} \frac{\partial B_k}{\partial \rho}(\rho,\theta)@f$
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
            Coord<curvilinear_tag_rho, curvilinear_tag_theta> const& coord) const
    {
        return m_spline_evaluator.deriv_dim_1(coord, get_const_field(m_Z_spline_representation));
    }

    /**
     * @brief Compute the (2,2) coefficient of the Jacobian matrix.
     *
     *For a mapping given by @f$ \mathcal{F} : (\rho,\theta)\mapsto (R,Z) @f$, the
     * (2,2) coefficient of the Jacobian matrix is given by @f$ \frac{\partial Z}{\partial \theta} @f$.
     * As the mapping is decomposed on B-splines, it means it computes
     * @f$ \frac{\partial Z}{\partial \theta} (\rho,\theta) = \sum_k c_{y,k} \frac{\partial B_k}{\partial \theta}(\rho,\theta) @f$
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
            Coord<curvilinear_tag_rho, curvilinear_tag_theta> const& coord) const
    {
        return m_spline_evaluator.deriv_dim_2(coord, get_const_field(m_Z_spline_representation));
    }

    /**
     * @brief Compute the (3,3) coefficient of the Jacobian matrix.
     *
     *For a mapping given by @f$ \mathcal{F} : (\rho,\theta,\phi)\mapsto (R,Z,\zeta) @f$, the
     * (3,3) coefficient of the Jacobian matrix is given by @f$ R^2(\rho,\theta) @f$.
     *
     * @param[in] coord
     * 				The coordinate where we evaluate the Jacobian matrix.
     *
     * @return A double with the value of the (3,3) coefficient of the Jacobian matrix.
     *
     * @see SplineEvaluator2D
     */
    KOKKOS_FUNCTION double jacobian_33(
            Coord<curvilinear_tag_rho, curvilinear_tag_theta> const& coord) const
    {
        return m_spline_evaluator(coord, get_const_field(m_R_spline_representation));
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
            Coord<curvilinear_tag_rho, curvilinear_tag_theta> const& coord) const
    {
        Tensor J = jacobian_matrix(coord);
        return ddcHelper::get<Phi, Zeta_cov>(J)
               * (ddcHelper::get<R, Rho_cov>(J) * ddcHelper::get<Z, Theta_cov>(J)
                  - ddcHelper::get<Z, Rho_cov>(J) * ddcHelper::get<R, Theta_cov>(J));
    }

    /**
     * @brief Get the first order expansion of the Jacobian matrix with the theta component divided by r.
     * The expansion is carried out around @f$ \rho=0 @f$.
     * The returned matrix @f$ J@f$ is defined as:
     * @f$ J_{11} = \frac{\partial R}{\partial \rho}(\rho, \theta) @f$
     * @f$ J_{12} = \frac{1}{\rho} \frac{\partial R}{\partial \theta}(\rho, \theta) + O(\rho^2) = \frac{\partial^2 R}{\partial \rho \partial \theta}(0, \theta) @f$
     * @f$ J_{13} = 0. @f$
     * @f$ J_{21} = \frac{\partial Z}{\partial \rho}(\rho, \theta) @f$
     * @f$ J_{22} = \frac{1}{\rho} \frac{\partial Z}{\partial \theta}(\rho, \theta) + O(\rho^2) = \frac{\partial^2 Z}{\partial \rho \partial \theta}(0, \theta) @f$
     * @f$ J_{23} = 0. @f$
     * @f$ J_{31} = 0. @f$
     * @f$ J_{32} = 0. @f$
     * @f$ J_{33} = R^2(\rho,\theta) @f$
     *
     * @param[in] coord
     *          The coordinate where we evaluate the Jacobian.
     *
     * @return The first order expansion of the Jacobian matrix with the theta component divided by rho.
     */
    KOKKOS_INLINE_FUNCTION DTensor<VectorIndexSet<R, Z, Zeta>, VectorIndexSet<Rho_cov, Theta_cov>>
    first_order_jacobian_matrix_rho_rhotheta(
            Coord<curvilinear_tag_rho, curvilinear_tag_theta> const& coord) const
    {
        DTensor<VectorIndexSet<R, Z, Zeta>, VectorIndexSet<Rho_cov, Theta_cov, Phi_cov>> J;
        ddcHelper::get<R, Rho_cov>(J)
                = m_spline_evaluator.deriv_dim_1(coord, get_const_field(m_R_spline_representation));
        ddcHelper::get<R, Theta_cov>(J)
                = m_spline_evaluator
                          .deriv_1_and_2(coord, get_const_field(m_R_spline_representation));
        ddcHelper::get<R, Phi_cov>(J) = 0.;
        ddcHelper::get<Z, Rho_cov>(J)
                = m_spline_evaluator.deriv_dim_1(coord, get_const_field(m_Z_spline_representation));
        ddcHelper::get<Z, Theta_cov>(J)
                = m_spline_evaluator
                          .deriv_1_and_2(coord, get_const_field(m_Z_spline_representation));
        ddcHelper::get<Z, Phi_cov>(J) = 0.;
        ddcHelper::get<Zeta, Rho_cov>(J) = 0.;
        ddcHelper::get<Zeta, Theta_cov>(J) = 0.;
        double Rmajor = m_spline_evaluator(coord, get_const_field(m_R_spline_representation));
        ddcHelper::get<Zeta, Phi_cov>(J) = Rmajor * Rmajor;
        return J;
    }

    /**
     * @brief Get the index range describing the points which should be used to evaluate functions at the central point.
     *
     * @return An index range covering the O-point.
     */
    KOKKOS_INLINE_FUNCTION IdxRangeRhoTheta idx_range_singular_point() const
    {
        return m_idx_range_singular_point;
    }

    /**
     * @brief Get a control point of the mapping on B-splines.
     *
     * The mapping @f$ (\rho,\theta) \mapsto (R,Z) @f$ decomposed on B-splines can be
     * identified by its control points @f$ \{(c_{x,k}, c_{y,k})\}_{k}  @f$ where
     * @f$ c_{x,k} @f$ and @f$ c_{y,k} @f$ are the B-splines coefficients:
     *
     * @f$ R(\rho,\theta) = \sum_{k=0}^{N_\rho\times N_{\theta}-1} c_{x, k} B_k{\rho,\theta} @f$,
     *
     * @f$ Z(\rho,\theta) = \sum_{k=0}^{N_\rho\times N_{\theta}-1} c_{y, k} B_k{\rho,\theta} @f$,
     *
     * where @f$ N_\rho\times N_{\theta} @f$ is the number of B-splines.
     *
     * The control points can be obtained by interpolating the mapping on interpolation
     * points (see GrevilleInterpolationPoints or KnotsAsInterpolationPoints).
     * We can also note that the first control points @f$ \{(c_{x,k}, c_{y,k})\}_{k=0}^{N_{\theta}} @f$
     * are equal to the pole @f$ (c_{x,k}, c_{y,k}) = (R_0,Z_0) @f$, @f$ \forall k = 0, ..., N_{\theta}-1 @f$
     * where @f$ R(0,\theta), Z(0,\theta) = (R_0,Z_0) @f$ @f$ \forall \theta @f$.
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
    KOKKOS_INLINE_FUNCTION const Coord<R, Z> control_point(
            Idx<BSplineRho, BSplineTheta> const& el) const
    {
        return Coord<R, Z>(m_R_spline_representation(el), m_Z_spline_representation(el));
    }
};


namespace mapping_detail {
template <
        class R,
        class Z,
        class SplineEvaluator,
        class Rho,
        class Theta,
        class MemorySpace,
        class ExecSpace>
struct MappingAccessibility<
        ExecSpace,
        DiscreteToCartesian<R, Z, SplineEvaluator, Rho, Theta, MemorySpace>>
{
    static constexpr bool value = Kokkos::SpaceAccessibility<ExecSpace, MemorySpace>::accessible;
};

template <class R, class Z, class SplineEvaluator, class Rho, class Theta, class MemorySpace>
struct IsCurvilinear2DMapping<DiscreteToCartesian<R, Z, SplineEvaluator, Rho, Theta, MemorySpace>>
    : std::true_type
{
};

template <class R, class Z, class SplineEvaluator, class Rho, class Theta, class MemorySpace>
struct SingularOPointInvJacobian<
        DiscreteToCartesian<R, Z, SplineEvaluator, Rho, Theta, MemorySpace>> : std::true_type
{
};

} // namespace mapping_detail
