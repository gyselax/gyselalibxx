// SPDX-License-Identifier: MIT
#pragma once

#include "cartesian_to_circular.hpp"
#include "circular_to_cartesian.hpp"
#include "czarny_to_cartesian.hpp"
#include "discrete_to_cartesian.hpp"
#include "indexed_tensor.hpp"
#include "mapping_tools.hpp"
#include "view.hpp"

/**
 * @brief An operator for calculating the inverse of the Jacobian at an O-point.
 * This class is used in CombinedMapping to calculate the inverse of the Jacobian
 * at an O-point when one of the mappings does not allow the evaluation of its
 * Jacobian/inverse Jacobian at the O-point.
 *
 * Specialisations of this class must implement:
 * - A constructor taking the Mapping as an argument
 * - An operator() returning the Jacobian matrix at the O-point.
 *
 * @tparam Mapping The mapping for which the inverse of the Jacobian is calculated.
 * @tparam CoordRTheta The coordinate system in which the inverse of the Jacobian
 *                  is calculated.
 */
template <class Mapping, class CoordRTheta>
class InvJacobianOPoint;

// Pre-declaration of CombinedMapping.
template <class Mapping1, class Mapping2>
class CombinedMapping;

/**
 * A specialisation of InvJacobianOPoint for a combined mapping 
 * @f$ \mathcal{F} \circ \mathcal{G} @f$  where @f$ \mathcal{F} @f$ is a circular
 * mapping from logical to physical, and @f$ \mathcal{G} @f$ is an inverse circular
 * mapping from physical to logical. The combined mapping @f$ \mathcal{F} \circ \mathcal{G} @f$ 
 * therefore maps from a physical domain @f$ (X_{pc}, Y_{pc}) @f$ to a physical domain 
 * @f$ (X, Y) @f$ (this mapping is equivalent to the identity).
 */
template <class X, class Y, class R, class Theta, class Xpc, class Ypc>
class InvJacobianOPoint<
        CombinedMapping<
                CircularToCartesian<R, Theta, X, Y>,
                CartesianToCircular<Xpc, Ypc, R, Theta>>,
        Coord<R, Theta>>
{
    /// The coordinate system in which the inverse of the Jacobian is calculated.
    using CoordRTheta = Coord<R, Theta>;

    /// @brief The covariant form of the first physical coordinate.
    using X_cov = typename X::Dual;
    /// @brief The covariant form of the second physical coordinate.
    using Y_cov = typename Y::Dual;

private:
    using Mapping = CombinedMapping<
            CircularToCartesian<R, Theta, X, Y>,
            CartesianToCircular<Xpc, Ypc, R, Theta>>;

private:
    Mapping m_mapping;

public:
    /**
     * @brief The constructor of InvJacobianOPoint.
     * @param[in] mapping The mapping for which the inverse of the Jacobian is calculated.
     */
    KOKKOS_FUNCTION explicit InvJacobianOPoint(Mapping const& mapping) : m_mapping(mapping) {}

    /**
     * @brief Compute the full inverse Jacobian matrix from a coordinate system (x_pc, y_pc) to a
     * coordinate system (x, y).
     *
     * Here, as @f$ \mathcal{G}^{-1} =  \mathcal{F} @f$, the Jacobian matrix of
     * @f$(\mathcal{F} \circ \mathcal{G}^{-1})^{-1} @f$ is the identity matrix.
     * So, the pseudo-Cartesian Jacobian matrix for a circular mapping is given by :
     * - @f$ (J_{\mathcal{F}}J_{\mathcal{G}}^{-1})^{-1}_{11}(0, \theta) = 1, @f$
     * - @f$ (J_{\mathcal{F}}J_{\mathcal{G}}^{-1})^{-1}_{12}(0, \theta) = 0, @f$
     * - @f$ (J_{\mathcal{F}}J_{\mathcal{G}}^{-1})^{-1}_{21}(0, \theta) = 0, @f$
     * - @f$ (J_{\mathcal{F}}J_{\mathcal{G}}^{-1})^{-1}_{22}(0, \theta) = 1. @f$
     *
     * @result The matrix evaluated at the central point.
     */
    KOKKOS_INLINE_FUNCTION DTensor<VectorIndexSet<Xpc, Ypc>, VectorIndexSet<X_cov, Y_cov>>
    operator()() const
    {
        DTensor<VectorIndexSet<Xpc, Ypc>, VectorIndexSet<X_cov, Y_cov>> J;
        ddcHelper::get<Xpc, X_cov>(J) = 1.;
        ddcHelper::get<Xpc, Y_cov>(J) = 0.;
        ddcHelper::get<Ypc, X_cov>(J) = 0.;
        ddcHelper::get<Ypc, Y_cov>(J) = 1.;
        return J;
    }
};

/**
 * A specialisation of InvJacobianOPoint for a combined mapping
 * @f$ \mathcal{F} \circ \mathcal{G} @f$ where @f$ \mathcal{F} @f$ is a Czarny mapping 
 * from logical to physical, and @f$ \mathcal{G} @f$ is an inverse circular mapping 
 * from physical to logical. The combined mapping @f$ \mathcal{F} \circ \mathcal{G} @f$ 
 * therefore maps from a physical domain @f$ (X_{pc}, Y_{pc}) @f$ to a physical domain
 * @f$ (X, Y) @f$.
 */
template <class X, class Y, class R, class Theta, class Xpc, class Ypc>
class InvJacobianOPoint<
        CombinedMapping<CzarnyToCartesian<R, Theta, X, Y>, CartesianToCircular<Xpc, Ypc, R, Theta>>,
        Coord<R, Theta>>
{
    /// The coordinate system in which the inverse of the Jacobian is calculated.
    using CoordRTheta = Coord<R, Theta>;

    /// @brief The covariant form of the first physical coordinate.
    using X_cov = typename X::Dual;
    /// @brief The covariant form of the second physical coordinate.
    using Y_cov = typename Y::Dual;

private:
    using CzarnyToCart = CzarnyToCartesian<R, Theta, X, Y>;
    using Mapping = CombinedMapping<CzarnyToCart, CartesianToCircular<Xpc, Ypc, R, Theta>>;

private:
    Mapping m_mapping;

public:
    /**
     * @brief The constructor of InvJacobianOPoint.
     * @param[in] mapping The mapping for which the inverse of the Jacobian is calculated.
     */
    KOKKOS_FUNCTION explicit InvJacobianOPoint(Mapping const& mapping) : m_mapping(mapping) {}

    /**
     * @brief Compute the full inverse Jacobian matrix from a coordinate system (x_pc, y_pc) to a
     * coordinate system (x, y) at the central point.
     *
     * The pseudo-Cartesian Jacobian matrix for a Czarny mapping is given by :
     * - @f$ (J_{\mathcal{F}}J_{\mathcal{G}}^{-1})^{-1}_{11}(0, \theta) = - \sqrt{1 + \varepsilon^2}, @f$
     * - @f$ (J_{\mathcal{F}}J_{\mathcal{G}}^{-1})^{-1}_{12}(0, \theta) = 0, @f$
     * - @f$ (J_{\mathcal{F}}J_{\mathcal{G}}^{-1})^{-1}_{21}(0, \theta) = 0, @f$
     * - @f$ (J_{\mathcal{F}}J_{\mathcal{G}}^{-1})^{-1}_{22}(0, \theta) = \frac{2 - \sqrt{1 + \varepsilon^2}}{e \xi}. @f$
     *
     * @result The matrix evaluated at the central point.
     */
    KOKKOS_INLINE_FUNCTION DTensor<VectorIndexSet<Xpc, Ypc>, VectorIndexSet<X_cov, Y_cov>>
    operator()() const
    {
        const double epsilon = m_mapping.template get<CzarnyToCart>().epsilon();
        const double e = m_mapping.template get<CzarnyToCart>().e();
        const double xi = Kokkos::sqrt(1. / (1. - epsilon * epsilon * 0.25));
        const double sqrt_eps_2 = Kokkos::sqrt(1. + epsilon * epsilon);

        DTensor<VectorIndexSet<Xpc, Ypc>, VectorIndexSet<X_cov, Y_cov>> J;
        ddcHelper::get<Xpc, X_cov>(J) = -sqrt_eps_2;
        ddcHelper::get<Xpc, Y_cov>(J) = 0.;
        ddcHelper::get<Ypc, X_cov>(J) = 0.;
        ddcHelper::get<Ypc, Y_cov>(J) = (2 - sqrt_eps_2) / e / xi;
        return J;
    }
};

/**
 * A specialisation of InvJacobianOPoint for a combined mapping @f$ \mathcal{F} \circ \mathcal{G} @f$ 
 * where @f$ \mathcal{F} @f$ is a discrete mapping from logical to physical, and 
 * @f$ \mathcal{G} @f$ is an inverse circular mapping from physical to logical. 
 * The combined mapping @f$ \mathcal{F} \circ \mathcal{G} @f$ therefore maps from a 
 * physical domain @f$ (X_{pc}, Y_{pc}) @f$ to a physical domain @f$ (X, Y) @f$.
 */
template <
        class X,
        class Y,
        class SplineEvaluator,
        class R,
        class Theta,
        class MemorySpace,
        class Xpc,
        class Ypc>
class InvJacobianOPoint<
        CombinedMapping<
                DiscreteToCartesian<X, Y, SplineEvaluator, R, Theta, MemorySpace>,
                CartesianToCircular<Xpc, Ypc, R, Theta>>,
        Coord<R, Theta>>
{
    /// The coordinate system in which the inverse of the Jacobian is calculated.
    using CoordRTheta = Coord<R, Theta>;

    /// @brief The covariant form of the first physical coordinate.
    using X_cov = typename X::Dual;
    /// @brief The covariant form of the second physical coordinate.
    using Y_cov = typename Y::Dual;
    using Xpc_cov = typename Xpc::Dual;
    using Ypc_cov = typename Ypc::Dual;
    using R_cov = typename R::Dual;
    using Theta_cov = typename Theta::Dual;

private:
    using Mapping = CombinedMapping<
            DiscreteToCartesian<X, Y, SplineEvaluator, R, Theta, MemorySpace>,
            CartesianToCircular<Xpc, Ypc, R, Theta>>;

    using IdxRangeR = typename SplineEvaluator::evaluation_domain_type1;
    using IdxRangeTheta = typename SplineEvaluator::evaluation_domain_type2;
    using IdxRangeRTheta = typename SplineEvaluator::evaluation_domain_type;
    using IdxR = typename IdxRangeR::discrete_element_type;
    using IdxTheta = typename IdxRangeTheta::discrete_element_type;

private:
    Mapping m_mapping;

public:
    /**
     * @brief The constructor of InvJacobianOPoint.
     * @param[in] mapping The mapping for which the inverse of the Jacobian is calculated.
     */
    KOKKOS_FUNCTION explicit InvJacobianOPoint(Mapping const& mapping) : m_mapping(mapping) {}

    /**
     * @brief Compute the full inverse Jacobian matrix from a coordinate system (x_pc, y_pc) to a
     * coordinate system (x, y) at the central point.
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
     * So the pseudo-Cartesian Jacobian matrix at the central point,
     * @f$ (J_{\mathcal{F}}J_{\mathcal{G}}^{-1})^{-1}(0, \theta) @f$, is obtained by inversing this matrix. 
     *
     *
     * @result The matrix evaluated at the central point.
     *
     * @see BslAdvection
     * @see AdvectionDomain
     */
    KOKKOS_FUNCTION DTensor<VectorIndexSet<Xpc, Ypc>, VectorIndexSet<X_cov, Y_cov>> operator()()
            const
    {
        DiscreteToCartesian<X, Y, SplineEvaluator, R, Theta, MemorySpace> const& discrete_mapping
                = m_mapping.template get<
                        DiscreteToCartesian<X, Y, SplineEvaluator, R, Theta, MemorySpace>>();
        DTensor<VectorIndexSet<Xpc, Ypc>, VectorIndexSet<X_cov, Y_cov>> J(0.0);
        IdxRangeRTheta idx_range_singular_point = discrete_mapping.idx_range_singular_point();
        // Average the values at (r = 0, theta):
        IdxR ir(idx_range_singular_point.front());
        for (IdxTheta itheta : IdxRangeTheta(idx_range_singular_point)) {
            Coord<R, Theta> coord(ddc::coordinate(ir), ddc::coordinate(itheta));
            DTensor<VectorIndexSet<X, Y>, VectorIndexSet<R_cov, Theta_cov>> J_first_order
                    = discrete_mapping.first_order_jacobian_matrix_r_rtheta(coord);

            double th = ddc::get<Theta>(coord);
            DTensor<VectorIndexSet<R, Theta>, VectorIndexSet<Xpc_cov, Ypc_cov>> J_circ_r_rtheta;
            ddcHelper::get<R, Xpc_cov>(J_circ_r_rtheta) = Kokkos::cos(th);
            ddcHelper::get<R, Ypc_cov>(J_circ_r_rtheta) = Kokkos::sin(th);
            ddcHelper::get<Theta, Xpc_cov>(J_circ_r_rtheta) = -Kokkos::sin(th);
            ddcHelper::get<Theta, Ypc_cov>(J_circ_r_rtheta) = Kokkos::cos(th);

            Tensor J_theta = inverse(
                    tensor_mul(index<'i', 'j'>(J_first_order), index<'j', 'k'>(J_circ_r_rtheta)));

            J += J_theta;
        }

        int const theta_size = idx_range_singular_point.size();
        J /= theta_size;
        return J;
    }
};
