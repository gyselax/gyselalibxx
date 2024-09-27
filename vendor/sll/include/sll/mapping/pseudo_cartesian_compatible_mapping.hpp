// SPDX-License-Identifier: MIT
#pragma once
#include <sll/view.hpp>

/**
 * @brief An operator to calculate the Jacobian matrix of the mapping from the physical domain to the pseudo-Cartesian domain at the central point.
 * All operators which can calculate the pseudo-Cartesian Jacobian at the central point should inherit from this class.
 *
 * Some mapping can be difficult to invert especially at the central point.
 * In the case of a non-analytically invertible mapping, we can work in another domain called the pseudo-Cartesian domain.
 * In this domain, it is easier to inverse the Jacobian matrix. The idea is detailed in Edoardo Zoni's article
 * (*Solving hyperbolic-elliptic problems on singular mapped disk-like domains with the method of characteristics and
 * spline finite elements*, https://doi.org/10.1016/j.jcp.2019.108889)
 *
 * The mapping maps from the logical domain to the physical domain @f$ \mathcal{F}: (r,\theta) \mapsto (x, y)@f$.
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
 * The pseudo-Cartesian Jacobian matrix at the central point,
 * @f$ (J_{\mathcal{F}}J_{\mathcal{G}}^{-1})^{-1}(0, \theta) @f$, is obtained by inversing this matrix.
 *
 * For analytical invertible mappings, this technique is not necessary, but it is still possible to define and use the
 * matrix in order to handle problems in a general manner.
 */
class PseudoCartesianCompatibleMapping
{
public:
    /// The type of the pseudo-Cartesian Jacobian matrix and its inverse
    using Matrix_2x2 = std::array<std::array<double, 2>, 2>;

public:
    /**
     * @brief  Compute the full Jacobian matrix from the mapping to the pseudo-Cartesian mapping at the central point.
     *
     * @param[out] matrix
     *      The pseudo-Cartesian matrix evaluated at the central point.
     *
     * @see DiscreteToCartesian
     * @see BslAdvection
     * @see AdvectionDomain
     */
    virtual KOKKOS_FUNCTION void to_pseudo_cartesian_jacobian_center_matrix(
            Matrix_2x2& matrix) const = 0;

    /**
     * @brief Compute the (1,1) coefficient of the pseudo-Cartesian Jacobian matrix at the central point.
     *
     * @f$ (J_{\mathcal{F}}J_{\mathcal{G}}^{-1})^{-1}_{11}(0, \theta) =
     *     \partial_r x (0, \theta) \cos(\theta) - \partial_{r \theta} x (0, \theta) \sin(\theta) @f$
     *
     * @param[in] grid
     *      The domain where the mapping is defined.
     *
     * @return A double with the (1,1) coefficient of the pseudo-Cartesian Jacobian matrix at the central point.
     *
     */
    virtual KOKKOS_FUNCTION double to_pseudo_cartesian_jacobian_11_center() const = 0;

    /**
     * @brief Compute the (1,2) coefficient of the pseudo-Cartesian Jacobian matrix at the central point.
     *
     * @f$ (J_{\mathcal{F}}J_{\mathcal{G}}^{-1})^{-1}_{12}(0, \theta) =
     *     \partial_r x (0, \theta) \sin(\theta) + \partial_{r \theta} x (0, \theta) \cos(\theta) @f$
     *
     * @param[in] grid
     *      The domain where the mapping is defined.
     *
     * @return A double with the (1,2) coefficient of the pseudo-Cartesian Jacobian matrix at the central point.
     *
     *
     */
    virtual KOKKOS_FUNCTION double to_pseudo_cartesian_jacobian_12_center() const = 0;

    /**
     * @brief Compute the (2,1) coefficient of the pseudo-Cartesian Jacobian matrix at the central point.
     *
     * @f$ (J_{\mathcal{F}}J_{\mathcal{G}}^{-1})^{-1}_{21}(0, \theta) =
     *     \partial_r y (0, \theta) \cos(\theta) - \partial_{r \theta} y (0, \theta) \sin(\theta)@f$
     *
     * @param[in] grid
     *      The domain where the mapping is defined.
     *
     * @return A double with the (2,1) coefficient of the pseudo-Cartesian Jacobian matrix at the central point.
     *
     *
     */
    virtual KOKKOS_FUNCTION double to_pseudo_cartesian_jacobian_21_center() const = 0;

    /**
     * @brief Compute the (2,2) coefficient of the pseudo-Cartesian Jacobian matrix at the central point.
     *
     * @f$ (J_{\mathcal{F}}J_{\mathcal{G}}^{-1})^{-1}_{22}(0, \theta) =
     *     \partial_r y (0, \theta) \sin(\theta) + \partial_{r \theta} y (0, \theta) \cos(\theta) @f$
     *
     * @param[in] grid
     *      The domain where the mapping is defined.
     *
     * @return A double with the (2,2) coefficient of the pseudo-Cartesian Jacobian matrix at the central point.
     *
     * 
     */
    virtual KOKKOS_FUNCTION double to_pseudo_cartesian_jacobian_22_center() const = 0;
};
