// SPDX-License-Identifier: MIT
#pragma once

#include <cassert>
#include <cmath>

#include <ddc/ddc.hpp>

#include "ddc_aliases.hpp"
#include "mapping_tools.hpp"
#include "tensor.hpp"
#include "view.hpp"

// Pre-declaration of analytical inverse
template <class R, class Z, class Zeta, class X, class Y>
class CylindricalToCartesian;

/**
 * @brief A class for describing the cylindrical 3D mapping.
 *
 * The mapping @f$ (x,y,z) \mapsto (r,z, \zeta) @f$ is defined as follow :
 *
 * @f$ R(x,y,z) = \sqrt{x^2+y^2} ,@f$
 *
 * @f$ \zeta(x,y,z) = atan2(\frac{y}{x}) .@f$
 *
 * @f$ Z(x,y,z) = z @f$
 *
 * It and its Jacobian matrix are invertible everywhere except for @f$ R = 0 @f$.
 *
 * The Jacobian matrix coefficients are defined as follow
 *
 * @f$ J^R_{\;x}(x,y,z)  =\frac{x}{\sqrt{x^2+y^2}}  @f$
 *
 * @f$ J^R_{\;y}(x,y,z)  =\frac{y}{\sqrt{x^2+y^2}}  @f$
 *
 * @f$ J^R_{\;z}(x,y,z)  =0  @f$
 *
 * @f$ J^Z_{\;x}(x,y,z)  =0  @f$
 *
 * @f$ J^Z_{\;y}(x,y,z)  =0  @f$
 *
 * @f$ J^Z_{\;z}(x,y,z)  =1  @f$
 *
 * @f$ J^{\zeta}_{\;x}(x,y,z)  =\frac{-y}{x^2+y^2}  @f$
 *
 * @f$ J^{\zeta}_{\;y}(x,y,z)  =\frac{x}{x^2+y^2}  @f$
 *
 * @f$ J^{\zeta}_{\;z}(x,y,z)  =0  @f$
 *
 * and the matrix determinant: @f$ det(J) = -\frac{1}{\sqrt{x^2+y^2}} @f$.
 *
 */
template <class X, class Y, class Z, class R, class Zeta>
class CartesianToCylindrical
{
public:
    /// @brief Indicate the first Cartesian coordinate.
    using cartesian_tag_x = X;
    /// @brief Indicate the second Cartesian coordinate.
    using cartesian_tag_y = Y;
    /// @brief Indicate the second Cartesian coordinate.
    using cartesian_tag_z = Z;
    /// @brief Indicate the radial cylindrical coordinate.
    using cylindrical_tag_R = R;
    /// @brief Indicate the longitudinal cylindrical coordinate.
    using cylindrical_tag_Z = Z;
    /// @brief Indicate the angular cylindrical coordinate.
    using cylindrical_tag_Zeta = Zeta;

    /// The type of the argument of the function described by this mapping
    using CoordArg = Coord<X, Y, Z>;
    /// The type of the result of the function described by this mapping
    using CoordResult = Coord<R, Z, Zeta>;

    /// @brief The covariant form of the first Cartesian coordinate.
    using X_cov = typename X::Dual;
    /// @brief The covariant form of the second Cartesian coordinate.
    using Y_cov = typename Y::Dual;
    /// @brief The covariant form of the third Cartesian coordinate.
    using Z_cov = typename Z::Dual;
    /// @brief The covariant form of the radial cylindrical coordinate.
    using R_cov = typename R::Dual;
    /// @brief The covariant form of the angular cylindrical coordinate.
    using Zeta_cov = typename Zeta::Dual;

public:
    CartesianToCylindrical() = default;

    /**
     * @brief Instantiate a CartesianToCylindrical from another CartesianToCylindrical (lvalue).
     *
     * @param[in] other
     * 		CartesianToCylindrical mapping used to instantiate the new one.
     */
    KOKKOS_FUNCTION CartesianToCylindrical(CartesianToCylindrical const& other) {}

    /**
     * @brief Instantiate a Curvilinear2DToCartesian from another temporary CartesianToCylindrical (rvalue).
     *
     * @param[in] x
     * 		Curvilinear2DToCartesian mapping used to instantiate the new one.
     */
    CartesianToCylindrical(CartesianToCylindrical&& x) = default;

    ~CartesianToCylindrical() = default;

    /**
     * @brief Assign a CartesianToCylindrical from another CartesianToCylindrical (lvalue).
     *
     * @param[in] x
     * 		CartesianToCylindrical mapping used to assign.
     *
     * @return The CartesianToCylindrical assigned.
     */
    CartesianToCylindrical& operator=(CartesianToCylindrical const& x) = default;

    /**
     * @brief Assign a CartesianToCylindrical from another temporary CartesianToCylindrical (rvalue).
     *
     * @param[in] x
     * 		CartesianToCylindrical mapping used to assign.
     *
     * @return The CartesianToCylindrical assigned.
     */
    CartesianToCylindrical& operator=(CartesianToCylindrical&& x) = default;

    /**
     * @brief Convert the coordinate (x,y,z) to the equivalent @f$ (r, z, \zeta) @f$ coordinate.
     *
     * @param[in] coord The coordinate to be converted.
     *
     * @return The equivalent coordinate.
     */
    KOKKOS_FUNCTION CoordResult operator()(CoordArg const& coord) const
    {
        const double x = ddc::get<X>(coord);
        const double y = ddc::get<Y>(coord);
        Coord<R> r(Kokkos::sqrt(x * x + y * y));
        const double zeta_pi_to_pi(Kokkos::atan2(y, x));
        Coord<Zeta> zeta(zeta_pi_to_pi + 2 * M_PI * (zeta_pi_to_pi < 0));

        return CoordResult(r, ddc::select<Z>(coord), zeta);
    }

    /**
     * @brief Compute the Jacobian, the determinant of the Jacobian matrix of the mapping.
     *
     * @param[in] coord
     *          The coordinate where we evaluate the Jacobian.
     *
     * @return A double with the value of the determinant of the Jacobian matrix.
     */
    KOKKOS_FUNCTION double jacobian(CoordArg const& coord)
    {
        const double x = ddc::get<X>(coord);
        const double y = ddc::get<Y>(coord);
        return -1. / Kokkos::sqrt(x * x + y * y);
    }

    /**
     * @brief Compute full Jacobian matrix.
     *
     * For some computations, we need the complete Jacobian matrix or just the
     * coefficients.
     * The coefficients can be given independently with the function jacobian_component.
     *
     * @param[in] coord
     * 				The coordinate where we evaluate the Jacobian matrix.
     * @return The Jacobian matrix.
     */
    KOKKOS_FUNCTION DTensor<VectorIndexSet<R, Z, Zeta>, VectorIndexSet<X_cov, Y_cov, Z_cov>>
    jacobian_matrix(CoordArg const& coord) const
    {
        const double x = ddc::get<X>(coord);
        const double y = ddc::get<Y>(coord);

        DTensor<VectorIndexSet<R, Z, Zeta>, VectorIndexSet<X_cov, Y_cov, Z_cov>> matrix;
        ddcHelper::get<R, X_cov>(matrix) = x / Kokkos::sqrt(x * x + y * y);
        ddcHelper::get<R, Y_cov>(matrix) = y / Kokkos::sqrt(x * x + y * y);
        ddcHelper::get<R, Z_cov>(matrix) = 0;
        ddcHelper::get<Z, X_cov>(matrix) = 0;
        ddcHelper::get<Z, Y_cov>(matrix) = 0;
        ddcHelper::get<Z, Z_cov>(matrix) = 1;
        ddcHelper::get<Zeta, X_cov>(matrix) = -y / (x * x + y * y);
        ddcHelper::get<Zeta, Y_cov>(matrix) = x / (x * x + y * y);
        ddcHelper::get<Zeta, Z_cov>(matrix) = 0;

        return matrix;
    }


    /**
     * @brief Compute the (i,j) coefficient of the Jacobian matrix.
     *
     * Be careful because not all mappings are invertible, especially at the centre point.
     *
     * @param[in] coord
     *              The coordinate where we evaluate the inverse Jacobian matrix.
     *
     * @return A double with the value of the (i,j) coefficient of the inverse Jacobian matrix.
     */
    template <class IndexTag1, class IndexTag2>
    KOKKOS_FUNCTION double jacobian_component(CoordArg const& coord) const
    {
        static_assert(ddc::in_tags_v<IndexTag1, VectorIndexSet<X, Y>>);
        static_assert(ddc::in_tags_v<IndexTag2, VectorIndexSet<R_cov, Zeta_cov>>);
        static_assert(ddc::in_tags_v<IndexTag1, VectorIndexSet<R, Z, Zeta>>);
        static_assert(ddc::in_tags_v<IndexTag2, VectorIndexSet<X_cov, Y_cov, Z_cov>>);

        const double x = ddc::get<X>(coord);
        const double y = ddc::get<Y>(coord);
        if constexpr (std::is_same_v<IndexTag1, X> && std::is_same_v<IndexTag2, R_cov>) {
            // Component (1,1), i.e dx/dr
            return x / Kokkos::sqrt(x * x + y * y);
        } else if constexpr (std::is_same_v<IndexTag1, X> && std::is_same_v<IndexTag2, Zeta_cov>) {
            // Component (1,2), i.e dx/dzeta
            return y / Kokkos::sqrt(x * x + y * y);
        } else if constexpr (std::is_same_v<IndexTag1, Y> && std::is_same_v<IndexTag2, R_cov>) {
            // Component (2,1), i.e dy/dr
            return -y / (x * x + y * y);
        } else if constexpr (std::is_same_v<IndexTag1, Y> && std::is_same_v<IndexTag2, Zeta_cov>) {
            // Component (2,2), i.e dy/dzeta
            return x / (x * x + y * y);
        } else if constexpr (std::is_same_v<IndexTag1, Z> && std::is_same_v<IndexTag2, Z_cov>) {
            return 1;
        } else {
            return 0;
        }
    }

    /**
     * @brief Get the inverse mapping.
     *
     * @return The inverse mapping.
     */
    CylindricalToCartesian<R, Z, Zeta, X, Y> get_inverse_mapping() const
    {
        return CylindricalToCartesian<R, Z, Zeta, X, Y>();
    }
};
