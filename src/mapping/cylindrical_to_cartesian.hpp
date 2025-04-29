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
template <class X, class Y, class Z, class R, class Zeta>
class CartesianToCylindrical;

/**
 * @brief A class for describing the cylindrical 3D mapping.
 *
 * The mapping @f$ (R, Z, \zeta)\mapsto (x,y,z) @f$ is defined as follow :
 *
 * @f$ x(R, Z, \zeta) = R \cos(\zeta), @f$
 *
 * @f$ y(R, Z, \zeta) = R \sin(\zeta). @f$
 *
 * @f$ z(R, Z, \zeta) = z. @f$
 *
 * It and its Jacobian matrix are invertible everywhere except for @f$ R = 0 @f$.
 *
 * The Jacobian matrix coefficients are defined as follow
 *
 * @f$ J^x_{\;R}(R, Z, \zeta) = \cos(\zeta) @f$
 *
 * @f$ J^x_{\;Z}(R, Z, \zeta) = 0 @f$
 *
 * @f$ J^x_{\;\zeta}(R, Z, \zeta) = - R \sin(\zeta) @f$
 *
 * @f$ J^y_{\;R}(R, Z, \zeta) = \sin(\zeta) @f$
 *
 * @f$ J^y_{\;Z}(R, Z, \zeta) = 0 @f$
 *
 * @f$ J^y_{\;\zeta}(R, Z, \zeta) = R \cos(\zeta) @f$
 *
 * @f$ J^z_{\;R}(R, Z, \zeta) = 0 @f$
 *
 * @f$ J^z_{\;Z}(R, Z, \zeta) = 1 @f$
 *
 * @f$ J^z_{\;\zeta}(R, Z, \zeta) = 0 @f$
 *
 * and the matrix determinant: @f$ \det(J) = -R @f$.
 *
 */
template <class R, class Z, class Zeta, class X, class Y>
class CylindricalToCartesian
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
    using CoordArg = Coord<R, Z, Zeta>;
    /// The type of the result of the function described by this mapping
    using CoordResult = Coord<X, Y, Z>;

    /// @brief The covariant form of the first Cartesian coordinate.
    using X_cov = typename X::Dual;
    /// @brief The covariant form of the second Cartesian coordinate.
    using Y_cov = typename Y::Dual;
    /// @brief The covariant form of the third Cartesian coordinate and the longitudinal cylindrical coordinate.
    using Z_cov = typename Z::Dual;
    /// @brief The covariant form of the radial cylindrical coordinate.
    using R_cov = typename R::Dual;
    /// @brief The covariant form of the angular cylindrical coordinate.
    using Zeta_cov = typename Zeta::Dual;

public:
    CylindricalToCartesian() = default;

    /**
     * @brief Instantiate a CylindricalToCartesian from another CylindricalToCartesian (lvalue).
     *
     * @param[in] other
     * 		CylindricalToCartesian mapping used to instantiate the new one.
     */
    KOKKOS_FUNCTION CylindricalToCartesian(CylindricalToCartesian const& other) {}

    /**
     * @brief Instantiate a CylindricalToCartesian from another temporary CylindricalToCartesian (rvalue).
     *
     * @param[in] x
     * 		CylindricalToCartesian mapping used to instantiate the new one.
     */
    CylindricalToCartesian(CylindricalToCartesian&& x) = default;

    ~CylindricalToCartesian() = default;

    /**
     * @brief Assign a CylindricalToCartesian from another CylindricalToCartesian (lvalue).
     *
     * @param[in] x
     * 		CylindricalToCartesian mapping used to assign.
     *
     * @return The CylindricalToCartesian assigned.
     */
    CylindricalToCartesian& operator=(CylindricalToCartesian const& x) = default;

    /**
     * @brief Assign a CylindricalToCartesian from another temporary CylindricalToCartesian (rvalue).
     *
     * @param[in] x
     * 		CylindricalToCartesian mapping used to assign.
     *
     * @return The CylindricalToCartesian assigned.
     */
    CylindricalToCartesian& operator=(CylindricalToCartesian&& x) = default;

    /**
     * @brief Convert the @f$ (r, \zeta) @f$ coordinate to the equivalent (x,y) coordinate.
     *
     * @param[in] coord The coordinate to be converted.
     *
     * @return The equivalent coordinate.
     */
    KOKKOS_FUNCTION CoordResult operator()(CoordArg const& coord) const
    {
        const double r = ddc::get<R>(coord);
        const double zeta = ddc::get<Zeta>(coord);
        Coord<X> x(r * Kokkos::cos(zeta));
        Coord<Y> y(r * Kokkos::sin(zeta));
        return CoordResult(x, y, ddc::select<Z>(coord));
    }

    /**
     * @brief Compute the Jacobian, the determinant of the Jacobian matrix of the mapping.
     *
     * @param[in] coord
     *          The coordinate where we evaluate the Jacobian.
     *
     * @return A double with the value of the determinant of the Jacobian matrix.
     */
    KOKKOS_FUNCTION double jacobian(CoordArg const& coord) const
    {
        double r = ddc::get<R>(coord);
        return -r;
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
    KOKKOS_FUNCTION DTensor<VectorIndexSet<X, Y, Z>, VectorIndexSet<R_cov, Z_cov, Zeta_cov>>
    jacobian_matrix(CoordArg const& coord) const
    {
        const double r = ddc::get<R>(coord);
        const double zeta = ddc::get<Zeta>(coord);
        DTensor<VectorIndexSet<X, Y, Z>, VectorIndexSet<R_cov, Z_cov, Zeta_cov>> jacobian_matrix;
        ddcHelper::get<X, R_cov>(jacobian_matrix) = Kokkos::cos(zeta);
        ddcHelper::get<X, Z_cov>(jacobian_matrix) = 0;
        ddcHelper::get<X, Zeta_cov>(jacobian_matrix) = -r * Kokkos::sin(zeta);
        ddcHelper::get<Y, R_cov>(jacobian_matrix) = Kokkos::sin(zeta);
        ddcHelper::get<Y, Z_cov>(jacobian_matrix) = 0;
        ddcHelper::get<Y, Zeta_cov>(jacobian_matrix) = r * Kokkos::cos(zeta);
        ddcHelper::get<Z, R_cov>(jacobian_matrix) = 0;
        ddcHelper::get<Z, Z_cov>(jacobian_matrix) = 1;
        ddcHelper::get<Z, Zeta_cov>(jacobian_matrix) = 0;
        return jacobian_matrix;
    }

    /**
     * @brief Compute the (i,j) coefficient of the Jacobian matrix.
     * 
     * @param[in] coord
     *              The coordinate where we evaluate the Jacobian matrix.
     *
     * @return A double with the value of the (i,j) coefficient of the Jacobian matrix.
     */
    template <class IndexTag1, class IndexTag2>
    KOKKOS_FUNCTION double jacobian_component(CoordArg const& coord) const
    {
        static_assert(ddc::in_tags_v<IndexTag1, VectorIndexSet<X, Y, Z>>);
        static_assert(ddc::in_tags_v<IndexTag2, VectorIndexSet<R_cov, Z_cov, Zeta_cov>>);

        const double zeta = ddc::get<Zeta>(coord);

        if constexpr (std::is_same_v<IndexTag1, X> && std::is_same_v<IndexTag2, R_cov>) {
            //Compute the (1,1) coefficient of the Jacobian matrix, i.e J^x_r.
            return Kokkos::cos(zeta);
        } else if constexpr (std::is_same_v<IndexTag1, X> && std::is_same_v<IndexTag2, Zeta_cov>) {
            //Compute the (1,2) coefficient of the Jacobian matrix, i.e J^x_theta.
            const double r = ddc::get<R>(coord);
            return -r * Kokkos::sin(zeta);
        } else if constexpr (std::is_same_v<IndexTag1, Y> && std::is_same_v<IndexTag2, R_cov>) {
            //Compute the (2,1) coefficient of the Jacobian matrix, i.e J^y_r.
            return Kokkos::sin(zeta);
        } else if constexpr (std::is_same_v<IndexTag1, Y> && std::is_same_v<IndexTag2, Zeta_cov>) {
            //Compute the (2,2) coefficient of the Jacobian matrix, i.e J^y_theta.
            const double r = ddc::get<R>(coord);
            return r * Kokkos::cos(zeta);
        } else if constexpr (std::is_same_v<IndexTag1, Z> && std::is_same_v<IndexTag2, Z_cov>) {
            return 1;
        } else {
            return 0;
        }
    }


    /**
     * @brief Compute full inverse Jacobian matrix.
     *
     * For some computations, we need the complete inverse Jacobian matrix or just the
     * coefficients.
     * The coefficients can be given independently with the function inv_jacobian_component.
     *
     * @param[in] coord
     * 				The coordinate where we evaluate the Jacobian matrix.
     * @return The inverse Jacobian matrix.
     */
    KOKKOS_FUNCTION DTensor<VectorIndexSet<R, Z, Zeta>, VectorIndexSet<X_cov, Y_cov, Z_cov>>
    inv_jacobian_matrix(CoordArg const& coord) const
    {
        const double r = ddc::get<R>(coord);
        const double zeta = ddc::get<Zeta>(coord);
        assert(fabs(r) >= 1e-15);

        DTensor<VectorIndexSet<R, Z, Zeta>, VectorIndexSet<X_cov, Y_cov, Z_cov>> matrix(0);
        ddcHelper::get<R, X_cov>(matrix) = Kokkos::cos(zeta);
        ddcHelper::get<R, Y_cov>(matrix) = Kokkos::sin(zeta);
        ddcHelper::get<Zeta, X_cov>(matrix) = -1 / r * Kokkos::sin(zeta);
        ddcHelper::get<Zeta, Y_cov>(matrix) = 1 / r * Kokkos::cos(zeta);
        ddcHelper::get<Z, Z_cov>(matrix) = 1;
        return matrix;
    }


    /**
     * @brief Compute the (i,j) coefficient of the inverse Jacobian matrix.
     *
     * Be careful because not all mappings are invertible, especially at the centre point.
     *
     * @param[in] coord
     *              The coordinate where we evaluate the inverse Jacobian matrix.
     *
     * @return A double with the value of the (i,j) coefficient of the inverse Jacobian matrix.
     */
    template <class IndexTag1, class IndexTag2>
    KOKKOS_FUNCTION double inv_jacobian_component(CoordArg const& coord) const
    {
        static_assert(ddc::in_tags_v<IndexTag1, VectorIndexSet<R, Z, Zeta>>);
        static_assert(ddc::in_tags_v<IndexTag2, VectorIndexSet<X_cov, Y_cov, Z_cov>>);

        const double zeta = ddc::get<Zeta>(coord);

        if constexpr (std::is_same_v<IndexTag1, R> && std::is_same_v<IndexTag2, X_cov>) {
            //Compute the (1,1) coefficient of the inverse Jacobian matrix.
            return Kokkos::cos(zeta);
        } else if constexpr (std::is_same_v<IndexTag1, R> && std::is_same_v<IndexTag2, Y_cov>) {
            //Compute the (1,2) coefficient of the inverse Jacobian matrix.
            return Kokkos::sin(zeta);
        } else if constexpr (std::is_same_v<IndexTag1, Z> && std::is_same_v<IndexTag2, Z_cov>) {
            return 1;
        } else if constexpr (std::is_same_v<IndexTag1, Z> || std::is_same_v<IndexTag2, Z_cov>) {
            return 0;
        } else {
            const double r = ddc::get<R>(coord);
            assert(fabs(r) >= 1e-15);
            if constexpr (std::is_same_v<IndexTag2, X_cov>) {
                //Compute the (2,1) coefficient of the inverse Jacobian matrix.
                return -1 / r * Kokkos::sin(zeta);
            } else {
                //Compute the (2,2) coefficient of the inverse Jacobian matrix.
                return 1 / r * Kokkos::cos(zeta);
            }
        }
    }


    /**
     * @brief Get the inverse mapping.
     *
     * @return The inverse mapping.
     */
    CartesianToCylindrical<X, Y, Z, R, Zeta> get_inverse_mapping() const
    {
        return CartesianToCylindrical<X, Y, Z, R, Zeta>();
    }
};

namespace mapping_detail {
template <class X, class Y, class Z, class R, class Zeta, class ExecSpace>
struct MappingAccessibility<ExecSpace, CylindricalToCartesian<R, Z, Zeta, X, Y>> : std::true_type
{
};

} // namespace mapping_detail
