// SPDX-License-Identifier: MIT
#pragma once

#include <sll/view.hpp>

#include "inverse_jacobian_matrix.hpp"
#include "mapping_tools.hpp"

/**
 * @brief A class describing a mapping from a Cartesian domain to a pseudo-Cartesian.
 * This operation is carried out using 2 curvilinear to Cartesian mappings.
 * The Cartesian to pseudo-Cartesian mapping is obtained by combining these 2 mappings.
 */
template <class MappingToCartesian, class MappingToPseudoCartesian>
class CartesianToPseudoCartesian
{
    static_assert(is_curvilinear_2d_mapping_v<MappingToCartesian>);
    static_assert(is_curvilinear_2d_mapping_v<MappingToPseudoCartesian>);
    static_assert(
            std::is_same_v<
                    typename MappingToCartesian::curvilinear_tag_r,
                    typename MappingToPseudoCartesian::curvilinear_tag_r>,
            "The Cartesian and pseudo-Cartesian mappings must share a logical space");
    static_assert(
            std::is_same_v<
                    typename MappingToCartesian::curvilinear_tag_theta,
                    typename MappingToPseudoCartesian::curvilinear_tag_theta>,
            "The Cartesian and pseudo-Cartesian mappings must share a logical space");

public:
    /// The type of a coordinate in the Cartesian geometry.
    using CartesianCoordinate = ddc::Coordinate<
            typename MappingToCartesian::cartesian_tag_x,
            typename MappingToCartesian::cartesian_tag_y>;

    /// The X dimension in the Cartesian geometry.
    using cartesian_tag_x = typename MappingToCartesian::cartesian_tag_x;
    /// The Y dimension in the Cartesian geometry.
    using cartesian_tag_y = typename MappingToCartesian::cartesian_tag_y;

    /// The X dimension in the pseudo-Cartesian geometry.
    using pseudo_cartesian_tag_x = typename MappingToPseudoCartesian::cartesian_tag_x;
    /// The Y dimension in the pseudo-Cartesian geometry.
    using pseudo_cartesian_tag_y = typename MappingToPseudoCartesian::cartesian_tag_y;

    /// The type of the argument of the function described by this mapping
    using CoordArg = ddc::Coordinate<cartesian_tag_x, cartesian_tag_y>;
    /// The type of the result of the function described by this mapping
    using CoordResult = ddc::Coordinate<pseudo_cartesian_tag_x, pseudo_cartesian_tag_y>;

private:
    using R = typename MappingToCartesian::curvilinear_tag_r;
    using Theta = typename MappingToCartesian::curvilinear_tag_theta;
    using CoordRTheta = ddc::Coordinate<R, Theta>;

    static_assert(is_mapping_v<MappingToCartesian>);
    static_assert(is_mapping_v<MappingToPseudoCartesian>);
    static_assert(has_2d_jacobian_v<MappingToCartesian, CoordRTheta>);
    static_assert(has_2d_jacobian_v<MappingToPseudoCartesian, CoordRTheta>);

private:
    MappingToCartesian m_mapping_to_cartesian;
    MappingToPseudoCartesian m_mapping_to_pseudo_cartesian;
    double m_epsilon;

public:
    /**
     * @brief Construct an instance of the class.
     * @param[in] mapping_to_cartesian A curvilinear to Cartesian mapping.
     * @param[in] mapping_to_pseudo_cartesian A curvilinear to Cartesian mapping where the
     *          role of the Cartesian geometry is fulfilled by the pseudo-Cartesian geometry.
     * @param[in] epsilon The parameter @f$ \varepsilon @f$ which determines when a point is
     *          close enough to the central O-point for linearization to be required when
     *          calculating the Jacobian.
     */
    CartesianToPseudoCartesian(
            MappingToPseudoCartesian mapping_to_pseudo_cartesian,
            MappingToCartesian mapping_to_cartesian,
            double epsilon)
        : m_mapping_to_cartesian(mapping_to_cartesian)
        , m_mapping_to_pseudo_cartesian(mapping_to_pseudo_cartesian)
        , m_epsilon(epsilon)
    {
    }

    /**
     * @brief Compute the full Jacobian matrix.
     * This is calculated by combining the Jacobian matrices of the 2 curvilinear mappings.
     * If a point is within @f$ \varepsilon @f$ of the O-point then a linearisation is carried
     * out between the Jacobian matrix at @f$ \varepsilon @f$ and the matrix at the O-point.
     *
     * @f$ J_{(x,y)->(x_pC, y_pC)} = J_{(r, \theta)->(x_pC,y_pC)}J_{(x, y)->(r, \theta)} @f$
     * @f$ = J_{(r, \theta)->(x_pC,y_pC)} [J_{(r, \theta)->(x, y)}]^{-1} @f$
     *
     * @param[in] coord_rtheta The coordinate where we evaluate the Jacobian matrix.
     * @param[out] J The Jacobian matrix returned.
     */
    KOKKOS_FUNCTION void jacobian_matrix(CoordRTheta const& coord_rtheta, Matrix_2x2& J) const
    {
        double r = ddc::get<R>(coord_rtheta);

        Matrix_2x2 inv_jacobian_from_cartesian;
        Matrix_2x2 jacobian_from_pseudo_cartesian;
        InverseJacobianMatrix inv_jacobian_cartesian(m_mapping_to_cartesian);
        if (r < m_epsilon) {
            Matrix_2x2 J_0;
            m_mapping_to_cartesian.to_pseudo_cartesian_jacobian_center_matrix(J_0);

            CoordRTheta coord_eps(m_epsilon, ddc::get<Theta>(coord_rtheta));

            inv_jacobian_from_cartesian = inv_jacobian_cartesian(coord_eps);
            m_mapping_to_pseudo_cartesian
                    .jacobian_matrix(coord_eps, jacobian_from_pseudo_cartesian);
            Matrix_2x2 J_eps = mat_mul(jacobian_from_pseudo_cartesian, inv_jacobian_from_cartesian);

            J[0][0] = (1 - r / m_epsilon) * J_0[0][0] + r / m_epsilon * J_eps[0][0];
            J[0][1] = (1 - r / m_epsilon) * J_0[0][1] + r / m_epsilon * J_eps[0][1];
            J[1][0] = (1 - r / m_epsilon) * J_0[1][0] + r / m_epsilon * J_eps[1][0];
            J[1][1] = (1 - r / m_epsilon) * J_0[1][1] + r / m_epsilon * J_eps[1][1];
        } else {
            inv_jacobian_from_cartesian = inv_jacobian_cartesian(coord_rtheta);
            m_mapping_to_pseudo_cartesian
                    .jacobian_matrix(coord_rtheta, jacobian_from_pseudo_cartesian);
            J = mat_mul(jacobian_from_pseudo_cartesian, inv_jacobian_from_cartesian);
        }
    }

    /**
     * @brief Compute the determinant of the Jacobian matrix.
     * @see jacobian_matrix
     * @param[in] coord_rtheta The coordinate where we evaluate the Jacobian matrix.
     * @returns The determinant of the Jacobian matrix.
     */
    KOKKOS_FUNCTION double jacobian(CoordRTheta const& coord_rtheta) const
    {
        Matrix_2x2 J;
        jacobian_matrix(coord_rtheta, J);
        return determinant(J);
    }

    /**
     * @brief Compute the (1,1) coefficient of the Jacobian matrix.
     * @see jacobian_matrix
     * @param[in] coord_rtheta The coordinate where we evaluate the Jacobian matrix.
     * @return The (1,1) coefficient of the Jacobian matrix.
     */
    KOKKOS_FUNCTION double jacobian_11(CoordRTheta const& coord_rtheta) const
    {
        Matrix_2x2 J;
        jacobian_matrix(coord_rtheta, J);
        return J[0][0];
    }

    /**
     * @brief Compute the (1,2) coefficient of the Jacobian matrix.
     * @see jacobian_matrix
     * @param[in] coord_rtheta The coordinate where we evaluate the Jacobian matrix.
     * @return The (1,2) coefficient of the Jacobian matrix.
     */
    KOKKOS_FUNCTION double jacobian_12(CoordRTheta const& coord_rtheta) const
    {
        Matrix_2x2 J;
        jacobian_matrix(coord_rtheta, J);
        return J[0][1];
    }

    /**
     * @brief Compute the (2,1) coefficient of the Jacobian matrix.
     * @see jacobian_matrix
     * @param[in] coord_rtheta The coordinate where we evaluate the Jacobian matrix.
     * @return The (2,1) coefficient of the Jacobian matrix.
     */
    KOKKOS_FUNCTION double jacobian_21(CoordRTheta const& coord_rtheta) const
    {
        Matrix_2x2 J;
        jacobian_matrix(coord_rtheta, J);
        return J[1][0];
    }

    /**
     * @brief Compute the (2,2) coefficient of the Jacobian matrix.
     * @see jacobian_matrix
     * @param[in] coord_rtheta The coordinate where we evaluate the Jacobian matrix.
     * @return The (2,2) coefficient of the Jacobian matrix.
     */
    KOKKOS_FUNCTION double jacobian_22(CoordRTheta const& coord_rtheta) const
    {
        Matrix_2x2 J;
        jacobian_matrix(coord_rtheta, J);
        return J[1][1];
    }

    /**
     * @brief Compute the full inverse of the Jacobian matrix.
     * This is calculated by combining the Jacobian matrices of the 2 curvilinear mappings.
     * If a point is within @f$ \varepsilon @f$ of the O-point then a linearisation is carried
     * out between the Jacobian matrix at @f$ \varepsilon @f$ and the matrix at the O-point.
     *
     * @f$ [J_{(x,y)->(x_pC, y_pC)}]^{-1} = [J_{(r, \theta)->(x_pC, y_pC)} J_{(x, y)->(r, \theta)}]^{-1} @f$
     * @f$ = [J_{(r, \theta)->(x_pC, y_pC)} [J_{(x, y)->(r, \theta)}]^{-1}]^{-1} @f$
     * @f$ = J_{(x, y)->(r, \theta)}  [J_{(r, \theta)->(x_pC, y_pC)}]^{-1} @f$
     *
     * @param[in] coord_rtheta The coordinate where we evaluate the Jacobian matrix.
     * @param[out] J The Jacobian matrix returned.
     */
    KOKKOS_FUNCTION void inv_jacobian_matrix(CoordRTheta const& coord_rtheta, Matrix_2x2& J) const
    {
        double r = ddc::get<R>(coord_rtheta);

        Matrix_2x2 jacobian_from_cartesian;
        Matrix_2x2 inv_jacobian_from_pseudo_cartesian;
        InverseJacobianMatrix inv_jacobian_pseudo_cartesian(m_mapping_to_pseudo_cartesian);
        if (r < m_epsilon) {
            Matrix_2x2 J_0;
            m_mapping_to_cartesian.to_pseudo_cartesian_jacobian_center_matrix(J_0);
            Matrix_2x2 J_0_inv = inverse(J_0);

            CoordRTheta coord_eps(m_epsilon, ddc::get<Theta>(coord_rtheta));

            m_mapping_to_cartesian.jacobian_matrix(coord_eps, jacobian_from_cartesian);
            inv_jacobian_from_pseudo_cartesian = inv_jacobian_pseudo_cartesian(coord_eps);
            Matrix_2x2 J_eps = mat_mul(inv_jacobian_from_pseudo_cartesian, jacobian_from_cartesian);

            J[0][0] = (1 - r / m_epsilon) * J_0_inv[0][0] + r / m_epsilon * J_eps[0][0];
            J[0][1] = (1 - r / m_epsilon) * J_0_inv[0][1] + r / m_epsilon * J_eps[0][1];
            J[1][0] = (1 - r / m_epsilon) * J_0_inv[1][0] + r / m_epsilon * J_eps[1][0];
            J[1][1] = (1 - r / m_epsilon) * J_0_inv[1][1] + r / m_epsilon * J_eps[1][1];
        } else {
            m_mapping_to_cartesian.jacobian_matrix(coord_rtheta, jacobian_from_cartesian);
            inv_jacobian_from_pseudo_cartesian = inv_jacobian_pseudo_cartesian(coord_rtheta);
            J = mat_mul(jacobian_from_cartesian, inv_jacobian_from_pseudo_cartesian);
        }
    }

    /**
     * @brief Compute the (1,1) coefficient of the inverse of the Jacobian matrix.
     * @see inv_jacobian_matrix
     * @param[in] coord_rtheta The coordinate where we evaluate the inverse of the Jacobian matrix.
     * @return The (1,1) coefficient of the inverse of the Jacobian matrix.
     */
    KOKKOS_FUNCTION double inv_jacobian_11(CoordRTheta const& coord_rtheta) const
    {
        Matrix_2x2 J;
        inv_jacobian_matrix(coord_rtheta, J);
        return J[0][0];
    }

    /**
     * @brief Compute the (1,2) coefficient of the inverse of the Jacobian matrix.
     * @see inv_jacobian_matrix
     * @param[in] coord_rtheta The coordinate where we evaluate the inverse of the Jacobian matrix.
     * @return The (1,2) coefficient of the inverse of the Jacobian matrix.
     */
    KOKKOS_FUNCTION double inv_jacobian_12(CoordRTheta const& coord_rtheta) const
    {
        Matrix_2x2 J;
        inv_jacobian_matrix(coord_rtheta, J);
        return J[0][1];
    }

    /**
     * @brief Compute the (2,1) coefficient of the inverse of the Jacobian matrix.
     * @see inv_jacobian_matrix
     * @param[in] coord_rtheta The coordinate where we evaluate the inverse of the Jacobian matrix.
     * @return The (2,1) coefficient of the inverse of the Jacobian matrix.
     */
    KOKKOS_FUNCTION double inv_jacobian_21(CoordRTheta const& coord_rtheta) const
    {
        Matrix_2x2 J;
        inv_jacobian_matrix(coord_rtheta, J);
        return J[1][0];
    }

    /**
     * @brief Compute the (2,2) coefficient of the inverse of the Jacobian matrix.
     * @see inv_jacobian_matrix
     * @param[in] coord_rtheta The coordinate where we evaluate the inverse of the Jacobian matrix.
     * @return The (2,2) coefficient of the inverse of the Jacobian matrix.
     */
    KOKKOS_FUNCTION double inv_jacobian_22(CoordRTheta const& coord_rtheta) const
    {
        Matrix_2x2 J;
        inv_jacobian_matrix(coord_rtheta, J);
        return J[1][1];
    }
};

namespace mapping_detail {
template <class ExecSpace, class MappingToPseudoCartesian, class MappingToCartesian>
struct MappingAccessibility<
        ExecSpace,
        CartesianToPseudoCartesian<MappingToPseudoCartesian, MappingToCartesian>>
{
    static bool constexpr value = MappingAccessibility<ExecSpace, MappingToPseudoCartesian>::value
                                  && MappingAccessibility<ExecSpace, MappingToCartesian>::value;
};

} // namespace mapping_detail
