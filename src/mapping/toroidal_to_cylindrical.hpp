// SPDX-License-Identifier: MIT
#pragma once

/**
 * @brief A class describing a coordinate change from a toroidal system of coordinates
 * to a cylindrical system of coordinates. The toroidal coordinates are described by a
 * polar plane @f$(\rho, \theta) \f$ and a perpendicular dimension @f$ \phi @f$. The
 * cylindrical coordinates are @f$ (R, Z, \zeta) @f$. @f$ (R, Z) @f$ describe a Cartesian
 * slice. @f$(\rho, \theta) \f$ are therefore defined from this slice with a 2D
 * coordinate change operator.
 *
 * @tparam Curvilinear2DToCartesian An operator describing the coordinate change from
 *      @f$(\rho, \theta) \f$ to @f$ (R, Z) @f$.
 * @tparam Zeta The angle of the cylindrical coordinates.
 * @tparam Phi The toroidal component of the toroidal coordinates.
 */
template <class Curvilinear2DToCartesian, class Zeta, class Phi>
class ToroidalToCylindrical
{
    static_assert(Zeta::IS_CONTRAVARIANT);
    static_assert(Phi::IS_CONTRAVARIANT);

public:
    /// @brief Indicate the first physical coordinate.
    using cylindrical_tag_R = typename Curvilinear2DToCartesian::cartesian_tag_x;
    /// @brief Indicate the second physical coordinate.
    using cylindrical_tag_Z = typename Curvilinear2DToCartesian::cartesian_tag_y;
    /// @brief Indicate the third physical coordinate.
    using cylindrical_tag_Zeta = Zeta;

    /// @brief Indicate the first logical coordinate.
    using toroidal_tag_rho = typename Curvilinear2DToCartesian::curvilinear_tag_r;
    /// @brief Indicate the second logical coordinate.
    using toroidal_tag_theta = typename Curvilinear2DToCartesian::curvilinear_tag_theta;
    /// @brief Indicate the third logical coordinate.
    using toroidal_tag_phi = Phi;

    /// The type of the argument of the function described by this mapping.
    using CoordArg = Coord<Rho, Theta, Zeta>;
    /// The type of the result of the function described by this mapping.
    using CoordResult = Coord<R, Z, Phi>;

    /// The type of the coordinate change on the polar plane.
    using PolarToCartesian = Curvilinear2DToCartesian;

private:
    using R = cylindrical_tag_R;
    using R_cov = typename R::Dual;
    using Z = cylindrical_tag_Z;
    using Z_cov = typename Z::Dual;
    using Zeta_cov = typename Zeta::Dual;
    using Rho = toroidal_tag_rho;
    using Rho_cov = typename Rho::Dual;
    using Theta = toroidal_tag_theta;
    using Theta_cov = typename Theta::Dual;
    using Phi_cov = typename Phi::Dual;

    using CoordArg2D = Coord<Rho, Theta>;

private:
    Curvilinear2DToCartesian m_mapping_2d;

public:
    explicit ToroidalToCylindrical(Curvilinear2DToCartesian const& mapping_2d)
        : m_mapping_2d(mapping_2d)
    {
    }

    /**
     * @brief Convert the @f$ (\rho, \theta, \phi) @f$ coordinate to the equivalent
     * @f$ (R, Z, \zeta) @f$ coordinate.
     *
     * @param[in] coord The coordinate to be converted.
     *
     * @return The equivalent coordinate.
     */
    KOKKOS_FUNCTION CoordResult operator()(CoordArg const& coord) const
    {
        return CoordResult(m_mapping_2d(CoordArg2D(coord)), Coord<Phi>(-ddc::get<Zeta>(coord)));
    }


    KOKKOS_FUNCTION double jacobian(CoordArg const& coord) const
    {
        return -m_mapping_2d.jacobian(CoordArg2D(coord));
    }

    KOKKOS_FUNCTION DTensor<VectorIndexSet<R, Z, Zeta>, VectorIndexSet<Rho_cov, Theta_cov, Phi_cov>>
    jacobian_matrix(CoordArg const& coord) const
    {
        DTensor<VectorIndexSet<R, Z, Zeta>, VectorIndexSet<Rho_cov, Theta_cov, Zeta_cov>> J_3d(0);
        DTensor<VectorIndexSet<R, Z>, VectorIndexSet<Rho_cov, Theta_cov>> J_2d
                = m_mapping_2d.jacobian_matrix(CoordArg2D(coord));
        ddcHelper::get<R, Rho_cov>(J_3d) = ddcHelper::get<R, Rho_cov>(J_2d);
        ddcHelper::get<R, Theta_cov>(J_3d) = ddcHelper::get<R, Theta_cov>(J_2d);
        ddcHelper::get<Z, Rho_cov>(J_3d) = ddcHelper::get<Z, Rho_cov>(J_2d);
        ddcHelper::get<Z, Theta_cov>(J_3d) = ddcHelper::get<Z, Theta_cov>(J_2d);
        ddcHelper::get<Zeta, Phi_cov>(J_3d) = -1;
        return J_3d;
    }

    template <class IndexTag1, class IndexTag2>
    KOKKOS_INLINE_FUNCTION double jacobian_component(CoordArg const& coord) const
    {
        static_assert(ddc::in_tags_v<IndexTag1, ddc::detail::TypeSeq<R, Z, Zeta>>);
        static_assert(ddc::in_tags_v<IndexTag2, ddc::detail::TypeSeq<Rho_cov, Theta_cov, Phi_cov>>);

        if constexpr (std::is_same_v<IndexTag1, Zeta> && std::is_same_v<IndexTag2, Phi_cov>) {
            return -1;
        } else if constexpr (
                (std::is_same_v<IndexTag1, Zeta>) || (std::is_same_v<IndexTag2, Phi_cov>)) {
            return 0;
        } else {
            return m_mapping_2d.template jacobian_component<IndexTag1, IndexTag2>(
                    CoordArg2D(coord));
        }
    }
};

namespace mapping_detail {
template <class Curvilinear2DToCartesian, class Zeta, class Phi>
struct MappingAccessibility<ExecSpace, ToroidalToCylindrical<Curvilinear2DToCartesian, Zeta, Phi>>
    : MappingAccessibility<ExecSpace, Curvilinear2DToCartesian>
{
};

} // namespace mapping_detail
