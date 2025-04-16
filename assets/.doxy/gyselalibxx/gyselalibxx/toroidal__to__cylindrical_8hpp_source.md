

# File toroidal\_to\_cylindrical.hpp

[**File List**](files.md) **>** [**mapping**](dir_5300298560c4bf255ab9f36681603d89.md) **>** [**toroidal\_to\_cylindrical.hpp**](toroidal__to__cylindrical_8hpp.md)

[Go to the documentation of this file](toroidal__to__cylindrical_8hpp.md)


```C++
// SPDX-License-Identifier: MIT
#pragma once

template <class Curvilinear2DToCartesian, class Zeta, class Phi>
class ToroidalToCylindrical
{
    static_assert(Zeta::IS_CONTRAVARIANT);
    static_assert(Phi::IS_CONTRAVARIANT);

public:
    using cylindrical_tag_R = typename Curvilinear2DToCartesian::cartesian_tag_x;
    using cylindrical_tag_Z = typename Curvilinear2DToCartesian::cartesian_tag_y;
    using cylindrical_tag_Zeta = Zeta;

    using toroidal_tag_rho = typename Curvilinear2DToCartesian::curvilinear_tag_r;
    using toroidal_tag_theta = typename Curvilinear2DToCartesian::curvilinear_tag_theta;
    using toroidal_tag_phi = Phi;

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

public:
    using CoordArg = Coord<Rho, Theta, Phi>;
    using CoordResult = Coord<R, Z, Zeta>;

private:
    Curvilinear2DToCartesian m_mapping_2d;

public:
    explicit ToroidalToCylindrical(Curvilinear2DToCartesian const& mapping_2d)
        : m_mapping_2d(mapping_2d)
    {
    }

    KOKKOS_FUNCTION CoordResult operator()(CoordArg const& coord) const
    {
        return CoordResult(m_mapping_2d(CoordArg2D(coord)), Coord<Zeta>(-ddc::get<Phi>(coord)));
    }

    KOKKOS_FUNCTION double jacobian(CoordArg const& coord) const
    {
        return -m_mapping_2d.jacobian(CoordArg2D(coord));
    }

    KOKKOS_FUNCTION DTensor<VectorIndexSet<R, Z, Zeta>, VectorIndexSet<Rho_cov, Theta_cov, Phi_cov>>
    jacobian_matrix(CoordArg const& coord) const
    {
        DTensor<VectorIndexSet<R, Z, Zeta>, VectorIndexSet<Rho_cov, Theta_cov, Phi_cov>> J_3d;
        DTensor<VectorIndexSet<R, Z>, VectorIndexSet<Rho_cov, Theta_cov>> J_2d
                = m_mapping_2d.jacobian_matrix(CoordArg2D(coord));
        ddcHelper::get<R, Rho_cov>(J_3d) = ddcHelper::get<R, Rho_cov>(J_2d);
        ddcHelper::get<R, Theta_cov>(J_3d) = ddcHelper::get<R, Theta_cov>(J_2d);
        ddcHelper::get<R, Phi_cov>(J_3d) = 0;
        ddcHelper::get<Z, Rho_cov>(J_3d) = ddcHelper::get<Z, Rho_cov>(J_2d);
        ddcHelper::get<Z, Theta_cov>(J_3d) = ddcHelper::get<Z, Theta_cov>(J_2d);
        ddcHelper::get<Z, Phi_cov>(J_3d) = 0;
        ddcHelper::get<Zeta, Rho_cov>(J_3d) = 0;
        ddcHelper::get<Zeta, Theta_cov>(J_3d) = 0;
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

    KOKKOS_FUNCTION DTensor<VectorIndexSet<Rho, Theta, Phi>, VectorIndexSet<R_cov, Z_cov, Zeta_cov>>
    inv_jacobian_matrix(CoordArg const& coord) const
    {
        InverseJacobianMatrix<Curvilinear2DToCartesian, CoordArg2D> inv_jacobian_matrix_2d(
                m_mapping_2d);
        DTensor<VectorIndexSet<Rho, Theta, Phi>, VectorIndexSet<R_cov, Z_cov, Zeta_cov>> inv_J_3d;
        DTensor<VectorIndexSet<Rho, Theta>, VectorIndexSet<R_cov, Z_cov>> inv_J_2d
                = inv_jacobian_matrix_2d(CoordArg2D(coord));
        ddcHelper::get<Rho, R_cov>(inv_J_3d) = ddcHelper::get<Rho, R_cov>(inv_J_2d);
        ddcHelper::get<Rho, Z_cov>(inv_J_3d) = ddcHelper::get<Rho, Z_cov>(inv_J_2d);
        ddcHelper::get<Rho, Zeta_cov>(inv_J_3d) = 0;
        ddcHelper::get<Theta, R_cov>(inv_J_3d) = ddcHelper::get<Theta, R_cov>(inv_J_2d);
        ddcHelper::get<Theta, Z_cov>(inv_J_3d) = ddcHelper::get<Theta, Z_cov>(inv_J_2d);
        ddcHelper::get<Theta, Zeta_cov>(inv_J_3d) = 0;
        ddcHelper::get<Phi, R_cov>(inv_J_3d) = 0;
        ddcHelper::get<Phi, Z_cov>(inv_J_3d) = 0;
        ddcHelper::get<Phi, Zeta_cov>(inv_J_3d) = -1;
        return inv_J_3d;
    }

    template <class IndexTag1, class IndexTag2>
    KOKKOS_INLINE_FUNCTION double inv_jacobian_component(CoordArg const& coord) const
    {
        static_assert(has_inv_jacobian_v<Curvilinear2DToCartesian, CoordArg>);
        static_assert(ddc::in_tags_v<IndexTag1, ddc::detail::TypeSeq<Rho, Theta, Phi>>);
        static_assert(ddc::in_tags_v<IndexTag2, ddc::detail::TypeSeq<R_cov, Z_cov, Zeta_cov>>);

        if constexpr (std::is_same_v<IndexTag1, Phi> && std::is_same_v<IndexTag2, Zeta_cov>) {
            return -1;
        } else if constexpr (
                (std::is_same_v<IndexTag1, Phi>) || (std::is_same_v<IndexTag2, Zeta_cov>)) {
            return 0;
        } else {
            return m_mapping_2d.template inv_jacobian_component<IndexTag1, IndexTag2>(
                    CoordArg2D(coord));
        }
    }

    Curvilinear2DToCartesian get_2d_polar_mapping() const
    {
        return m_mapping_2d;
    }
};

namespace mapping_detail {
template <class Curvilinear2DToCartesian, class Zeta, class Phi, class ExecSpace>
struct MappingAccessibility<ExecSpace, ToroidalToCylindrical<Curvilinear2DToCartesian, Zeta, Phi>>
    : MappingAccessibility<ExecSpace, Curvilinear2DToCartesian>
{
};

} // namespace mapping_detail
```


