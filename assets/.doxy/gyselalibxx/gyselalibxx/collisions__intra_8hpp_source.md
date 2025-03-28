

# File collisions\_intra.hpp

[**File List**](files.md) **>** [**geometryXVx**](dir_e51b496b46dd687775e46e0826614574.md) **>** [**rhs**](dir_53474cb30a3389ee74cb3186cae99ac0.md) **>** [**collisions\_intra.hpp**](collisions__intra_8hpp.md)

[Go to the documentation of this file](collisions__intra_8hpp.md)


```C++
// SPDX-License-Identifier: MIT
#pragma once
#include <cassert>
#include <cmath>

#include <ddc/ddc.hpp>

#include "ddc_aliases.hpp"
#include "geometry.hpp"
#include "irighthandside.hpp"
#include "matrix_banded.hpp"
#include "quadrature.hpp"
#include "trapezoid_quadrature.hpp"

class CollisionsIntra : public IRightHandSide
{
private:
    static constexpr bool uniform_edge_v = ddc::is_uniform_point_sampling_v<GridVx>;

public:
    struct GhostedVx
        : std::conditional_t<uniform_edge_v, UniformGridBase<Vx>, NonUniformGridBase<Vx>>
    {
    };

    struct GhostedVxStaggered
        : std::conditional_t<uniform_edge_v, UniformGridBase<Vx>, NonUniformGridBase<Vx>>
    {
    };

    using IdxRangeSpXVx_ghosted = IdxRange<Species, GridX, GhostedVx>;

    using IdxRangeSpXVx_ghosted_staggered = IdxRange<Species, GridX, GhostedVxStaggered>;

    using IdxVx_ghosted = Idx<GhostedVx>;

    using IdxVx_ghosted_staggered = Idx<GhostedVxStaggered>;

    using IdxSpXVx_ghosted = Idx<Species, GridX, GhostedVx>;

    using IdxSpXVx_ghosted_staggered = Idx<Species, GridX, GhostedVxStaggered>;


private:
    template <class TargetDim>
    KOKKOS_FUNCTION static Idx<TargetDim> to_index(Idx<GridVx> const& index);

    // The "Spoof" variables will be identical to the non-spoof versions. They are simply used
    // to prevent the compiler from trying to compile code for the non-uniform case when splines
    // are uniform.
    template <
            class GridVxSpoof,
            class GhostedVxSpoof = GhostedVx,
            class GhostedVxStaggeredSpoof = GhostedVxStaggered>
    std::enable_if_t<!ddc::is_uniform_point_sampling_v<GridVxSpoof>>
    build_ghosted_staggered_vx_point_sampling(IdxRange<GridVxSpoof> const& idx_range);

    // The "Spoof" variables will be identical to the non-spoof versions. They are simply used
    // to prevent the compiler from trying to compile code for the uniform case when splines
    // are non-uniform.
    template <
            class GridVxSpoof,
            class GhostedVxSpoof = GhostedVx,
            class GhostedVxStaggeredSpoof = GhostedVxStaggered>
    std::enable_if_t<ddc::is_uniform_point_sampling_v<GridVxSpoof>>
    build_ghosted_staggered_vx_point_sampling(IdxRange<GridVxSpoof> const& idx_range);

    double m_nustar0;
    double m_fthresh;
    DFieldMemSpX m_nustar_profile_alloc;
    DFieldSpX m_nustar_profile;

    IdxRange<GhostedVx> m_gridvx_ghosted;
    IdxRange<GhostedVxStaggered> m_gridvx_ghosted_staggered;

    IdxRangeSpXVx_ghosted m_mesh_ghosted;
    IdxRangeSpXVx_ghosted_staggered m_mesh_ghosted_staggered;

public:
    CollisionsIntra(IdxRangeSpXVx const& mesh, double nustar0);

    DFieldSpXVx operator()(DFieldSpXVx allfdistribu, double dt) const override;

    double get_nustar0() const;

    IdxRange<GhostedVx> const& get_gridvx_ghosted() const;

    IdxRange<GhostedVxStaggered> const& get_gridvx_ghosted_staggered() const;

    IdxRange<Species, GridX, GhostedVx> const& get_mesh_ghosted() const;

    void compute_rhs_vector(
            DFieldSpXVx RR,
            DConstFieldSpXVx AA,
            DConstFieldSpXVx BB,
            DConstFieldSpXVx CC,
            DConstFieldSpXVx allfdistribu,
            double fthresh) const;
    void compute_matrix_coeff(
            DFieldSpXVx AA,
            DFieldSpXVx BB,
            DFieldSpXVx CC,
            Field<double, IdxRangeSpXVx_ghosted> Dcoll,
            Field<double, IdxRangeSpXVx_ghosted_staggered> Dcoll_staggered,
            Field<double, IdxRangeSpXVx_ghosted> Nucoll,
            double deltat) const;

    void fill_matrix_with_coeff(
            Matrix_Banded& matrix,
            host_t<DConstFieldVx> AA,
            host_t<DConstFieldVx> BB,
            host_t<DConstFieldVx> CC) const;
};
```


