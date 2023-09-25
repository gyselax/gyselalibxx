#pragma once

#include <cassert>
#include <cmath>

#include <ddc/ddc.hpp>

#include <sll/matrix_banded.hpp>

#include <geometry.hpp>
#include <irighthandside.hpp>
#include <quadrature.hpp>
#include <trapezoid_quadrature.hpp>


class CollisionsIntra : public IRightHandSide
{
public:
    struct GhostedVx
    {
    };

    struct GhostedVxStaggered
    {
    };

private:
    static ddc::Coordinate<GhostedVx> ghosted_from_coord(ddc::Coordinate<RDimVx> const& coord)
    {
        return ddc::Coordinate<GhostedVx>(ddc::get<RDimVx>(coord));
    }
    static ddc::Coordinate<RDimVx> coord_from_ghosted(ddc::Coordinate<GhostedVx> const& coord)
    {
        return ddc::Coordinate<RDimVx>(ddc::get<GhostedVx>(coord));
    }
    static ddc::Coordinate<GhostedVxStaggered> ghosted_staggered_from_coord(
            ddc::Coordinate<RDimVx> const& coord)
    {
        return ddc::Coordinate<GhostedVxStaggered>(ddc::get<RDimVx>(coord));
    }
    static ddc::Coordinate<RDimVx> coord_from_ghosted_staggered(
            ddc::Coordinate<GhostedVxStaggered> const& coord)
    {
        return ddc::Coordinate<RDimVx>(ddc::get<GhostedVxStaggered>(coord));
    }

    static IndexVx index_from_ghosted(ddc::DiscreteElement<GhostedVx> const& index_ghosted)
    {
        return IndexVx(index_ghosted.uid() - 1);
    }
    static ddc::DiscreteElement<GhostedVx> ghosted_from_index(IndexVx const& index)
    {
        return ddc::DiscreteElement<GhostedVx>(index.uid() + 1);
    }

    static constexpr bool uniform_edge_v
            = std::is_same_v<IDimVx, ddc::UniformPointSampling<RDimVx>>;

public:
    using ghosted_vx_point_sampling = std::conditional_t<
            uniform_edge_v,
            ddc::UniformPointSampling<GhostedVx>,
            ddc::NonUniformPointSampling<GhostedVx>>;

    using ghosted_vx_staggered_point_sampling = std::conditional_t<
            uniform_edge_v,
            ddc::UniformPointSampling<GhostedVxStaggered>,
            ddc::NonUniformPointSampling<GhostedVxStaggered>>;

    using IDomainSpXVx_ghosted = ddc::DiscreteDomain<IDimSp, IDimX, ghosted_vx_point_sampling>;

    using IDomainSpXVx_ghosted_staggered
            = ddc::DiscreteDomain<IDimSp, IDimX, ghosted_vx_staggered_point_sampling>;

    using IndexVx_ghosted = ddc::DiscreteElement<ghosted_vx_point_sampling>;

    using IndexVx_ghosted_staggered = ddc::DiscreteElement<ghosted_vx_staggered_point_sampling>;

    using IndexSpXVx_ghosted = ddc::DiscreteElement<IDimSp, IDimX, ghosted_vx_point_sampling>;

    using IndexSpXVx_ghosted_staggered
            = ddc::DiscreteElement<IDimSp, IDimX, ghosted_vx_staggered_point_sampling>;

private:
    double m_nustar0;
    double m_fthresh;
    DFieldSpX m_nustar_profile;

    ddc::DiscreteDomain<ghosted_vx_point_sampling> m_gridvx_ghosted;
    ddc::DiscreteDomain<ghosted_vx_staggered_point_sampling> m_gridvx_ghosted_staggered;

    IDomainSpXVx_ghosted m_mesh_ghosted;
    IDomainSpXVx_ghosted_staggered m_mesh_ghosted_staggered;

public:
    CollisionsIntra(IDomainSpXVx const& mesh, double nustar0);

    ~CollisionsIntra() = default;

    DSpanSpXVx operator()(DSpanSpXVx allfdistribu, double dt) const override;

    double get_nustar0() const;

    ddc::DiscreteDomain<ghosted_vx_point_sampling> const& get_gridvx_ghosted() const;

    ddc::DiscreteDomain<ghosted_vx_staggered_point_sampling> const& get_gridvx_ghosted_staggered()
            const;

    ddc::DiscreteDomain<IDimSp, IDimX, ghosted_vx_point_sampling> const& get_mesh_ghosted() const;

    void compute_rhs_vector(
            DSpanSpXVx RR,
            DViewSpXVx AA,
            DViewSpXVx BB,
            DViewSpXVx CC,
            DViewSpXVx allfdistribu,
            double fthresh) const;

    void compute_matrix_coeff(
            DSpanSpXVx AA,
            DSpanSpXVx BB,
            DSpanSpXVx CC,
            ddc::ChunkSpan<double, IDomainSpXVx_ghosted> Dcoll,
            ddc::ChunkSpan<double, IDomainSpXVx_ghosted_staggered> Dcoll_staggered,
            ddc::ChunkSpan<double, IDomainSpXVx_ghosted> Nucoll,
            double deltat) const;

    void fill_matrix_with_coeff(Matrix_Banded& matrix, DViewVx AA, DViewVx BB, DViewVx CC) const;
};
