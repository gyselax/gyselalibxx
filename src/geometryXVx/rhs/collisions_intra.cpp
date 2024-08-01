#include <iomanip>

#include <fluid_moments.hpp>
#include <pdi.h>

#include "sll/matrix_batch_tridiag.hpp"

#include "collisions_intra.hpp"
#include "collisions_utils.hpp"

template <class TargetDim>
KOKKOS_FUNCTION Idx<TargetDim> CollisionsIntra::to_index(Idx<GridVx> const& index)
{
    static_assert(
            std::is_same_v<TargetDim, GhostedVx> || std::is_same_v<TargetDim, GhostedVxStaggered>);
    if constexpr (std::is_same_v<TargetDim, GhostedVx>) {
        return Idx<GhostedVx>(index.uid() + 1);
    } else {
        return Idx<GhostedVxStaggered>(index.uid() + 1);
    }
}

template <class VDim>
std::enable_if_t<!ddc::is_uniform_point_sampling_v<VDim>> CollisionsIntra::
        build_ghosted_staggered_vx_point_sampling(IdxRange<VDim> const& dom)
{
    static_assert(
            std::is_same_v<VDim, GridVx>,
            "The function is only designed to work with the GridVx dimension");

    CoordVx const v0 = ddc::coordinate(dom.front());
    CoordVx const v1 = ddc::coordinate(dom.front() + 1);
    CoordVx const vN = ddc::coordinate(dom.back());
    CoordVx const vNm1 = ddc::coordinate(dom.back() - 1);
    int const ncells(dom.size() - 1);

    // ghosted points
    int const npoints(ncells + 3);
    std::vector<CoordVx> breaks(npoints);
    breaks[0] = v0 - (v1 - v0);
    breaks[npoints - 1] = vN + (vN - vNm1);
    ddc::for_each(dom, [&](IdxVx const iv) {
        breaks[to_index<GhostedVx>(iv).uid()] = ddc::coordinate(iv);
    });
    ddc::init_discrete_space<GhostedVx>(breaks);

    // ghosted staggered points
    int const npoints_stag(ncells + 2);
    std::vector<CoordVx> breaks_stag(npoints_stag);
    breaks_stag[0] = v0 - (v1 - v0) / 2.;
    breaks_stag[npoints_stag - 1] = vN + (vN - vNm1) / 2.;
    IdxRangeVx const gridv_less(dom.remove_last(IdxStepVx(1)));
    ddc::for_each(gridv_less, [&](IdxVx const iv) {
        breaks_stag[iv.uid() + 1] = CoordVx((ddc::coordinate(iv) + ddc::coordinate(iv + 1)) / 2.);
    });
    ddc::init_discrete_space<GhostedVxStaggered>(breaks_stag);
}

template <class VDim>
std::enable_if_t<ddc::is_uniform_point_sampling_v<VDim>> CollisionsIntra::
        build_ghosted_staggered_vx_point_sampling(IdxRange<VDim> const& dom)
{
    static_assert(
            std::is_same_v<VDim, GridVx>,
            "The function is only designed to work with the GridVx dimension");

    CoordVx const v0 = ddc::coordinate(dom.front());
    CoordVx const vN = ddc::coordinate(dom.back());
    int const ncells(dom.size() - 1);
    double const step(ddc::step<VDim>());

    // ghosted points
    ddc::init_discrete_space<GhostedVx>(
            GhostedVx::init(v0 - step, vN + step, IdxStep<GhostedVx>(ncells + 3)));

    // ghosted staggered points
    ddc::init_discrete_space<GhostedVxStaggered>(
            GhostedVxStaggered::
                    init(v0 - step / 2, vN + step / 2, IdxStep<GhostedVxStaggered>(ncells + 2)));
}

CollisionsIntra::CollisionsIntra(IdxRangeSpXVx const& mesh, double nustar0)
    : m_nustar0(nustar0)
    , m_fthresh(1.e-30)
    , m_nustar_profile_alloc(ddc::select<Species, GridX>(mesh))
    , m_gridvx_ghosted(Idx<GhostedVx>(0), IdxStep<GhostedVx>(ddc::select<GridVx>(mesh).size() + 2))
    , m_gridvx_ghosted_staggered(
              Idx<GhostedVxStaggered>(0),
              IdxStep<GhostedVxStaggered>(ddc::select<GridVx>(mesh).size() + 1))
    , m_mesh_ghosted(ddc::select<Species>(mesh), ddc::select<GridX>(mesh), m_gridvx_ghosted)
    , m_mesh_ghosted_staggered(
              ddc::select<Species>(mesh),
              ddc::select<GridX>(mesh),
              m_gridvx_ghosted_staggered)

{
    // validity checks
    if (ddc::select<Species>(mesh).size() != 2) {
        throw std::runtime_error("Intra species collisions requires two kinetic species.");
    }
    if (m_nustar0 == 0.) {
        throw std::invalid_argument("Collision operator should not be used with nustar0=0.");
    }

    build_ghosted_staggered_vx_point_sampling(ddc::select<GridVx>(mesh));

    m_nustar_profile = get_field(m_nustar_profile_alloc);
    compute_nustar_profile(m_nustar_profile, m_nustar0);
    ddc::expose_to_pdi("collintra_nustar0", m_nustar0);
}

IdxRange<CollisionsIntra::GhostedVx> const& CollisionsIntra::get_gridvx_ghosted() const
{
    return m_gridvx_ghosted;
}

IdxRange<CollisionsIntra::GhostedVxStaggered> const& CollisionsIntra::get_gridvx_ghosted_staggered()
        const
{
    return m_gridvx_ghosted_staggered;
}

IdxRange<Species, GridX, CollisionsIntra::GhostedVx> const& CollisionsIntra::get_mesh_ghosted()
        const
{
    return m_mesh_ghosted;
}

void CollisionsIntra::compute_matrix_coeff(
        DFieldSpXVx AA,
        DFieldSpXVx BB,
        DFieldSpXVx CC,
        Field<double, IDomainSpXVx_ghosted> Dcoll,
        Field<double, IDomainSpXVx_ghosted_staggered> Dcoll_staggered,
        Field<double, IDomainSpXVx_ghosted> Nucoll,
        double deltat) const
{
    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            get_idx_range(AA),
            KOKKOS_LAMBDA(IdxSpXVx const ispxvx) {
                IdxSp const isp = ddc::select<Species>(ispxvx);
                IdxX const ix = ddc::select<GridX>(ispxvx);
                IdxVx const ivx = ddc::select<GridVx>(ispxvx);

                IndexVx_ghosted ivx_ghosted(to_index<GhostedVx>(ivx));
                IndexVx_ghosted_staggered ivx_ghosted_staggered(to_index<GhostedVxStaggered>(ivx));
                IndexVx_ghosted ivx_next_ghosted(ivx_ghosted + 1);
                IndexVx_ghosted ivx_prev_ghosted(ivx_ghosted - 1);
                IndexVx_ghosted_staggered ivx_prev_ghosted_staggered(ivx_ghosted_staggered - 1);

                double const dv_i
                        = ddc::coordinate(ivx_next_ghosted) - ddc::coordinate(ivx_ghosted);
                double const delta_i
                        = dv_i / (ddc::coordinate(ivx_ghosted) - ddc::coordinate(ivx_prev_ghosted));

                double const alpha_i = deltat / (dv_i * dv_i * (1. + delta_i));
                double const beta_i = deltat / (2. * dv_i * (1. + delta_i));

                double const coeffa
                        = alpha_i
                                  * (Dcoll_staggered(isp, ix, ivx_prev_ghosted_staggered) * delta_i
                                             * delta_i * delta_i
                                     - Dcoll(isp, ix, ivx_ghosted) * delta_i * delta_i
                                               * (delta_i - 1.))
                          + beta_i * Nucoll(isp, ix, ivx_prev_ghosted) * delta_i * delta_i;

                double const coeffb
                        = -alpha_i
                                  * (-Dcoll_staggered(isp, ix, ivx_ghosted_staggered)
                                     - Dcoll_staggered(isp, ix, ivx_prev_ghosted_staggered)
                                               * delta_i * delta_i * delta_i
                                     + Dcoll(isp, ix, ivx_ghosted) * (delta_i - 1.)
                                               * (delta_i * delta_i - 1.))
                          + beta_i * Nucoll(isp, ix, ivx_ghosted) * (delta_i * delta_i - 1.);

                double const coeffc = alpha_i
                                              * (Dcoll_staggered(isp, ix, ivx_ghosted_staggered)
                                                 + Dcoll(isp, ix, ivx_ghosted) * (delta_i - 1.))
                                      - beta_i * Nucoll(isp, ix, ivx_next_ghosted);

                AA(ispxvx) = -coeffa;
                BB(ispxvx) = 1. + coeffb;
                CC(ispxvx) = -coeffc;
            });
}

void CollisionsIntra::fill_matrix_with_coeff(
        Matrix_Banded& matrix,
        host_t<DConstFieldVx> AA,
        host_t<DConstFieldVx> BB,
        host_t<DConstFieldVx> CC) const
{
    matrix.set_element(0, 0, BB(IdxVx(0)));
    matrix.set_element(0, 1, CC(IdxVx(0)));

    int const npoints(get_idx_range<GridVx>(AA).size());
    matrix.set_element(npoints - 1, npoints - 1, BB(IdxVx(npoints - 1)));
    matrix.set_element(npoints - 1, npoints - 2, AA(IdxVx(npoints - 1)));

    IdxRangeVx const gridvx_inner(
            get_idx_range<GridVx>(AA).remove_first(IdxStepVx(1)).remove_last(IdxStepVx(1)));
    ddc::for_each(gridvx_inner, [&](IdxVx const ivx) {
        matrix.set_element(ivx.uid(), ivx.uid() - 1, AA(ivx));
        matrix.set_element(ivx.uid(), ivx.uid(), BB(ivx));
        matrix.set_element(ivx.uid(), ivx.uid() + 1, CC(ivx));
    });
}

void CollisionsIntra::compute_rhs_vector(
        DFieldSpXVx RR,
        DConstFieldSpXVx AA,
        DConstFieldSpXVx BB,
        DConstFieldSpXVx CC,
        DConstFieldSpXVx allfdistribu,
        double fthresh) const
{
    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            get_idx_range(RR),
            KOKKOS_LAMBDA(IdxSpXVx const ispxvx) {
                IdxSp const isp = ddc::select<Species>(ispxvx);
                IdxX const ix = ddc::select<GridX>(ispxvx);
                IdxVx const ivx = ddc::select<GridVx>(ispxvx);

                IdxVx const ivx_next = ivx + 1;
                IdxVx const ivx_prev = ivx - 1;

                if (ivx == get_idx_range<GridVx>(AA).front()) {
                    RR(isp, ix, ivx) = (2. - BB(isp, ix, ivx)) * allfdistribu(isp, ix, ivx)
                                       + (-CC(isp, ix, ivx)) * allfdistribu(isp, ix, ivx_next)
                                       - 2. * AA(isp, ix, ivx) * fthresh;

                } else if (ivx == get_idx_range<GridVx>(AA).back()) {
                    RR(isp, ix, ivx) = (2. - BB(isp, ix, ivx)) * allfdistribu(isp, ix, ivx)
                                       + (-AA(isp, ix, ivx)) * allfdistribu(isp, ix, ivx_prev)
                                       - 2. * CC(isp, ix, ivx) * fthresh;

                } else {
                    RR(isp, ix, ivx) = -AA(isp, ix, ivx) * allfdistribu(isp, ix, ivx_prev)
                                       + (2. - BB(isp, ix, ivx)) * allfdistribu(isp, ix, ivx)
                                       - CC(isp, ix, ivx) * allfdistribu(isp, ix, ivx_next);
                }
            });
}



DFieldSpXVx CollisionsIntra::operator()(DFieldSpXVx allfdistribu, double dt) const
{
    Kokkos::Profiling::pushRegion("CollisionsIntra");
    auto allfdistribu_alloc = ddc::create_mirror_view_and_copy(allfdistribu);
    ddc::ChunkSpan allfdistribu_host = get_field(allfdistribu_alloc);

    IdxRangeSpX grid_sp_x(get_idx_range<Species, GridX>(allfdistribu));
    // density and temperature
    DFieldMemSpX density_alloc(grid_sp_x);
    DFieldMemSpX fluid_velocity_alloc(grid_sp_x);
    DFieldMemSpX temperature_alloc(grid_sp_x);
    auto density = get_field(density_alloc);
    auto fluid_velocity = get_field(fluid_velocity_alloc);
    auto temperature = get_field(temperature_alloc);

    DFieldMemVx const quadrature_coeffs_alloc(
            trapezoid_quadrature_coefficients<Kokkos::DefaultExecutionSpace>(
                    ddc::get_domain<GridVx>(allfdistribu)));
    auto quadrature_coeffs = get_field(quadrature_coeffs_alloc);

    //Moments computation
    ddc::parallel_fill(density, 0.);
    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            grid_sp_x,
            KOKKOS_LAMBDA(IdxSpX const ispx) {
                double particle_flux(0);
                double momentum_flux(0);
                for (IdxVx const ivx : get_idx_range<GridVx>(allfdistribu)) {
                    CoordVx const coordv = ddc::coordinate(ivx);
                    double const val(quadrature_coeffs(ivx) * allfdistribu(ispx, ivx));
                    density(ispx) += val;
                    particle_flux += val * coordv;
                    momentum_flux += val * coordv * coordv;
                }
                fluid_velocity(ispx) = particle_flux / density(ispx);
                temperature(ispx)
                        = (momentum_flux - particle_flux * fluid_velocity(ispx)) / density(ispx);
            });

    // collision frequency
    DFieldMemSpX collfreq_alloc(grid_sp_x);
    auto collfreq = get_field(collfreq_alloc);
    DFieldMemSpX nustar_profile(grid_sp_x);
    ddc::parallel_deepcopy(nustar_profile, m_nustar_profile);
    compute_collfreq(collfreq, nustar_profile, density, temperature);

    // diffusion coefficient
    FieldMem<double, IDomainSpXVx_ghosted> Dcoll_alloc(m_mesh_ghosted);
    auto Dcoll = get_field(Dcoll_alloc);
    compute_Dcoll<GhostedVx>(Dcoll, collfreq, density, temperature);

    FieldMem<double, IDomainSpXVx_ghosted> dvDcoll_alloc(m_mesh_ghosted);
    auto dvDcoll = get_field(dvDcoll_alloc);
    compute_dvDcoll<GhostedVx>(dvDcoll, collfreq, density, temperature);

    FieldMem<double, IDomainSpXVx_ghosted_staggered> Dcoll_staggered_alloc(
            m_mesh_ghosted_staggered);
    auto Dcoll_staggered = get_field(Dcoll_staggered_alloc);
    compute_Dcoll<GhostedVxStaggered>(Dcoll_staggered, collfreq, density, temperature);

    // kernel maxwellian fluid moments
    DFieldMemSpX Vcoll_alloc(grid_sp_x);
    DFieldMemSpX Tcoll_alloc(grid_sp_x);
    auto Vcoll = get_field(Vcoll_alloc);
    auto Tcoll = get_field(Tcoll_alloc);
    compute_Vcoll_Tcoll<GhostedVx>(Vcoll, Tcoll, allfdistribu, Dcoll, dvDcoll);

    // convection coefficient Nucoll
    FieldMem<double, IDomainSpXVx_ghosted> Nucoll_alloc(m_mesh_ghosted);
    auto Nucoll = get_field(Nucoll_alloc);
    compute_Nucoll<GhostedVx>(Nucoll, Dcoll, Vcoll, Tcoll);

    // matrix coefficients
    DFieldMemSpXVx AA_alloc(get_idx_range(allfdistribu));
    DFieldMemSpXVx BB_alloc(get_idx_range(allfdistribu));
    DFieldMemSpXVx CC_alloc(get_idx_range(allfdistribu));
    auto AA = get_field(AA_alloc);
    auto BB = get_field(BB_alloc);
    auto CC = get_field(CC_alloc);
    compute_matrix_coeff(AA, BB, CC, Dcoll, Dcoll_staggered, Nucoll, dt);


    // rhs vector coefficient
    DFieldMemSpXVx RR_alloc(get_idx_range(allfdistribu));
    auto RR = get_field(RR_alloc);
    compute_rhs_vector(RR, AA, BB, CC, allfdistribu, m_fthresh);

    int const batch_size = get_idx_range<Species, GridX>(allfdistribu).size();
    int const mat_size = get_idx_range<GridVx>(allfdistribu).size();
    /* Here we do not use allocation_kokkos_view() ddc function since we change the shape 
       from (Sp,X,Vx)-->(batch_dim,Vx)*/
    Kokkos::View<double**, Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace>
            AA_view(AA.data_handle(), batch_size, mat_size);
    Kokkos::View<double**, Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace>
            BB_view(BB.data_handle(), batch_size, mat_size);
    Kokkos::View<double**, Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace>
            CC_view(CC.data_handle(), batch_size, mat_size);
    Kokkos::View<double**, Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace>
            RR_view(RR.data_handle(), batch_size, mat_size);


    MatrixBatchTridiag<Kokkos::DefaultExecutionSpace>
            matrix(batch_size, mat_size, AA_view, BB_view, CC_view);
    matrix.setup_solver();
    matrix.solve(RR_view);

    ddc::parallel_deepcopy(allfdistribu, RR);
    Kokkos::Profiling::popRegion();
    return allfdistribu;
}
