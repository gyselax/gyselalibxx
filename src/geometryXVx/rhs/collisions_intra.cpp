#include <iomanip>

#include <fluid_moments.hpp>
#include <pdi.h>

#include "collisions_intra.hpp"
#include "collisions_utils.hpp"


CollisionsIntra::CollisionsIntra(IDomainSpXVx const& mesh, double nustar0)
    : m_nustar0(nustar0)
    , m_fthresh(1.e-30)
    , m_nustar_profile_alloc(ddc::select<IDimSp, IDimX>(mesh))
    , m_gridvx_ghosted(
              ddc::DiscreteElement<ghosted_vx_point_sampling>(0),
              ddc::DiscreteVector<ghosted_vx_point_sampling>(ddc::select<IDimVx>(mesh).size() + 2))
    , m_gridvx_ghosted_staggered(
              ddc::DiscreteElement<ghosted_vx_staggered_point_sampling>(0),
              ddc::DiscreteVector<ghosted_vx_staggered_point_sampling>(
                      ddc::select<IDimVx>(mesh).size() + 1))
    , m_mesh_ghosted(ddc::select<IDimSp>(mesh), ddc::select<IDimX>(mesh), m_gridvx_ghosted)
    , m_mesh_ghosted_staggered(
              ddc::select<IDimSp>(mesh),
              ddc::select<IDimX>(mesh),
              m_gridvx_ghosted_staggered)

{
    // validity checks
    if (ddc::select<IDimSp>(mesh).size() != 2) {
        throw std::runtime_error("Intra species collisions requires two kinetic species.");
    }
    if (m_nustar0 == 0.) {
        throw std::invalid_argument("Collision operator should not be used with nustar0=0.");
    }

    double const vx0 = ddc::coordinate(ddc::select<IDimVx>(mesh).front());
    double const vx1 = ddc::coordinate(ddc::select<IDimVx>(mesh).front() + 1);
    double const vxN = ddc::coordinate(ddc::select<IDimVx>(mesh).back());
    double const vxNm1 = ddc::coordinate(ddc::select<IDimVx>(mesh).back() - 1);
    int const ncells(ddc::select<IDimVx>(mesh).size() - 1);
    if constexpr (uniform_edge_v) {
        double const step(ddc::step<IDimVx>());
        ddc::init_discrete_space(
                ddc::UniformPointSampling<GhostedVx>::
                        init(ddc::Coordinate<GhostedVx>(vx0 - step),
                             ddc::Coordinate<GhostedVx>(vxN + step),
                             ddc::DiscreteVector<ddc::UniformPointSampling<GhostedVx>>(
                                     ncells + 3)));
    } else {
        int const npoints(ncells + 3);
        std::vector<ddc::Coordinate<GhostedVx>> breaks(npoints);
        breaks[0] = ddc::Coordinate<GhostedVx>(vx0 - (vx1 - vx0));
        breaks[npoints - 1] = ddc::Coordinate<GhostedVx>(vxN + (vxN - vxNm1));
        ddc::for_each(ddc::select<IDimVx>(mesh), [&](IndexVx const ivx) {
            breaks[ghosted_from_index(ivx).uid()] = ghosted_from_coord(ddc::coordinate(ivx));
        });
        ddc::init_discrete_space<ddc::NonUniformPointSampling<GhostedVx>>(breaks);
    }

    if constexpr (uniform_edge_v) {
        double const step(ddc::step<IDimVx>());
        ddc::init_discrete_space(
                ddc::UniformPointSampling<GhostedVxStaggered>::
                        init(ddc::Coordinate<GhostedVxStaggered>(vx0 - step / 2),
                             ddc::Coordinate<GhostedVxStaggered>(vxN + step / 2),
                             ddc::DiscreteVector<ddc::UniformPointSampling<GhostedVxStaggered>>(
                                     ncells + 2)));
    } else {
        int const npoints(ncells + 2);
        std::vector<ddc::Coordinate<GhostedVxStaggered>> breaks(npoints);
        breaks[0] = ddc::Coordinate<GhostedVxStaggered>(vx0 - (vx1 - vx0) / 2.);
        breaks[npoints - 1] = ddc::Coordinate<GhostedVxStaggered>(vxN + (vxN - vxNm1) / 2.);
        IDomainVx const gridvx_less(ddc::select<IDimVx>(mesh).remove_last(IVectVx(1)));
        ddc::for_each(gridvx_less, [&](IndexVx const ivx) {
            breaks[ivx.uid() + 1] = CollisionsIntra::ghosted_staggered_from_coord(
                    CoordVx((ddc::coordinate(ivx) + ddc::coordinate(ivx + 1)) / 2.));
        });
        ddc::init_discrete_space<ddc::NonUniformPointSampling<GhostedVxStaggered>>(breaks);
    }

    m_nustar_profile = m_nustar_profile_alloc.span_view();
    compute_nustar_profile(m_nustar_profile, m_nustar0);
    ddc::expose_to_pdi("collintra_nustar0", m_nustar0);
}

ddc::DiscreteDomain<CollisionsIntra::ghosted_vx_point_sampling> const& CollisionsIntra::
        get_gridvx_ghosted() const
{
    return m_gridvx_ghosted;
}

ddc::DiscreteDomain<CollisionsIntra::ghosted_vx_staggered_point_sampling> const& CollisionsIntra::
        get_gridvx_ghosted_staggered() const
{
    return m_gridvx_ghosted_staggered;
}

ddc::DiscreteDomain<IDimSp, IDimX, CollisionsIntra::ghosted_vx_point_sampling> const&
CollisionsIntra::get_mesh_ghosted() const
{
    return m_mesh_ghosted;
}

void CollisionsIntra::compute_matrix_coeff(
        device_t<DSpanSpXVx> AA,
        device_t<DSpanSpXVx> BB,
        device_t<DSpanSpXVx> CC,
        device_t<ddc::ChunkSpan<double, IDomainSpXVx_ghosted>> Dcoll,
        device_t<ddc::ChunkSpan<double, IDomainSpXVx_ghosted_staggered>> Dcoll_staggered,
        device_t<ddc::ChunkSpan<double, IDomainSpXVx_ghosted>> Nucoll,
        double deltat) const
{
    ddc::for_each(
            ddc::policies::parallel_device,
            AA.domain(),
            KOKKOS_LAMBDA(IndexSpXVx const ispxvx) {
                IndexSp const isp = ddc::select<IDimSp>(ispxvx);
                IndexX const ix = ddc::select<IDimX>(ispxvx);
                IndexVx const ivx = ddc::select<IDimVx>(ispxvx);

                IndexVx_ghosted ivx_ghosted(ghosted_from_index(ivx));
                IndexVx_ghosted_staggered ivx_ghosted_staggered(ghosted_staggered_from_index(ivx));
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
        DViewVx AA,
        DViewVx BB,
        DViewVx CC) const
{
    matrix.set_element(0, 0, BB(IndexVx(0)));
    matrix.set_element(0, 1, CC(IndexVx(0)));

    int const npoints(ddc::get_domain<IDimVx>(AA).size());
    matrix.set_element(npoints - 1, npoints - 1, BB(IndexVx(npoints - 1)));
    matrix.set_element(npoints - 1, npoints - 2, AA(IndexVx(npoints - 1)));

    IDomainVx const gridvx_inner(
            ddc::get_domain<IDimVx>(AA).remove_first(IVectVx(1)).remove_last(IVectVx(1)));
    ddc::for_each(gridvx_inner, [&](IndexVx const ivx) {
        matrix.set_element(ivx.uid(), ivx.uid() - 1, AA(ivx));
        matrix.set_element(ivx.uid(), ivx.uid(), BB(ivx));
        matrix.set_element(ivx.uid(), ivx.uid() + 1, CC(ivx));
    });
}

void CollisionsIntra::compute_rhs_vector(
        device_t<DSpanSpXVx> RR,
        device_t<DViewSpXVx> AA,
        device_t<DViewSpXVx> BB,
        device_t<DViewSpXVx> CC,
        device_t<DViewSpXVx> allfdistribu,
        double fthresh) const
{
    ddc::for_each(
            ddc::policies::parallel_device,
            RR.domain(),
            KOKKOS_LAMBDA(IndexSpXVx const ispxvx) {
                IndexSp const isp = ddc::select<IDimSp>(ispxvx);
                IndexX const ix = ddc::select<IDimX>(ispxvx);
                IndexVx const ivx = ddc::select<IDimVx>(ispxvx);

                IndexVx const ivx_next = ivx + 1;
                IndexVx const ivx_prev = ivx - 1;

                if (ivx == AA.domain<IDimVx>().front()) {
                    RR(isp, ix, ivx) = (2. - BB(isp, ix, ivx)) * allfdistribu(isp, ix, ivx)
                                       + (-CC(isp, ix, ivx)) * allfdistribu(isp, ix, ivx_next)
                                       - 2. * AA(isp, ix, ivx) * fthresh;

                } else if (ivx == AA.domain<IDimVx>().back()) {
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



device_t<DSpanSpXVx> CollisionsIntra::operator()(device_t<DSpanSpXVx> allfdistribu, double dt) const
{
    Kokkos::Profiling::pushRegion("CollisionsIntra");
    auto allfdistribu_alloc = ddc::create_mirror_view_and_copy(allfdistribu);
    ddc::ChunkSpan allfdistribu_host = allfdistribu_alloc.span_view();

    IDomainSpX grid_sp_x(allfdistribu.domain<IDimSp, IDimX>());
    // density and temperature
    device_t<DFieldSpX> density_alloc(grid_sp_x);
    device_t<DFieldSpX> fluid_velocity_alloc(grid_sp_x);
    device_t<DFieldSpX> temperature_alloc(grid_sp_x);
    auto density = density_alloc.span_view();
    auto fluid_velocity = fluid_velocity_alloc.span_view();
    auto temperature = temperature_alloc.span_view();

    DFieldVx const quadrature_coeffs_host(
            trapezoid_quadrature_coefficients(ddc::get_domain<IDimVx>(allfdistribu)));
    auto quadrature_coeffs_alloc = ddc::create_mirror_view_and_copy(
            Kokkos::DefaultExecutionSpace(),
            quadrature_coeffs_host.span_view());
    auto quadrature_coeffs = quadrature_coeffs_alloc.span_view();

    //Moments computation
    ddc::fill(density, 0.);
    ddc::for_each(
            ddc::policies::parallel_device,
            grid_sp_x,
            KOKKOS_LAMBDA(IndexSpX const ispx) {
                double particle_flux(0);
                double momentum_flux(0);
                for (IndexVx const ivx : allfdistribu.domain<IDimVx>()) {
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
    device_t<DFieldSpX> collfreq_alloc(grid_sp_x);
    auto collfreq = collfreq_alloc.span_view();
    device_t<DFieldSpX> nustar_profile(grid_sp_x);
    ddc::deepcopy(nustar_profile, m_nustar_profile);
    compute_collfreq(collfreq, nustar_profile, density, temperature);

    // diffusion coefficient
    device_t<ddc::Chunk<double, IDomainSpXVx_ghosted>> Dcoll_alloc(m_mesh_ghosted);
    auto Dcoll = Dcoll_alloc.span_view();
    compute_Dcoll<ghosted_vx_point_sampling>(Dcoll, collfreq, density, temperature);

    device_t<ddc::Chunk<double, IDomainSpXVx_ghosted>> dvDcoll_alloc(m_mesh_ghosted);
    auto dvDcoll = dvDcoll_alloc.span_view();
    compute_dvDcoll<ghosted_vx_point_sampling>(dvDcoll, collfreq, density, temperature);

    device_t<ddc::Chunk<double, IDomainSpXVx_ghosted_staggered>> Dcoll_staggered_alloc(
            m_mesh_ghosted_staggered);
    auto Dcoll_staggered = Dcoll_staggered_alloc.span_view();
    compute_Dcoll<
            ghosted_vx_staggered_point_sampling>(Dcoll_staggered, collfreq, density, temperature);

    // kernel maxwellian fluid moments
    device_t<DFieldSpX> Vcoll_alloc(grid_sp_x);
    device_t<DFieldSpX> Tcoll_alloc(grid_sp_x);
    auto Vcoll = Vcoll_alloc.span_view();
    auto Tcoll = Tcoll_alloc.span_view();
    compute_Vcoll_Tcoll<ghosted_vx_point_sampling>(Vcoll, Tcoll, allfdistribu, Dcoll, dvDcoll);

    // convection coefficient Nucoll
    device_t<ddc::Chunk<double, IDomainSpXVx_ghosted>> Nucoll_alloc(m_mesh_ghosted);
    auto Nucoll = Nucoll_alloc.span_view();
    compute_Nucoll<ghosted_vx_point_sampling>(Nucoll, Dcoll, Vcoll, Tcoll);

    // matrix coefficients
    device_t<DFieldSpXVx> AA_alloc(allfdistribu.domain());
    device_t<DFieldSpXVx> BB_alloc(allfdistribu.domain());
    device_t<DFieldSpXVx> CC_alloc(allfdistribu.domain());
    auto AA = AA_alloc.span_view();
    auto BB = BB_alloc.span_view();
    auto CC = CC_alloc.span_view();
    compute_matrix_coeff(AA, BB, CC, Dcoll, Dcoll_staggered, Nucoll, dt);


    // rhs vector coefficient
    device_t<DFieldSpXVx> RR_alloc(allfdistribu.domain());
    auto RR = RR_alloc.span_view();
    compute_rhs_vector(RR, AA, BB, CC, allfdistribu, m_fthresh);

    DFieldSpXVx AA_host(allfdistribu.domain());
    DFieldSpXVx BB_host(allfdistribu.domain());
    DFieldSpXVx CC_host(allfdistribu.domain());
    DFieldSpXVx RR_host(allfdistribu.domain());
    ddc::deepcopy(AA_host, AA);
    ddc::deepcopy(BB_host, BB);
    ddc::deepcopy(CC_host, CC);
    ddc::deepcopy(RR_host, RR);

    ddc::for_each(grid_sp_x, [&](IndexSpX const ispx) {
        Matrix_Banded matrix(ddc::get_domain<IDimVx>(allfdistribu).size(), 1, 1);
        fill_matrix_with_coeff(
                matrix,
                AA_host[ispx].span_cview(),
                BB_host[ispx].span_cview(),
                CC_host[ispx].span_cview());

        DSpan1D RR_Span1D(
                RR_host[ispx].allocation_mdspan().data_handle(),
                ddc::get_domain<IDimVx>(allfdistribu).size());
        matrix.factorize();
        matrix.solve_inplace(RR_Span1D);
        ddc::deepcopy(allfdistribu_host[ispx], RR_host[ispx]);
    });

    ddc::deepcopy(allfdistribu, allfdistribu_host);
    Kokkos::Profiling::popRegion();
    return allfdistribu;
}
