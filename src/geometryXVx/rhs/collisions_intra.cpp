#include <iomanip>

#include <fluid_moments.hpp>
#include <pdi.h>

#include "sll/matrix_batch_tridiag.hpp"

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
        ddc::init_discrete_space<ghosted_vx_point_sampling>(
                ghosted_vx_point_sampling::
                        init(ddc::Coordinate<GhostedVx>(vx0 - step),
                             ddc::Coordinate<GhostedVx>(vxN + step),
                             ddc::DiscreteVector<ghosted_vx_point_sampling>(ncells + 3)));
    } else {
        int const npoints(ncells + 3);
        std::vector<ddc::Coordinate<GhostedVx>> breaks(npoints);
        breaks[0] = ddc::Coordinate<GhostedVx>(vx0 - (vx1 - vx0));
        breaks[npoints - 1] = ddc::Coordinate<GhostedVx>(vxN + (vxN - vxNm1));
        ddc::for_each(ddc::select<IDimVx>(mesh), [&](IndexVx const ivx) {
            breaks[ghosted_from_index(ivx).uid()] = ghosted_from_coord(ddc::coordinate(ivx));
        });
        ddc::init_discrete_space<ghosted_vx_point_sampling>(breaks);
    }

    if constexpr (uniform_edge_v) {
        double const step(ddc::step<IDimVx>());
        ddc::init_discrete_space<ghosted_vx_staggered_point_sampling>(
                ghosted_vx_staggered_point_sampling::
                        init(ddc::Coordinate<GhostedVxStaggered>(vx0 - step / 2),
                             ddc::Coordinate<GhostedVxStaggered>(vxN + step / 2),
                             ddc::DiscreteVector<ghosted_vx_staggered_point_sampling>(ncells + 2)));
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
        ddc::init_discrete_space<ghosted_vx_staggered_point_sampling>(breaks);
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
        DSpanSpXVx AA,
        DSpanSpXVx BB,
        DSpanSpXVx CC,
        device_t<ddc::ChunkSpan<double, IDomainSpXVx_ghosted>> Dcoll,
        device_t<ddc::ChunkSpan<double, IDomainSpXVx_ghosted_staggered>> Dcoll_staggered,
        device_t<ddc::ChunkSpan<double, IDomainSpXVx_ghosted>> Nucoll,
        double deltat) const
{
    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
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
        host_t<DViewVx> AA,
        host_t<DViewVx> BB,
        host_t<DViewVx> CC) const
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
        DSpanSpXVx RR,
        DViewSpXVx AA,
        DViewSpXVx BB,
        DViewSpXVx CC,
        DViewSpXVx allfdistribu,
        double fthresh) const
{
    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
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



DSpanSpXVx CollisionsIntra::operator()(DSpanSpXVx allfdistribu, double dt) const
{
    Kokkos::Profiling::pushRegion("CollisionsIntra");
    auto allfdistribu_alloc = ddc::create_mirror_view_and_copy(allfdistribu);
    ddc::ChunkSpan allfdistribu_host = allfdistribu_alloc.span_view();

    IDomainSpX grid_sp_x(allfdistribu.domain<IDimSp, IDimX>());
    // density and temperature
    DFieldSpX density_alloc(grid_sp_x);
    DFieldSpX fluid_velocity_alloc(grid_sp_x);
    DFieldSpX temperature_alloc(grid_sp_x);
    auto density = density_alloc.span_view();
    auto fluid_velocity = fluid_velocity_alloc.span_view();
    auto temperature = temperature_alloc.span_view();

    host_t<DFieldVx> const quadrature_coeffs_host(
            trapezoid_quadrature_coefficients(ddc::get_domain<IDimVx>(allfdistribu)));
    auto quadrature_coeffs_alloc = ddc::create_mirror_view_and_copy(
            Kokkos::DefaultExecutionSpace(),
            quadrature_coeffs_host.span_view());
    auto quadrature_coeffs = quadrature_coeffs_alloc.span_view();

    //Moments computation
    ddc::parallel_fill(density, 0.);
    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
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
    DFieldSpX collfreq_alloc(grid_sp_x);
    auto collfreq = collfreq_alloc.span_view();
    DFieldSpX nustar_profile(grid_sp_x);
    ddc::parallel_deepcopy(nustar_profile, m_nustar_profile);
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
    DFieldSpX Vcoll_alloc(grid_sp_x);
    DFieldSpX Tcoll_alloc(grid_sp_x);
    auto Vcoll = Vcoll_alloc.span_view();
    auto Tcoll = Tcoll_alloc.span_view();
    compute_Vcoll_Tcoll<ghosted_vx_point_sampling>(Vcoll, Tcoll, allfdistribu, Dcoll, dvDcoll);

    // convection coefficient Nucoll
    device_t<ddc::Chunk<double, IDomainSpXVx_ghosted>> Nucoll_alloc(m_mesh_ghosted);
    auto Nucoll = Nucoll_alloc.span_view();
    compute_Nucoll<ghosted_vx_point_sampling>(Nucoll, Dcoll, Vcoll, Tcoll);

    // matrix coefficients
    DFieldSpXVx AA_alloc(allfdistribu.domain());
    DFieldSpXVx BB_alloc(allfdistribu.domain());
    DFieldSpXVx CC_alloc(allfdistribu.domain());
    auto AA = AA_alloc.span_view();
    auto BB = BB_alloc.span_view();
    auto CC = CC_alloc.span_view();
    compute_matrix_coeff(AA, BB, CC, Dcoll, Dcoll_staggered, Nucoll, dt);


    // rhs vector coefficient
    DFieldSpXVx RR_alloc(allfdistribu.domain());
    auto RR = RR_alloc.span_view();
    compute_rhs_vector(RR, AA, BB, CC, allfdistribu, m_fthresh);

    int const batch_size = ddc::get_domain<IDimSp, IDimX>(allfdistribu).size();
    int const mat_size = ddc::get_domain<IDimVx>(allfdistribu).size();
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
    matrix.factorize();
    matrix.solve_inplace(RR_view);

    ddc::parallel_deepcopy(allfdistribu, RR);
    Kokkos::Profiling::popRegion();
    return allfdistribu;
}
