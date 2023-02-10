#include <iomanip>

#include <fluid_moments.hpp>
#include <pdi.h>

#include "collisions_intra.hpp"
#include "collisions_utils.hpp"

/**
 *  Treatment of the collision operator
 *  We solve the following equation:
 *     df_a/dt = C_aa + C_ab
 *
 *  C_aa= d/dv (Dcoll df_a/dv - Nucoll f_a)
 *    is the single species collision operator
 *  - Dcoll and Nucoll depend on the module of the velocity
 *    and can have a spatial profile
 *  - Nucoll = Dcoll * d/dv(Ln(fMcoll))
 *  - fMcoll is a Maxwellian of velocity and temperature Vcoll and Tcoll
 *    so that d/dv(Ln(fMcoll)) = -(v-Vcoll)/Tcoll
 *  - Tcoll(x,t) and Vcoll(x,t) are computed at each time step
 *  - nustar (hence the collision frequency) depends on space
 *
 *  'Vcoll' and 'Tcoll' are defined as follows:
 *  - Tcoll = Pcoll^{-1}[Imean0*Imean2 - Imean1*Imean1]
 *  - Vcoll = Pcoll^{-1}[Imean4*Imean1 - Imean3*Imean2]
 *  - Pcoll = Imean0*Imean4 - Imean1*Imean3
 *  where the 5 integrals are defined as:
 *     Imean0=<Dcoll> ;
 *     Imean1=<v*Dcoll> ;
 *     Imean2=<v^2*Dcoll> ;
 *     Imean3=<d/dv(Dcoll)>
 *     Imean4=<d/dv(v*Dcoll>
 *  The brackets <.> represent the integral in velocity: <.> = \int . dv
 *
 *  C_ab = { 2 Q_ab [m_a(v-V_a)^2/2T_a - 1/2]
 *         + R_ab (v-V_a) } * FM_a / (n_a Ta)
 *    accounts for total energy Q_ab+V_a*R_ab & momentum R_ab exchange between species
 *    Q_ab and R_ab are computed from the equivalent Maxwellians
 *    => order 0 of the collision operator
 */
CollisionsIntra::CollisionsIntra(IDomainSpXVx const& mesh, double nustar0)
    : m_nustar0(nustar0)
    , m_fthresh(1.e-30)
    , m_gridvx_ghosted(
              DiscreteElement<ghosted_vx_point_sampling>(0),
              DiscreteVector<ghosted_vx_point_sampling>(select<IDimVx>(mesh).size() + 2))
    , m_gridvx_ghosted_staggered(
              DiscreteElement<ghosted_vx_staggered_point_sampling>(0),
              DiscreteVector<ghosted_vx_staggered_point_sampling>(select<IDimVx>(mesh).size() + 1))
    , m_mesh_ghosted(select<IDimSp>(mesh), select<IDimX>(mesh), m_gridvx_ghosted)
    , m_mesh_ghosted_staggered(
              select<IDimSp>(mesh),
              select<IDimX>(mesh),
              m_gridvx_ghosted_staggered)

{
    // validity checks
    if (select<IDimSp>(mesh).size() != 2) {
        throw std::runtime_error("Intra species collisions requires two kinetic species.");
    }
    if (m_nustar0 == 0.) {
        throw std::invalid_argument("Collision operator should not be used with nustar0=0.");
    }

    double const vx0 = coordinate(select<IDimVx>(mesh).front());
    double const vx1 = coordinate(select<IDimVx>(mesh).front() + 1);
    double const vxN = coordinate(select<IDimVx>(mesh).back());
    double const vxNm1 = coordinate(select<IDimVx>(mesh).back() - 1);
    int const ncells(select<IDimVx>(mesh).size() - 1);
    if constexpr (uniform_edge_v) {
        double const step(ddc::step<IDimVx>());
        init_discrete_space<UniformPointSampling<GhostedVx>>(vx0 - step, vxN + step, ncells + 3);
    } else {
        int const npoints(ncells + 3);
        std::vector<Coordinate<GhostedVx>> breaks(npoints);
        breaks[0] = Coordinate<GhostedVx>(vx0 - (vx1 - vx0));
        breaks[npoints - 1] = Coordinate<GhostedVx>(vxN + (vxN - vxNm1));
        for_each(policies::parallel_host, select<IDimVx>(mesh), [&](IndexVx const ivx) {
            breaks[ghosted_from_index(ivx).uid()] = ghosted_from_coord(coordinate(ivx));
        });
        init_discrete_space<NonUniformPointSampling<GhostedVx>>(breaks);
    }

    if constexpr (uniform_edge_v) {
        double const step(ddc::step<IDimVx>());
        init_discrete_space<
                UniformPointSampling<GhostedVxStaggered>>(vx0 - step, vxN + step, ncells + 2);
    } else {
        int const npoints(ncells + 2);
        std::vector<Coordinate<GhostedVxStaggered>> breaks(npoints);
        breaks[0] = Coordinate<GhostedVxStaggered>(vx0 - (vx1 - vx0) / 2.);
        breaks[npoints - 1] = Coordinate<GhostedVxStaggered>(vxN + (vxN - vxNm1) / 2.);
        IDomainVx const gridvx_less(select<IDimVx>(mesh).remove_last(IVectVx(1)));
        for_each(policies::parallel_host, gridvx_less, [&](IndexVx const ivx) {
            breaks[ivx.uid() + 1] = CollisionsIntra::ghosted_staggered_from_coord(
                    CoordVx((coordinate(ivx) + coordinate(ivx + 1)) / 2.));
        });
        init_discrete_space<NonUniformPointSampling<GhostedVxStaggered>>(breaks);
    }
}

void CollisionsIntra::expose_rhs_to_pdi() const
{
    expose_to_pdi("collintra_nustar0", m_nustar0);
}

DiscreteDomain<CollisionsIntra::ghosted_vx_point_sampling> const& CollisionsIntra::
        get_gridvx_ghosted() const
{
    return m_gridvx_ghosted;
}

DiscreteDomain<CollisionsIntra::ghosted_vx_staggered_point_sampling> const& CollisionsIntra::
        get_gridvx_ghosted_staggered() const
{
    return m_gridvx_ghosted_staggered;
}

DiscreteDomain<IDimSp, IDimX, CollisionsIntra::ghosted_vx_point_sampling> const& CollisionsIntra::
        get_mesh_ghosted() const
{
    return m_mesh_ghosted;
}

/**
 * Computes the matrix coefficients for the non equidistant Crank Nicolson scheme
 */
void CollisionsIntra::compute_matrix_coeff(
        DSpanSpXVx AA,
        DSpanSpXVx BB,
        DSpanSpXVx CC,
        ddc::ChunkSpan<double, IDomainSpXVx_ghosted> Dcoll,
        ddc::ChunkSpan<double, IDomainSpXVx_ghosted_staggered> Dcoll_staggered,
        ddc::ChunkSpan<double, IDomainSpXVx_ghosted> Nucoll,
        double deltat) const
{
    for_each(policies::parallel_host, AA.domain(), [&](IndexSpXVx const ispxvx) {
        IndexSp const isp = select<IDimSp>(ispxvx);
        IndexX const ix = select<IDimX>(ispxvx);

        IndexVx_ghosted ivx_ghosted(select<IDimVx>(ispxvx).uid() + 1);
        IndexVx_ghosted_staggered ivx_ghosted_staggered(select<IDimVx>(ispxvx).uid() + 1);
        IndexVx_ghosted ivx_next_ghosted(ivx_ghosted + 1);
        IndexVx_ghosted ivx_prev_ghosted(ivx_ghosted - 1);
        IndexVx_ghosted_staggered ivx_prev_ghosted_staggered(ivx_ghosted_staggered - 1);

        double const dv_i = coordinate(ivx_next_ghosted) - coordinate(ivx_ghosted);
        double const delta_i = dv_i / (coordinate(ivx_ghosted) - coordinate(ivx_prev_ghosted));

        double const alpha_i = deltat / (dv_i * dv_i * (1. + delta_i));
        double const beta_i = deltat / (2. * dv_i * (1. + delta_i));

        double const coeffa
                = alpha_i
                          * (Dcoll_staggered(isp, ix, ivx_prev_ghosted_staggered) * delta_i
                                     * delta_i * delta_i
                             - Dcoll(isp, ix, ivx_ghosted) * delta_i * delta_i * (delta_i - 1.))
                  + beta_i * Nucoll(isp, ix, ivx_prev_ghosted) * delta_i * delta_i;

        double const coeffb = -alpha_i
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

    int const npoints(get_domain<IDimVx>(AA).size());
    matrix.set_element(npoints - 1, npoints - 1, BB(IndexVx(npoints - 1)));
    matrix.set_element(npoints - 1, npoints - 2, AA(IndexVx(npoints - 1)));

    IDomainVx const gridvx_inner(
            get_domain<IDimVx>(AA).remove_first(IVectVx(1)).remove_last(IVectVx(1)));
    for_each(gridvx_inner, [&](IndexVx const ivx) {
        matrix.set_element(ivx.uid(), ivx.uid() - 1, AA(ivx));
        matrix.set_element(ivx.uid(), ivx.uid(), BB(ivx));
        matrix.set_element(ivx.uid(), ivx.uid() + 1, CC(ivx));
    });
}

/**
 * Computes the rhs vector coefficients for the non equidistant Crank Nicolson scheme
 */
void CollisionsIntra::compute_rhs_vector(
        DSpanSpXVx RR,
        DViewSpXVx AA,
        DViewSpXVx BB,
        DViewSpXVx CC,
        DViewSpXVx allfdistribu,
        double fthresh) const
{
    for_each(policies::parallel_host, RR.domain(), [&](IndexSpXVx const ispxvx) {
        IndexSp const isp = select<IDimSp>(ispxvx);
        IndexX const ix = select<IDimX>(ispxvx);
        IndexVx const ivx = select<IDimVx>(ispxvx);

        IndexVx const ivx_next = ivx + 1;
        IndexVx const ivx_prev = ivx - 1;

        if (ivx == get_domain<IDimVx>(AA).front()) {
            RR(isp, ix, ivx) = (2. - BB(isp, ix, ivx)) * allfdistribu(isp, ix, ivx)
                               + (-CC(isp, ix, ivx)) * allfdistribu(isp, ix, ivx_next)
                               - 2. * AA(isp, ix, ivx) * fthresh;

        } else if (ivx == get_domain<IDimVx>(AA).back()) {
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
    // collisionality profile
    DFieldSpX nustar_profile(get_domain<IDimSp, IDimX>(allfdistribu));
    compute_nustar_profile(nustar_profile.span_view(), m_nustar0);

    // density and temperature
    DFieldSpX density(get_domain<IDimSp, IDimX>(allfdistribu));
    DFieldSpX mean_velocity(get_domain<IDimSp, IDimX>(allfdistribu));
    DFieldSpX temperature(get_domain<IDimSp, IDimX>(allfdistribu));
    FluidMoments moments(Quadrature<IDimVx>(
            trapezoid_quadrature_coefficients(get_domain<IDimVx>(allfdistribu))));

    moments(density.span_view(), allfdistribu.span_cview(), FluidMoments::s_density);
    moments(mean_velocity.span_view(),
            allfdistribu.span_cview(),
            density.span_cview(),
            FluidMoments::s_velocity);
    moments(temperature.span_view(),
            allfdistribu.span_cview(),
            density.span_cview(),
            mean_velocity.span_cview(),
            FluidMoments::s_temperature);

    // collision frequency
    DFieldSpX collfreq(get_domain<IDimSp, IDimX>(allfdistribu));
    compute_collfreq(
            collfreq.span_view(),
            nustar_profile.span_cview(),
            density.span_cview(),
            temperature.span_cview());

    // diffusion coefficient
    ddc::Chunk<double, IDomainSpXVx_ghosted> Dcoll(m_mesh_ghosted);
    compute_Dcoll<ghosted_vx_point_sampling>(
            Dcoll.span_view(),
            collfreq.span_cview(),
            density.span_cview(),
            temperature.span_cview());

    ddc::Chunk<double, IDomainSpXVx_ghosted> dvDcoll(m_mesh_ghosted);
    compute_dvDcoll<ghosted_vx_point_sampling>(
            dvDcoll.span_view(),
            collfreq.span_cview(),
            density.span_cview(),
            temperature.span_cview());

    ddc::Chunk<double, IDomainSpXVx_ghosted_staggered> Dcoll_staggered(m_mesh_ghosted_staggered);
    compute_Dcoll<ghosted_vx_staggered_point_sampling>(
            Dcoll_staggered.span_view(),
            collfreq.span_cview(),
            density.span_cview(),
            temperature.span_cview());

    // kernel maxwellian fluid moments
    DFieldSpX Vcoll(get_domain<IDimSp, IDimX>(allfdistribu));
    DFieldSpX Tcoll(get_domain<IDimSp, IDimX>(allfdistribu));
    compute_Vcoll_Tcoll(
            Vcoll.span_view(),
            Tcoll.span_view(),
            allfdistribu.span_cview(),
            Dcoll.span_cview(),
            dvDcoll.span_cview());

    // convection coefficient Nucoll
    ddc::Chunk<double, IDomainSpXVx_ghosted> Nucoll(m_mesh_ghosted);
    compute_Nucoll<ghosted_vx_point_sampling>(
            Nucoll.span_view(),
            Dcoll.span_view(),
            Vcoll.span_cview(),
            Tcoll.span_cview());

    // matrix coefficients
    DFieldSpXVx AA(allfdistribu.domain());
    DFieldSpXVx BB(allfdistribu.domain());
    DFieldSpXVx CC(allfdistribu.domain());
    compute_matrix_coeff(
            AA.span_view(),
            BB.span_view(),
            CC.span_view(),
            Dcoll.span_view(),
            Dcoll_staggered.span_view(),
            Nucoll.span_view(),
            dt);

    // rhs vector coefficient
    DFieldSpXVx RR(allfdistribu.domain());
    compute_rhs_vector(
            RR.span_view(),
            AA.span_cview(),
            BB.span_cview(),
            CC.span_cview(),
            allfdistribu.span_cview(),
            m_fthresh);


    for_each(
            policies::parallel_host,
            get_domain<IDimSp, IDimX>(allfdistribu),
            [&](IndexSpX const ispx) {
                Matrix_Banded matrix(get_domain<IDimVx>(allfdistribu).size(), 1, 1);
                fill_matrix_with_coeff(
                        matrix,
                        AA[ispx].span_cview(),
                        BB[ispx].span_cview(),
                        CC[ispx].span_cview());

                DSpan1D RR_Span1D(
                        RR[ispx].allocation_mdspan().data(),
                        get_domain<IDimVx>(allfdistribu).size());
                matrix.factorize();
                matrix.solve_inplace(RR_Span1D);
                ddc::deepcopy(allfdistribu[ispx], RR[ispx]);
            });

    return allfdistribu;
}
