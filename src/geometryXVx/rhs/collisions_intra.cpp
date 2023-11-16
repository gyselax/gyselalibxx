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
 *     @f$Imean0=\langle Dcoll \rangle@f$ ;
 *     @f$Imean1=\langle v*Dcoll \rangle@f$ ;
 *     @f$Imean2=\langle v^2*Dcoll \rangle@f$ ;
 *     @f$Imean3=\langle d/dv(Dcoll) \rangle@f$
 *     @f$Imean4=\langle d/dv(v*Dcoll) \rangle@f$
 *  The brackets @f$\langle\cdot\rangle@f$ represent the integral in velocity: @f$\langle\cdot\rangle = \int \cdot dv@f$
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
    , m_nustar_profile(ddc::select<IDimSp, IDimX>(mesh))
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

    compute_nustar_profile(m_nustar_profile.span_view(), m_nustar0);
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
    ddc::for_each(AA.domain(), [&](IndexSpXVx const ispxvx) {
        IndexSp const isp = ddc::select<IDimSp>(ispxvx);
        IndexX const ix = ddc::select<IDimX>(ispxvx);

        IndexVx_ghosted ivx_ghosted(ddc::select<IDimVx>(ispxvx).uid() + 1);
        IndexVx_ghosted_staggered ivx_ghosted_staggered(ddc::select<IDimVx>(ispxvx).uid() + 1);
        IndexVx_ghosted ivx_next_ghosted(ivx_ghosted + 1);
        IndexVx_ghosted ivx_prev_ghosted(ivx_ghosted - 1);
        IndexVx_ghosted_staggered ivx_prev_ghosted_staggered(ivx_ghosted_staggered - 1);

        double const dv_i = ddc::coordinate(ivx_next_ghosted) - ddc::coordinate(ivx_ghosted);
        double const delta_i
                = dv_i / (ddc::coordinate(ivx_ghosted) - ddc::coordinate(ivx_prev_ghosted));

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
    ddc::for_each(RR.domain(), [&](IndexSpXVx const ispxvx) {
        IndexSp const isp = ddc::select<IDimSp>(ispxvx);
        IndexX const ix = ddc::select<IDimX>(ispxvx);
        IndexVx const ivx = ddc::select<IDimVx>(ispxvx);

        IndexVx const ivx_next = ivx + 1;
        IndexVx const ivx_prev = ivx - 1;

        if (ivx == ddc::get_domain<IDimVx>(AA).front()) {
            RR(isp, ix, ivx) = (2. - BB(isp, ix, ivx)) * allfdistribu(isp, ix, ivx)
                               + (-CC(isp, ix, ivx)) * allfdistribu(isp, ix, ivx_next)
                               - 2. * AA(isp, ix, ivx) * fthresh;

        } else if (ivx == ddc::get_domain<IDimVx>(AA).back()) {
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
    // density and temperature
    DFieldSpX density(ddc::get_domain<IDimSp, IDimX>(allfdistribu));
    DFieldSpX mean_velocity(ddc::get_domain<IDimSp, IDimX>(allfdistribu));
    DFieldSpX temperature(ddc::get_domain<IDimSp, IDimX>(allfdistribu));

    DFieldVx const quadrature_coeffs
            = trapezoid_quadrature_coefficients(ddc::get_domain<IDimVx>(allfdistribu));
    Quadrature<IDimVx> integrate(quadrature_coeffs);
    FluidMoments moments(integrate);

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
    DFieldSpX collfreq(ddc::get_domain<IDimSp, IDimX>(allfdistribu));
    compute_collfreq(
            collfreq.span_view(),
            m_nustar_profile.span_cview(),
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
    DFieldSpX Vcoll(ddc::get_domain<IDimSp, IDimX>(allfdistribu));
    DFieldSpX Tcoll(ddc::get_domain<IDimSp, IDimX>(allfdistribu));
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


    ddc::for_each(ddc::get_domain<IDimSp, IDimX>(allfdistribu), [&](IndexSpX const ispx) {
        Matrix_Banded matrix(ddc::get_domain<IDimVx>(allfdistribu).size(), 1, 1);
        fill_matrix_with_coeff(
                matrix,
                AA[ispx].span_cview(),
                BB[ispx].span_cview(),
                CC[ispx].span_cview());

        DSpan1D RR_Span1D(
                RR[ispx].allocation_mdspan().data_handle(),
                ddc::get_domain<IDimVx>(allfdistribu).size());
        matrix.factorize();
        matrix.solve_inplace(RR_Span1D);
        ddc::deepcopy(allfdistribu[ispx], RR[ispx]);
    });

    return allfdistribu;
}
