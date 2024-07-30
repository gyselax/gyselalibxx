
#include <ddc/ddc.hpp>

#include <gtest/gtest.h>

#include "Lagrange_interpolator.hpp"
#include "bsl_advection_vx.hpp"
#include "geometry.hpp"
#include "spline_interpolator.hpp"


CoordX const x_min(0);
CoordX const x_max(2 * M_PI);
IVectX const x_size(50);
CoordVx const vx_min(-10);
CoordVx const vx_max(10);
IVectVx const vx_size(100);
std::pair<ddc::DiscreteDomain<IDimX>, ddc::DiscreteDomain<IDimVx>> Init_domain_velocity_adv()
{
    ddc::init_discrete_space<BSplinesX>(x_min, x_max, x_size);
    ddc::init_discrete_space<BSplinesVx>(vx_min, vx_max, vx_size);
    ddc::init_discrete_space<IDimX>(SplineInterpPointsX::get_sampling<IDimX>());
    ddc::init_discrete_space<IDimVx>(SplineInterpPointsVx::get_sampling<IDimVx>());
    ddc::DiscreteDomain<IDimX> x_dom = SplineInterpPointsX::get_domain<IDimX>();
    ddc::DiscreteDomain<IDimVx> vx_dom = SplineInterpPointsVx::get_domain<IDimVx>();
    return {x_dom, vx_dom};
}


template <class Geometry, class DDimV>
double VelocityAdvection(
        IAdvectionVelocity<Geometry, DDimV> const& advection_vx,
        ddc::DiscreteDomain<IDimX> x_dom,
        ddc::DiscreteDomain<IDimVx> vx_dom)
{
    //kinetic species
    IdxStepSp const nb_species(2);
    IdxRangeSp const dom_allsp(IdxSp(0), nb_species);
    IdxSp const i_elec = dom_allsp.front();
    IdxSp const i_ion = dom_allsp.back();
    //Mesh Initialization
    IDomainSpXVx const meshSpXVx(dom_allsp, x_dom, vx_dom);
    IDomainX const gridx = ddc::select<IDimX>(meshSpXVx);
    IDomainVx const gridvx = ddc::select<IDimVx>(meshSpXVx);
    // Charge Initialization
    host_t<DFieldMemSp> masses_host(dom_allsp);
    host_t<DFieldMemSp> charges_host(dom_allsp);
    host_t<DFieldMemSp> init_perturb_amplitude_host(dom_allsp);
    // Mass ratio is fixed to one
    masses_host(i_elec) = 1.;
    init_perturb_amplitude_host(i_elec) = 0.;
    charges_host(i_elec) = -1.;

    masses_host(i_ion) = 1.;
    init_perturb_amplitude_host(i_ion) = 0.;
    charges_host(i_ion) = 1.;

    // Initialization Species domain
    ddc::init_discrete_space<Species>(std::move(charges_host), std::move(masses_host));

    // Initialization of the distribution function
    host_t<DFieldSpXVx> allfdistribu_host(meshSpXVx);
    IDomainXVx const gridxvx = allfdistribu_host.domain<IDimX, IDimVx>();
    ddc::for_each(meshSpXVx, [&](IndexSpXVx const ispxvx) {
        IndexVx const ivx = ddc::select<IDimVx>(ispxvx);
        allfdistribu_host(ispxvx) = exp(-0.5 * ddc::coordinate(ivx) * ddc::coordinate(ivx));
    });
    // Initialization of transport coefficient
    host_t<DFieldX> adv_speed_host(gridx);
    ddc::for_each(gridx, [&](IndexX const ix) { adv_speed_host(ix) = ddc::distance_at_right(ix); });
    double const timestep = .1;
    std::vector<double> Error;
    Error.reserve(allfdistribu_host.size());
    DFieldSpXVx allfdistribu(meshSpXVx);
    DFieldX electric_field(gridx);

    ddc::parallel_deepcopy(allfdistribu, allfdistribu_host);
    ddc::parallel_deepcopy(electric_field, adv_speed_host);
    advection_vx(allfdistribu, electric_field, timestep);
    ddc::parallel_deepcopy(allfdistribu_host, allfdistribu);

    double const m_advection_error = ddc::transform_reduce(
            meshSpXVx,
            0.0,
            ddc::reducer::max<double>(),
            [&](IndexSpXVx const ispxvx) {
                IdxSp const isp = ddc::select<Species>(ispxvx);
                IndexX const ix = ddc::select<IDimX>(ispxvx);
                IndexVx const ivx = ddc::select<IDimVx>(ispxvx);
                return std::abs(
                        allfdistribu_host(ispxvx)
                        - exp(-0.5
                              * std::
                                      pow(ddc::coordinate(ivx)
                                                  - charge(isp) * adv_speed_host(ix) * timestep,
                                          2)));
            });
    return m_advection_error;
}

TEST(VelocityAdvection, BatchedLagrange)
{
    auto [x_dom, vx_dom] = Init_domain_velocity_adv();
    IVectVx static constexpr n_ghosts_vx {0};
    LagrangeInterpolator<IDimVx, BCond::DIRICHLET, BCond::DIRICHLET, IDimX, IDimVx> const
            lagrange_vx_non_preallocatable_interpolator(3, n_ghosts_vx);
    PreallocatableLagrangeInterpolator<
            IDimVx,
            BCond::DIRICHLET,
            BCond::DIRICHLET,
            IDimX,
            IDimVx> const lagrange_vx_interpolator(lagrange_vx_non_preallocatable_interpolator);
    BslAdvectionVelocity<GeometryXVx, IDimVx> const lag_advection_vx(lagrange_vx_interpolator);
    double const err = VelocityAdvection<GeometryXVx, IDimVx>(lag_advection_vx, x_dom, vx_dom);
    EXPECT_LE(err, 1e-3);
}

TEST(VelocityAdvection, SplineBatched)
{
    auto [x_dom, vx_dom] = Init_domain_velocity_adv();
    IDomainXVx meshXVx(x_dom, vx_dom);

    SplineVxBuilder const builder_vx(meshXVx);
    ddc::ConstantExtrapolationRule<RDimVx> bv_v_min(vx_min);
    ddc::ConstantExtrapolationRule<RDimVx> bv_v_max(vx_max);
    SplineVxEvaluator const spline_vx_evaluator(bv_v_min, bv_v_max);
    PreallocatableSplineInterpolator const spline_vx_interpolator(builder_vx, spline_vx_evaluator);
    BslAdvectionVelocity<GeometryXVx, IDimVx> const spline_advection_vx(spline_vx_interpolator);
    double const err = VelocityAdvection<GeometryXVx, IDimVx>(spline_advection_vx, x_dom, vx_dom);
    EXPECT_LE(err, 1e-5);
}