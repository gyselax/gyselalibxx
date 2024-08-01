
#include <ddc/ddc.hpp>

#include <gtest/gtest.h>

#include "Lagrange_interpolator.hpp"
#include "bsl_advection_vx.hpp"
#include "geometry.hpp"
#include "spline_interpolator.hpp"


CoordX const x_min(0);
CoordX const x_max(2 * M_PI);
IdxStepX const x_size(50);
CoordVx const vx_min(-10);
CoordVx const vx_max(10);
IdxStepVx const vx_size(100);
std::pair<IdxRange<GridX>, IdxRange<GridVx>> Init_idx_range_velocity_adv()
{
    ddc::init_discrete_space<BSplinesX>(x_min, x_max, x_size);
    ddc::init_discrete_space<BSplinesVx>(vx_min, vx_max, vx_size);
    ddc::init_discrete_space<GridX>(SplineInterpPointsX::get_sampling<GridX>());
    ddc::init_discrete_space<GridVx>(SplineInterpPointsVx::get_sampling<GridVx>());
    IdxRange<GridX> x_dom = SplineInterpPointsX::get_domain<GridX>();
    IdxRange<GridVx> vx_dom = SplineInterpPointsVx::get_domain<GridVx>();
    return {x_dom, vx_dom};
}


template <class Geometry, class GridV>
double VelocityAdvection(
        IAdvectionVelocity<Geometry, GridV> const& advection_vx,
        IdxRange<GridX> x_dom,
        IdxRange<GridVx> vx_dom)
{
    //kinetic species
    IdxStepSp const nb_species(2);
    IdxRangeSp const dom_allsp(IdxSp(0), nb_species);
    IdxSp const i_elec = dom_allsp.front();
    IdxSp const i_ion = dom_allsp.back();
    //Mesh Initialization
    IdxRangeSpXVx const meshSpXVx(dom_allsp, x_dom, vx_dom);
    IdxRangeX const gridx = ddc::select<GridX>(meshSpXVx);
    IdxRangeVx const gridvx = ddc::select<GridVx>(meshSpXVx);
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

    // Initialization Species index range
    ddc::init_discrete_space<Species>(std::move(charges_host), std::move(masses_host));

    // Initialization of the distribution function
    host_t<DFieldMemSpXVx> allfdistribu_host(meshSpXVx);
    IdxRangeXVx const gridxvx = get_idx_range<GridX, GridVx>(allfdistribu_host);
    ddc::for_each(meshSpXVx, [&](IdxSpXVx const ispxvx) {
        IdxVx const ivx = ddc::select<GridVx>(ispxvx);
        allfdistribu_host(ispxvx) = exp(-0.5 * ddc::coordinate(ivx) * ddc::coordinate(ivx));
    });
    // Initialization of transport coefficient
    host_t<DFieldMemX> adv_speed_host(gridx);
    ddc::for_each(gridx, [&](IdxX const ix) { adv_speed_host(ix) = ddc::distance_at_right(ix); });
    double const timestep = .1;
    std::vector<double> Error;
    Error.reserve(allfdistribu_host.size());
    DFieldMemSpXVx allfdistribu(meshSpXVx);
    DFieldMemX electric_field(gridx);

    ddc::parallel_deepcopy(allfdistribu, allfdistribu_host);
    ddc::parallel_deepcopy(electric_field, adv_speed_host);
    advection_vx(allfdistribu, electric_field, timestep);
    ddc::parallel_deepcopy(allfdistribu_host, allfdistribu);

    double const m_advection_error = ddc::transform_reduce(
            meshSpXVx,
            0.0,
            ddc::reducer::max<double>(),
            [&](IdxSpXVx const ispxvx) {
                IdxSp const isp = ddc::select<Species>(ispxvx);
                IdxX const ix = ddc::select<GridX>(ispxvx);
                IdxVx const ivx = ddc::select<GridVx>(ispxvx);
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
    auto [x_dom, vx_dom] = Init_idx_range_velocity_adv();
    IdxStepVx static constexpr n_ghosts_vx {0};
    LagrangeInterpolator<GridVx, BCond::DIRICHLET, BCond::DIRICHLET, GridX, GridVx> const
            lagrange_vx_non_preallocatable_interpolator(3, n_ghosts_vx);
    PreallocatableLagrangeInterpolator<
            GridVx,
            BCond::DIRICHLET,
            BCond::DIRICHLET,
            GridX,
            GridVx> const lagrange_vx_interpolator(lagrange_vx_non_preallocatable_interpolator);
    BslAdvectionVelocity<GeometryXVx, GridVx> const lag_advection_vx(lagrange_vx_interpolator);
    double const err = VelocityAdvection<GeometryXVx, GridVx>(lag_advection_vx, x_dom, vx_dom);
    EXPECT_LE(err, 1e-3);
}

TEST(VelocityAdvection, SplineBatched)
{
    auto [x_dom, vx_dom] = Init_idx_range_velocity_adv();
    IdxRangeXVx meshXVx(x_dom, vx_dom);

    SplineVxBuilder const builder_vx(meshXVx);
    ddc::ConstantExtrapolationRule<Vx> bv_v_min(vx_min);
    ddc::ConstantExtrapolationRule<Vx> bv_v_max(vx_max);
    SplineVxEvaluator const spline_vx_evaluator(bv_v_min, bv_v_max);
    PreallocatableSplineInterpolator const spline_vx_interpolator(builder_vx, spline_vx_evaluator);
    BslAdvectionVelocity<GeometryXVx, GridVx> const spline_advection_vx(spline_vx_interpolator);
    double const err = VelocityAdvection<GeometryXVx, GridVx>(spline_advection_vx, x_dom, vx_dom);
    EXPECT_LE(err, 1e-5);
}
