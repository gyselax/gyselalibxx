#include <ddc/ddc.hpp>

#include <gtest/gtest.h>

#include "bsl_advection_x.hpp"
#include "geometry.hpp"
#include "spline_interpolator.hpp"



CoordX const x_min(-M_PI);
CoordX const x_max(M_PI);
IdxStepX const x_size(100);
CoordVx const vx_min(-6);
CoordVx const vx_max(6);
IdxStepVx const vx_size(50);

std::pair<IdxRange<GridX>, IdxRange<GridVx>> Init_idx_range_spatial_adv()
{
    ddc::init_discrete_space<BSplinesX>(x_min, x_max, x_size);
    ddc::init_discrete_space<BSplinesVx>(vx_min, vx_max, vx_size);
    ddc::init_discrete_space<GridX>(SplineInterpPointsX::get_sampling<GridX>());
    ddc::init_discrete_space<GridVx>(SplineInterpPointsVx::get_sampling<GridVx>());
    IdxRange<GridX> x_dom = SplineInterpPointsX::get_domain<GridX>();
    IdxRange<GridVx> vx_dom = SplineInterpPointsVx::get_domain<GridVx>();

    return {x_dom, vx_dom};
}


template <class Geometry, class GridX>
double SpatialAdvection(
        IAdvectionSpatial<Geometry, GridX> const& advection_x,
        IdxRange<GridX> x_dom,
        IdxRange<GridVx> vx_dom)
{
    //kinetic species
    IdxStepSp const nb_kinsp(1);
    IdxStepSp const nb_species(2);
    IdxRangeSp const dom_allsp(IdxSp(0), nb_species);
    IdxSp const i_elec = dom_allsp.front();
    IdxSp const i_ion = dom_allsp.back();
    //Mesh Initialization
    IdxRangeSpXVx const meshSpXVx(dom_allsp, x_dom, vx_dom);
    // Charge Initialization
    host_t<DFieldMemSp> masses_host(dom_allsp);
    host_t<DFieldMemSp> charges_host(dom_allsp);

    masses_host(i_elec) = 1.;
    masses_host(i_ion) = 1.;
    charges_host(i_elec) = -1.;
    charges_host(i_ion) = 1.;
    // Initialization Species index range
    ddc::init_discrete_space<Species>(std::move(charges_host), std::move(masses_host));
    // Initialization of the distribution function
    host_t<DFieldMemSpXVx> allfdistribu_host(meshSpXVx);
    ddc::for_each(meshSpXVx, [&](IdxSpXVx const ispxvx) {
        IdxX const ix = ddc::select<GridX>(ispxvx);
        allfdistribu_host(ispxvx) = cos(ddc::coordinate(ix));
    });

    double const timestep = .1;
    std::vector<double> Error;
    Error.reserve(allfdistribu_host.size());

    DFieldMemSpXVx allfdistribu(meshSpXVx);

    ddc::parallel_deepcopy(allfdistribu, allfdistribu_host);
    advection_x(allfdistribu, timestep);
    ddc::parallel_deepcopy(allfdistribu_host, allfdistribu);

    ddc::for_each(meshSpXVx, [&](IdxSpXVx const ispxvx) {
        IdxX const ix = ddc::select<GridX>(ispxvx);
        IdxVx const ivx = ddc::select<GridVx>(ispxvx);
        double const err = std::abs(
                allfdistribu_host(ispxvx)
                - cos(ddc::coordinate(ix) - ddc::coordinate(ivx) * timestep));
        Error.push_back(err);
    });

    double const m_advection_error = ddc::transform_reduce(
            meshSpXVx,
            0.0,
            ddc::reducer::max<double>(),
            [&](IdxSpXVx const ispxvx) {
                IdxX const ix = ddc::select<GridX>(ispxvx);
                IdxVx const ivx = ddc::select<GridVx>(ispxvx);
                return std::abs(
                        allfdistribu_host(ispxvx)
                        - cos(ddc::coordinate(ix) - ddc::coordinate(ivx) * timestep));
            });
    ;
    return m_advection_error;
}

TEST(SpatialAdvection, SplineBatched)
{
    auto [x_dom, vx_dom] = Init_idx_range_spatial_adv();
    IdxRangeXVx meshXVx(x_dom, vx_dom);
    SplineXBuilder const builder_x(meshXVx);
    ddc::PeriodicExtrapolationRule<X> bv_x_min;
    ddc::PeriodicExtrapolationRule<X> bv_x_max;
    SplineXEvaluator const spline_x_evaluator(bv_x_min, bv_x_max);
    PreallocatableSplineInterpolator const spline_x_interpolator(builder_x, spline_x_evaluator);
    BslAdvectionSpatial<GeometryXVx, GridX> const spline_advection_x(spline_x_interpolator);
    double const err = SpatialAdvection<GeometryXVx, GridX>(spline_advection_x, x_dom, vx_dom);
    EXPECT_LE(err, 1.e-6);
}
