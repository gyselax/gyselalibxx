#include <ddc/ddc.hpp>

#include <gtest/gtest.h>

#include "bsl_advection_x.hpp"
#include "geometry.hpp"
#include "spline_interpolator.hpp"



CoordX const x_min(-M_PI);
CoordX const x_max(M_PI);
IVectX const x_size(100);
CoordVx const vx_min(-6);
CoordVx const vx_max(6);
IVectVx const vx_size(50);

std::pair<ddc::DiscreteDomain<IDimX>, ddc::DiscreteDomain<IDimVx>> Init_domain_spatial_adv()
{
    ddc::init_discrete_space<BSplinesX>(x_min, x_max, x_size);
    ddc::init_discrete_space<BSplinesVx>(vx_min, vx_max, vx_size);
    ddc::init_discrete_space<IDimX>(SplineInterpPointsX::get_sampling<IDimX>());
    ddc::init_discrete_space<IDimVx>(SplineInterpPointsVx::get_sampling<IDimVx>());
    ddc::DiscreteDomain<IDimX> x_dom = SplineInterpPointsX::get_domain<IDimX>();
    ddc::DiscreteDomain<IDimVx> vx_dom = SplineInterpPointsVx::get_domain<IDimVx>();

    return {x_dom, vx_dom};
}


template <class Geometry, class DDimX>
double SpatialAdvection(
        IAdvectionSpatial<Geometry, DDimX> const& advection_x,
        ddc::DiscreteDomain<IDimX> x_dom,
        ddc::DiscreteDomain<IDimVx> vx_dom)
{
    //kinetic species
    IVectSp const nb_kinsp(1);
    IVectSp const nb_species(2);
    IDomainSp const dom_allsp(IndexSp(0), nb_species);
    IndexSp const i_elec = dom_allsp.front();
    IndexSp const i_ion = dom_allsp.back();
    //Mesh Initialization
    IDomainSpXVx const meshSpXVx(dom_allsp, x_dom, vx_dom);
    // Charge Initialization
    host_t<DFieldSp> masses_host(dom_allsp);
    host_t<FieldSp<int>> charges_host(dom_allsp);

    masses_host(i_elec) = 1;
    masses_host(i_ion) = 1;
    charges_host(i_elec) = -1;
    charges_host(i_ion) = 1;
    // Initialization Species domain
    ddc::init_discrete_space<IDimSp>(std::move(charges_host), std::move(masses_host));
    // Initialization of the distribution function
    host_t<DFieldSpXVx> allfdistribu_host(meshSpXVx);
    ddc::for_each(meshSpXVx, [&](IndexSpXVx const ispxvx) {
        IndexX const ix = ddc::select<IDimX>(ispxvx);
        allfdistribu_host(ispxvx) = cos(ddc::coordinate(ix));
    });

    double const timestep = .1;
    std::vector<double> Error;
    Error.reserve(allfdistribu_host.size());

    DFieldSpXVx allfdistribu(meshSpXVx);

    ddc::parallel_deepcopy(allfdistribu, allfdistribu_host);
    advection_x(allfdistribu, timestep);
    ddc::parallel_deepcopy(allfdistribu_host, allfdistribu);

    ddc::for_each(meshSpXVx, [&](IndexSpXVx const ispxvx) {
        IndexX const ix = ddc::select<IDimX>(ispxvx);
        IndexVx const ivx = ddc::select<IDimVx>(ispxvx);
        double const err = std::abs(
                allfdistribu_host(ispxvx)
                - cos(ddc::coordinate(ix) - ddc::coordinate(ivx) * timestep));
        Error.push_back(err);
    });

    double const m_advection_error = ddc::transform_reduce(
            meshSpXVx,
            0.0,
            ddc::reducer::max<double>(),
            [&](IndexSpXVx const ispxvx) {
                IndexX const ix = ddc::select<IDimX>(ispxvx);
                IndexVx const ivx = ddc::select<IDimVx>(ispxvx);
                return std::abs(
                        allfdistribu_host(ispxvx)
                        - cos(ddc::coordinate(ix) - ddc::coordinate(ivx) * timestep));
            });
    ;
    return m_advection_error;
}

TEST(SpatialAdvection, SplineBatched)
{
    auto [x_dom, vx_dom] = Init_domain_spatial_adv();
    IDomainXVx meshXVx(x_dom, vx_dom);
    SplineXBuilder const builder_x(meshXVx);
    ddc::PeriodicExtrapolationRule<RDimX> bv_x_min;
    ddc::PeriodicExtrapolationRule<RDimX> bv_x_max;
    SplineXEvaluator const spline_x_evaluator(bv_x_min, bv_x_max);
    PreallocatableSplineInterpolator const spline_x_interpolator(builder_x, spline_x_evaluator);
    BslAdvectionSpatial<GeometryXVx, IDimX> const spline_advection_x(spline_x_interpolator);
    double const err = SpatialAdvection<GeometryXVx, IDimX>(spline_advection_x, x_dom, vx_dom);
    EXPECT_LE(err, 1.e-6);
}