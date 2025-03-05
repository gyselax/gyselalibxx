
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
    IdxRange<GridX> idx_range_x = SplineInterpPointsX::get_domain<GridX>();
    IdxRange<GridVx> idx_range_vx = SplineInterpPointsVx::get_domain<GridVx>();
    return {idx_range_x, idx_range_vx};
}


template <class Geometry, class GridV>
double VelocityAdvection(
        IAdvectionVelocity<Geometry, GridV> const& advection_vx,
        IdxRange<GridX> idx_range_x,
        IdxRange<GridVx> idx_range_vx)
{
    //kinetic species
    IdxStepSp const nb_species(2);
    IdxRangeSp const idx_range_allsp(IdxSp(0), nb_species);
    IdxSp const i_elec = idx_range_allsp.front();
    IdxSp const i_ion = idx_range_allsp.back();
    //Mesh Initialisation
    IdxRangeSpXVx const meshSpXVx(idx_range_allsp, idx_range_x, idx_range_vx);
    IdxRangeX const gridx = ddc::select<GridX>(meshSpXVx);
    // Charge Initialisation
    host_t<DFieldMemSp> masses_host(idx_range_allsp);
    host_t<DFieldMemSp> charges_host(idx_range_allsp);
    host_t<DFieldMemSp> init_perturb_amplitude_host(idx_range_allsp);
    // Mass ratio is fixed to one
    masses_host(i_elec) = 1.;
    init_perturb_amplitude_host(i_elec) = 0.;
    charges_host(i_elec) = -1.;

    masses_host(i_ion) = 1.;
    init_perturb_amplitude_host(i_ion) = 0.;
    charges_host(i_ion) = 1.;

    // Initialisation Species index range
    ddc::init_discrete_space<Species>(std::move(charges_host), std::move(masses_host));

    // Initialisation of the distribution function
    host_t<DFieldMemSpXVx> allfdistribu_host(meshSpXVx);
    ddc::for_each(meshSpXVx, [&](IdxSpXVx const ispxvx) {
        IdxVx const ivx = ddc::select<GridVx>(ispxvx);
        allfdistribu_host(ispxvx) = exp(-0.5 * ddc::coordinate(ivx) * ddc::coordinate(ivx));
    });
    // Initialisation of transport coefficient
    host_t<DFieldMemX> adv_speed_host(gridx);
    ddc::for_each(gridx, [&](IdxX const ix) { adv_speed_host(ix) = ddc::distance_at_right(ix); });
    double const timestep = .1;
    std::vector<double> Error;
    Error.reserve(allfdistribu_host.size());
    DFieldMemSpXVx allfdistribu(meshSpXVx);
    DFieldMemX electric_field(gridx);

    ddc::parallel_deepcopy(allfdistribu, allfdistribu_host);
    ddc::parallel_deepcopy(electric_field, adv_speed_host);
    advection_vx(get_field(allfdistribu), get_const_field(electric_field), timestep);
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
    auto [idx_range_x, idx_range_vx] = Init_idx_range_velocity_adv();
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
    double const err
            = VelocityAdvection<GeometryXVx, GridVx>(lag_advection_vx, idx_range_x, idx_range_vx);
    EXPECT_LE(err, 1e-3);
}

TEST(VelocityAdvection, SplineBatched)
{
    auto [idx_range_x, idx_range_vx] = Init_idx_range_velocity_adv();
    IdxRangeXVx meshXVx(idx_range_x, idx_range_vx);

    SplineVxBuilder const builder_vx(meshXVx);
    ddc::ConstantExtrapolationRule<Vx> bv_v_min(vx_min);
    ddc::ConstantExtrapolationRule<Vx> bv_v_max(vx_max);
    SplineVxEvaluator const spline_vx_evaluator(bv_v_min, bv_v_max);
    PreallocatableSplineInterpolator const spline_vx_interpolator(builder_vx, spline_vx_evaluator);
    BslAdvectionVelocity<GeometryXVx, GridVx> const spline_advection_vx(spline_vx_interpolator);
    double const err = VelocityAdvection<
            GeometryXVx,
            GridVx>(spline_advection_vx, idx_range_x, idx_range_vx);
    EXPECT_LE(err, 1e-5);
}
