// SPDX-License-Identifier: MIT
#include <cmath>
#include <iostream>
#include <string>

#include <ddc/ddc.hpp>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <paraconf.h>
#include <pdi.h>

#include "ddc_alias_inline_functions.hpp"
#include "geometry.hpp"
#include "irighthandside.hpp"
#include "krook_source_adaptive.hpp"
#include "krook_source_constant.hpp"
#include "mask_tanh.hpp"
#include "maxwellianequilibrium.hpp"
#include "quadrature.hpp"
#include "species_info.hpp"
#include "splitrighthandsidesolver.hpp"
#include "trapezoid_quadrature.hpp"

TEST(KrookSource, Adaptive)
{
    CoordX const x_min(0.0);
    CoordX const x_max(1.0);
    IdxStepX const x_size(10);

    CoordVx const vx_min(-6);
    CoordVx const vx_max(6);
    IdxStepVx const vx_size(50);

    IdxStepSp const nb_kinspecies(2);

    IdxRangeSp const idx_range_sp(IdxSp(0), nb_kinspecies);

    PC_tree_t conf_pdi = PC_parse_string("");
    PDI_init(conf_pdi);

    // Creating mesh & supports
    ddc::init_discrete_space<BSplinesX>(x_min, x_max, x_size);
    ddc::init_discrete_space<BSplinesVx>(vx_min, vx_max, vx_size);

    ddc::init_discrete_space<GridX>(SplineInterpPointsX::get_sampling<GridX>());
    ddc::init_discrete_space<GridVx>(SplineInterpPointsVx::get_sampling<GridVx>());

    IdxRangeX gridx(SplineInterpPointsX::get_domain<GridX>());
    IdxRangeVx gridvx(SplineInterpPointsVx::get_domain<GridVx>());

    SplineXBuilder_1d const builder_x(gridx);
    SplineVxBuilder_1d const builder_vx(gridvx);

    IdxRangeSp const gridsp = idx_range_sp;
    IdxRangeSpXVx const mesh(gridsp, gridx, gridvx);

    DFieldMemVx const quadrature_coeffs_vx(
            trapezoid_quadrature_coefficients<Kokkos::DefaultExecutionSpace>(gridvx));
    Quadrature<IdxRangeVx> const integrate_v(get_const_field(quadrature_coeffs_vx));

    host_t<DFieldMemSp> charges(idx_range_sp);
    host_t<DFieldMemSp> masses(idx_range_sp);
    IdxSp my_iion(idx_range_sp.front());
    IdxSp my_ielec(idx_range_sp.back());
    charges(my_iion) = 1.;
    charges(my_ielec) = -1.;
    ddc::for_each(idx_range_sp, [&](IdxSp const isp) { masses(isp) = 1.0; });

    // Initialization of the distribution function
    ddc::init_discrete_space<Species>(std::move(charges), std::move(masses));

    double const extent = 0.5;
    double const stiffness = 0.01;
    double const amplitude = 0.1;
    double const density_target = 0.5;
    double const temperature_target = 0.5;

    KrookSourceAdaptive const rhs_krook(
            gridx,
            gridvx,
            RhsType::Sink,
            extent,
            stiffness,
            amplitude,
            density_target,
            temperature_target);

    // Initialization of the distribution function : maxwellian
    double const density_init_ion = 1.;
    double const density_init_elec = 2.;
    double const temperature_init = 1.;
    DFieldMemSpXVx allfdistribu(mesh);
    ddc::for_each(ddc::select<Species, GridX>(mesh), [&](IdxSpX const ispx) {
        DFieldMemVx finit(gridvx);
        if (charge(ddc::select<Species>(ispx)) >= 0.) {
            MaxwellianEquilibrium::
                    compute_maxwellian(finit, density_init_ion, temperature_init, 0.);
        } else {
            MaxwellianEquilibrium::
                    compute_maxwellian(finit, density_init_elec, temperature_init, 0.);
        }
        auto finit_host = ddc::create_mirror_view_and_copy(get_field(finit));
        ddc::parallel_deepcopy(allfdistribu[ispx], finit_host);
    });

    // error with a given deltat
    double const deltat = 0.1;
    rhs_krook(allfdistribu, deltat);
    auto allfdistribu_host = ddc::create_mirror_view_and_copy(get_field(allfdistribu));

    host_t<DFieldMemSpX> densities(get_idx_range<Species, GridX>(allfdistribu));
    ddc::for_each(get_idx_range<Species, GridX>(allfdistribu), [&](IdxSpX const ispx) {
        densities(ispx) = integrate_v(Kokkos::DefaultExecutionSpace(), allfdistribu[ispx]);
    });

    // the charge should be conserved by the operator
    host_t<DFieldMemX> error(get_idx_range<GridX>(allfdistribu_host));
    ddc::for_each(get_idx_range<GridX>(allfdistribu), [&](IdxX const ix) {
        error(ix) = std::fabs(
                charge(my_iion) * (densities(my_iion, ix) - density_init_ion)
                + charge(my_ielec) * (densities(my_ielec, ix) - density_init_elec));
    });

    // reinitialization of the distribution function
    ddc::for_each(get_idx_range<Species, GridX>(allfdistribu), [&](IdxSpX const ispx) {
        DFieldMemVx finit(gridvx);
        if (charge(ddc::select<Species>(ispx)) >= 0.) {
            MaxwellianEquilibrium::
                    compute_maxwellian(finit, density_init_ion, temperature_init, 0.);
        } else {
            MaxwellianEquilibrium::
                    compute_maxwellian(finit, density_init_elec, temperature_init, 0.);
        }
        auto finit_host = ddc::create_mirror_view_and_copy(get_field(finit));
        ddc::parallel_deepcopy(allfdistribu[ispx], finit_host);
    });

    // error with a deltat 10 times smaller
    rhs_krook(allfdistribu, 0.01);
    ddc::parallel_deepcopy(allfdistribu_host, allfdistribu);
    ddc::for_each(get_idx_range<Species, GridX>(allfdistribu), [&](IdxSpX const ispx) {
        densities(ispx) = integrate_v(Kokkos::DefaultExecutionSpace(), allfdistribu[ispx]);
    });

    // the rk2 scheme used in the krook operator should be of order 2
    // hence the error should be divided by at least 100 when dt is divided by 10
    double const order = 2;
    ddc::for_each(get_idx_range<GridX>(allfdistribu), [&](IdxX const ix) {
        double const error_smalldt = std::fabs(
                charge(my_iion) * (densities(my_iion, ix) - density_init_ion)
                + charge(my_ielec) * (densities(my_ielec, ix) - density_init_elec));
        std::cout << "ix " << ix << "error " << error(ix) << " error_smalldt " << error_smalldt
                  << std::endl;
        EXPECT_LE(error_smalldt, error(ix) / std::pow(10, order));
    });

    PC_tree_destroy(&conf_pdi);
    PDI_finalize();
}

TEST(KrookSource, Constant)
{
    CoordX const x_min(0.0);
    CoordX const x_max(1.0);
    IdxStepX const x_size(10);

    CoordVx const vx_min(-6);
    CoordVx const vx_max(6);
    IdxStepVx const vx_size(10);

    IdxStepSp const nb_kinspecies(2);

    IdxRangeSp const idx_range_sp(IdxSp(0), nb_kinspecies);

    PC_tree_t conf_pdi = PC_parse_string("");
    PDI_init(conf_pdi);

    // Creating mesh & supports
    ddc::init_discrete_space<BSplinesX>(x_min, x_max, x_size);

    ddc::init_discrete_space<BSplinesVx>(vx_min, vx_max, vx_size);

    ddc::init_discrete_space<GridX>(SplineInterpPointsX::get_sampling<GridX>());
    ddc::init_discrete_space<GridVx>(SplineInterpPointsVx::get_sampling<GridVx>());

    IdxRangeX gridx(SplineInterpPointsX::get_domain<GridX>());
    IdxRangeVx gridvx(SplineInterpPointsVx::get_domain<GridVx>());

    SplineXBuilder_1d const builder_x(gridx);
    SplineVxBuilder_1d const builder_vx(gridvx);

    IdxRangeSp const gridsp = idx_range_sp;

    IdxRangeSpXVx const mesh(gridsp, gridx, gridvx);

    host_t<DFieldMemSp> charges(idx_range_sp);
    host_t<DFieldMemSp> masses(idx_range_sp);
    charges(idx_range_sp.front()) = 1.;
    charges(idx_range_sp.back()) = -1.;
    ddc::for_each(idx_range_sp, [&](IdxSp const isp) { masses(isp) = 1.0; });

    // Initialization of the distribution function
    ddc::init_discrete_space<Species>(std::move(charges), std::move(masses));

    double const extent = 0.25;
    double const stiffness = 0.01;
    double const amplitude = 0.1;
    double const density_target = 0.5;
    double const temperature_target = 0.5;

    KrookSourceConstant const rhs_krook(
            gridx,
            gridvx,
            RhsType::Sink,
            extent,
            stiffness,
            amplitude,
            density_target,
            temperature_target);

    // compute the krook mask (spatial extent)
    host_t<DFieldMemX> mask = mask_tanh(gridx, extent, stiffness, MaskType::Inverted, false);

    // simulation
    double const deltat = 1.;

    // Initialization of the distribution function : maxwellian
    double const density_init = 1.;
    double const temperature_init = 1.;
    DFieldMemVx finit(gridvx);
    MaxwellianEquilibrium::compute_maxwellian(finit, density_init, temperature_init, 0.);
    auto finit_host = ddc::create_mirror_view_and_copy(get_field(finit));
    DFieldMemSpXVx allfdistribu(mesh);
    ddc::for_each(ddc::select<Species, GridX>(mesh), [&](IdxSpX const ispx) {
        ddc::parallel_deepcopy(allfdistribu[ispx], finit_host);
    });

    int const nbsteps = 100;
    for (int iter = 0; iter < nbsteps; ++iter) {
        rhs_krook(allfdistribu, deltat);
    };
    auto allfdistribu_host = ddc::create_mirror_view_and_copy(get_field(allfdistribu));

    // tests if distribution function matches theoretical prediction
    DFieldMemVx ftarget(gridvx);
    MaxwellianEquilibrium::compute_maxwellian(ftarget, density_target, temperature_target, 0.);
    auto ftarget_host = ddc::create_mirror_view_and_copy(get_field(ftarget));

    ddc::for_each(get_idx_range(allfdistribu), [&](IdxSpXVx const ispxvx) {
        // predicted distribution function value
        double const allfdistribu_pred
                = ftarget_host(ddc::select<GridVx>(ispxvx))
                  + (finit_host(ddc::select<GridVx>(ispxvx))
                     - ftarget_host(ddc::select<GridVx>(ispxvx)))
                            * std::exp(
                                    -amplitude * mask(ddc::select<GridX>(ispxvx)) * deltat
                                    * nbsteps);
        double const error = std::fabs(allfdistribu_host(ispxvx) - allfdistribu_pred);

        EXPECT_LE(error, 1e-13);
    });

    PC_tree_destroy(&conf_pdi);
    PDI_finalize();
}
