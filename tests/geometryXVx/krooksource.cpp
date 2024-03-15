// SPDX-License-Identifier: MIT

#include <cmath>
#include <iostream>
#include <string>

#include <ddc/ddc.hpp>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <femperiodicpoissonsolver.hpp>
#include <geometry.hpp>
#include <irighthandside.hpp>
#include <krook_source_adaptive.hpp>
#include <krook_source_constant.hpp>
#include <mask_tanh.hpp>
#include <maxwellianequilibrium.hpp>
#include <paraconf.h>
#include <pdi.h>
#include <quadrature.hpp>
#include <species_info.hpp>
#include <splitrighthandsidesolver.hpp>
#include <trapezoid_quadrature.hpp>

TEST(KrookSource, Adaptive)
{
    CoordX const x_min(0.0);
    CoordX const x_max(1.0);
    IVectX const x_size(10);

    CoordVx const vx_min(-6);
    CoordVx const vx_max(6);
    IVectVx const vx_size(50);

    IVectSp const nb_kinspecies(2);

    IDomainSp const dom_sp(IndexSp(0), nb_kinspecies);

    PC_tree_t conf_pdi = PC_parse_string("");
    PDI_init(conf_pdi);

    // Creating mesh & supports
    ddc::init_discrete_space<BSplinesX>(x_min, x_max, x_size);
    ddc::init_discrete_space<BSplinesVx>(vx_min, vx_max, vx_size);

    ddc::init_discrete_space<IDimX>(SplineInterpPointsX::get_sampling());
    ddc::init_discrete_space<IDimVx>(SplineInterpPointsVx::get_sampling());

    IDomainX gridx(SplineInterpPointsX::get_domain());
    IDomainVx gridvx(SplineInterpPointsVx::get_domain());

    SplineXBuilder_1d const builder_x(gridx);
    SplineVxBuilder_1d const builder_vx(gridvx);

    IDomainSp const gridsp = dom_sp;
    IDomainSpXVx const mesh(gridsp, gridx, gridvx);

    host_t<DFieldX> const quadrature_coeffs_x = trapezoid_quadrature_coefficients(gridx);
    host_t<DFieldVx> const quadrature_coeffs_vx = trapezoid_quadrature_coefficients(gridvx);
    Quadrature<IDimX> const integrate_x(quadrature_coeffs_x);
    Quadrature<IDimVx> const integrate_v(quadrature_coeffs_vx);

    host_t<FieldSp<int>> charges(dom_sp);
    host_t<DFieldSp> masses(dom_sp);
    IndexSp my_iion(dom_sp.front());
    IndexSp my_ielec(dom_sp.back());
    charges(my_iion) = 1;
    charges(my_ielec) = -1;
    ddc::for_each(dom_sp, [&](IndexSp const isp) { masses(isp) = 1.0; });

    // Initialization of the distribution function
    ddc::init_discrete_space<IDimSp>(std::move(charges), std::move(masses));

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
    DFieldSpXVx allfdistribu(mesh);
    ddc::for_each(ddc::select<IDimSp, IDimX>(mesh), [&](IndexSpX const ispx) {
        DFieldVx finit(gridvx);
        if (charge(ddc::select<IDimSp>(ispx)) >= 0.) {
            MaxwellianEquilibrium::
                    compute_maxwellian(finit, density_init_ion, temperature_init, 0.);
        } else {
            MaxwellianEquilibrium::
                    compute_maxwellian(finit, density_init_elec, temperature_init, 0.);
        }
        auto finit_host = ddc::create_mirror_view_and_copy(finit.span_view());
        ddc::parallel_deepcopy(allfdistribu[ispx], finit_host);
    });

    // error with a given deltat
    double const deltat = 0.1;
    rhs_krook(allfdistribu, deltat);
    auto allfdistribu_host = ddc::create_mirror_view_and_copy(allfdistribu.span_view());

    host_t<DFieldSpX> densities(ddc::get_domain<IDimSp, IDimX>(allfdistribu));
    ddc::for_each(ddc::get_domain<IDimSp, IDimX>(allfdistribu), [&](IndexSpX const ispx) {
        densities(ispx) = integrate_v(allfdistribu_host[ispx]);
    });

    // the charge should be conserved by the operator
    host_t<DFieldX> error(ddc::get_domain<IDimX>(allfdistribu_host));
    ddc::for_each(ddc::get_domain<IDimX>(allfdistribu), [&](IndexX const ix) {
        error(ix) = std::fabs(
                charge(my_iion) * (densities(my_iion, ix) - density_init_ion)
                + charge(my_ielec) * (densities(my_ielec, ix) - density_init_elec));
    });

    // reinitialization of the distribution function
    ddc::for_each(ddc::get_domain<IDimSp, IDimX>(allfdistribu), [&](IndexSpX const ispx) {
        DFieldVx finit(gridvx);
        if (charge(ddc::select<IDimSp>(ispx)) >= 0.) {
            MaxwellianEquilibrium::
                    compute_maxwellian(finit, density_init_ion, temperature_init, 0.);
        } else {
            MaxwellianEquilibrium::
                    compute_maxwellian(finit, density_init_elec, temperature_init, 0.);
        }
        auto finit_host = ddc::create_mirror_view_and_copy(finit.span_view());
        ddc::parallel_deepcopy(allfdistribu[ispx], finit_host);
    });

    // error with a deltat 10 times smaller
    rhs_krook(allfdistribu, 0.01);
    ddc::parallel_deepcopy(allfdistribu_host, allfdistribu);
    ddc::for_each(ddc::get_domain<IDimSp, IDimX>(allfdistribu), [&](IndexSpX const ispx) {
        densities(ispx) = integrate_v(allfdistribu_host[ispx]);
    });

    // the rk2 scheme used in the krook operator should be of order 2
    // hence the error should be divided by at least 100 when dt is divided by 10
    double const order = 2;
    ddc::for_each(ddc::get_domain<IDimX>(allfdistribu), [&](IndexX const ix) {
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
    IVectX const x_size(10);

    CoordVx const vx_min(-6);
    CoordVx const vx_max(6);
    IVectVx const vx_size(10);

    IVectSp const nb_kinspecies(2);

    IDomainSp const dom_sp(IndexSp(0), nb_kinspecies);

    PC_tree_t conf_pdi = PC_parse_string("");
    PDI_init(conf_pdi);

    // Creating mesh & supports
    ddc::init_discrete_space<BSplinesX>(x_min, x_max, x_size);

    ddc::init_discrete_space<BSplinesVx>(vx_min, vx_max, vx_size);

    ddc::init_discrete_space<IDimX>(SplineInterpPointsX::get_sampling());
    ddc::init_discrete_space<IDimVx>(SplineInterpPointsVx::get_sampling());

    IDomainX gridx(SplineInterpPointsX::get_domain());
    IDomainVx gridvx(SplineInterpPointsVx::get_domain());

    SplineXBuilder_1d const builder_x(gridx);
    SplineVxBuilder_1d const builder_vx(gridvx);

    IDomainSp const gridsp = dom_sp;

    IDomainSpXVx const mesh(gridsp, gridx, gridvx);

    host_t<FieldSp<int>> charges(dom_sp);
    host_t<DFieldSp> masses(dom_sp);
    charges(dom_sp.front()) = 1;
    charges(dom_sp.back()) = -1;
    ddc::for_each(dom_sp, [&](IndexSp const isp) { masses(isp) = 1.0; });

    // Initialization of the distribution function
    ddc::init_discrete_space<IDimSp>(std::move(charges), std::move(masses));

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
    host_t<DFieldX> mask = mask_tanh(gridx, extent, stiffness, MaskType::Inverted, false);

    // simulation
    double const deltat = 1.;

    // Initialization of the distribution function : maxwellian
    double const density_init = 1.;
    double const temperature_init = 1.;
    DFieldVx finit(gridvx);
    MaxwellianEquilibrium::compute_maxwellian(finit, density_init, temperature_init, 0.);
    auto finit_host = ddc::create_mirror_view_and_copy(finit.span_view());
    DFieldSpXVx allfdistribu(mesh);
    ddc::for_each(ddc::select<IDimSp, IDimX>(mesh), [&](IndexSpX const ispx) {
        ddc::parallel_deepcopy(allfdistribu[ispx], finit_host);
    });

    int const nbsteps = 100;
    for (int iter = 0; iter < nbsteps; ++iter) {
        rhs_krook(allfdistribu, deltat);
    };
    auto allfdistribu_host = ddc::create_mirror_view_and_copy(allfdistribu.span_view());

    // tests if distribution function matches theoretical prediction
    DFieldVx ftarget(gridvx);
    MaxwellianEquilibrium::compute_maxwellian(ftarget, density_target, temperature_target, 0.);
    auto ftarget_host = ddc::create_mirror_view_and_copy(ftarget.span_view());

    ddc::for_each(allfdistribu.domain(), [&](IndexSpXVx const ispxvx) {
        // predicted distribution function value
        double const allfdistribu_pred
                = ftarget_host(ddc::select<IDimVx>(ispxvx))
                  + (finit_host(ddc::select<IDimVx>(ispxvx))
                     - ftarget_host(ddc::select<IDimVx>(ispxvx)))
                            * std::exp(
                                    -amplitude * mask(ddc::select<IDimX>(ispxvx)) * deltat
                                    * nbsteps);
        double const error = std::fabs(allfdistribu_host(ispxvx) - allfdistribu_pred);

        EXPECT_LE(error, 1e-13);
    });

    PC_tree_destroy(&conf_pdi);
    PDI_finalize();
}
