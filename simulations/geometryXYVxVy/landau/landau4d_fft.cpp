// SPDX-License-Identifier: MIT

#include <chrono>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string_view>

#include <ddc/ddc.hpp>

#include <paraconf.h>
#include <pdi.h>

#include "bsl_advection_vx.hpp"
#include "bsl_advection_x.hpp"
#include "fftpoissonsolver.hpp"
#include "geometry.hpp"
#include "maxwellianequilibrium.hpp"
#include "neumann_spline_quadrature.hpp"
#include "paraconfpp.hpp"
#include "params.yaml.hpp"
#include "pdi_out.yml.hpp"
#include "predcorr.hpp"
#include "singlemodeperturbinitialization.hpp"
//#include "species_info.hpp"
#include "spline_interpolator.hpp"
#include "splitvlasovsolver.hpp"

using std::cerr;
using std::endl;
using std::chrono::steady_clock;
namespace fs = std::filesystem;

int main(int argc, char** argv)
{
    // Environments variables for profiling
    setenv("KOKKOS_TOOLS_LIBS", KP_KERNEL_TIMER_PATH, false);
    setenv("KOKKOS_TOOLS_TIMER_JSON", "true", false);

    Kokkos::ScopeGuard kokkos_scope(argc, argv);
    ddc::ScopeGuard ddc_scope(argc, argv);

    PC_tree_t conf_voicexx;
    if (argc == 2) {
        conf_voicexx = PC_parse_path(fs::path(argv[1]).c_str());
    } else if (argc == 3) {
        if (argv[1] == std::string_view("--dump-config")) {
            std::fstream file(argv[2], std::fstream::out);
            file << params_yaml;
            return EXIT_SUCCESS;
        }
    } else {
        cerr << "usage: " << argv[0] << " [--dump-config] <config_file.yml>" << endl;
        return EXIT_FAILURE;
    }
    PC_tree_t conf_pdi = PC_parse_string(PDI_CFG);
    PC_errhandler(PC_NULL_HANDLER);
    PDI_init(conf_pdi);

    // Reading config
    // --> Mesh info
    CoordX const x_min(PCpp_double(conf_voicexx, ".SplineMesh.x_min"));
    CoordX const x_max(PCpp_double(conf_voicexx, ".SplineMesh.x_max"));
    IVectX const x_ncells(PCpp_int(conf_voicexx, ".SplineMesh.x_ncells"));

    CoordY const y_min(PCpp_double(conf_voicexx, ".SplineMesh.y_min"));
    CoordY const y_max(PCpp_double(conf_voicexx, ".SplineMesh.y_max"));
    IVectY const y_ncells(PCpp_int(conf_voicexx, ".SplineMesh.y_ncells"));

    CoordVx const vx_min(PCpp_double(conf_voicexx, ".SplineMesh.vx_min"));
    CoordVx const vx_max(PCpp_double(conf_voicexx, ".SplineMesh.vx_max"));
    IVectVx const vx_ncells(PCpp_int(conf_voicexx, ".SplineMesh.vx_ncells"));

    CoordVy const vy_min(PCpp_double(conf_voicexx, ".SplineMesh.vy_min"));
    CoordVy const vy_max(PCpp_double(conf_voicexx, ".SplineMesh.vy_max"));
    IVectVy const vy_ncells(PCpp_int(conf_voicexx, ".SplineMesh.vy_ncells"));

    // Creating mesh & supports
    ddc::init_discrete_space<BSplinesX>(x_min, x_max, x_ncells);
    ddc::init_discrete_space<BSplinesY>(y_min, y_max, y_ncells);
    ddc::init_discrete_space<BSplinesVx>(vx_min, vx_max, vx_ncells);
    ddc::init_discrete_space<BSplinesVy>(vy_min, vy_max, vy_ncells);
    ddc::init_discrete_space<IDimX>(SplineInterpPointsX::get_sampling<IDimX>());
    ddc::init_discrete_space<IDimY>(SplineInterpPointsY::get_sampling<IDimY>());
    ddc::init_discrete_space<IDimVx>(SplineInterpPointsVx::get_sampling<IDimVx>());
    ddc::init_discrete_space<IDimVy>(SplineInterpPointsVy::get_sampling<IDimVy>());

    IVectSp const nb_kinspecies(PCpp_len(conf_voicexx, ".SpeciesInfo"));
    IDomainSp const dom_kinsp(IndexSp(0), nb_kinspecies);

    IDomainX interpolation_domain_x(SplineInterpPointsX::get_domain<IDimX>());
    IDomainY interpolation_domain_y(SplineInterpPointsY::get_domain<IDimY>());
    IDomainVx interpolation_domain_vx(SplineInterpPointsVx::get_domain<IDimVx>());
    IDomainVy interpolation_domain_vy(SplineInterpPointsVy::get_domain<IDimVy>());
    IDomainVxVy interpolation_domain_vxvy(interpolation_domain_vx, interpolation_domain_vy);

    IDomainXYVxVy meshXYVxVy(
            interpolation_domain_x,
            interpolation_domain_y,
            interpolation_domain_vx,
            interpolation_domain_vy);
    IDomainSpVxVy const meshSpVxVy(dom_kinsp, interpolation_domain_vx, interpolation_domain_vy);
    IDomainSpXYVxVy const meshSpXYVxVy(dom_kinsp, meshXYVxVy);

    SplineXBuilder const builder_x(meshXYVxVy);
    SplineYBuilder const builder_y(meshXYVxVy);
    SplineVxBuilder const builder_vx(meshXYVxVy);
    SplineVyBuilder const builder_vy(meshXYVxVy);
    SplineVxBuilder_1d const builder_vx_1d(interpolation_domain_vx);
    SplineVyBuilder_1d const builder_vy_1d(interpolation_domain_vy);

    host_t<FieldSp<int>> kinetic_charges(dom_kinsp);
    host_t<DFieldSp> masses(dom_kinsp);
    host_t<DFieldSp> density_eq(dom_kinsp);
    host_t<DFieldSp> temperature_eq(dom_kinsp);
    host_t<DFieldSp> mean_velocity_eq(dom_kinsp);
    host_t<DFieldSp> init_perturb_amplitude(dom_kinsp);
    host_t<FieldSp<int>> init_perturb_mode(dom_kinsp);

    int nb_elec_adiabspecies = 1;
    int nb_ion_adiabspecies = 1;

    for (IndexSp const isp : dom_kinsp) {
        // --> SpeciesInfo info
        PC_tree_t const conf_isp = PCpp_get(conf_voicexx, ".SpeciesInfo[%d]", isp.uid());

        kinetic_charges(isp) = static_cast<int>(PCpp_int(conf_isp, ".charge"));
        if (kinetic_charges(isp) == -1) {
            nb_elec_adiabspecies = 0;
        } else {
            nb_ion_adiabspecies = 0;
        }

        masses(isp) = PCpp_double(conf_isp, ".mass");
        density_eq(isp) = PCpp_double(conf_isp, ".density_eq");
        temperature_eq(isp) = PCpp_double(conf_isp, ".temperature_eq");
        mean_velocity_eq(isp) = PCpp_double(conf_isp, ".mean_velocity_eq");
        init_perturb_amplitude(isp) = PCpp_double(conf_isp, ".perturb_amplitude");
        init_perturb_mode(isp) = PCpp_double(conf_isp, ".perturb_mode");
    }

    // Create the domain of all species including kinetic species + adiabatic species (if existing)
    IDomainSp const
            dom_allsp(IndexSp(0), nb_kinspecies + nb_elec_adiabspecies + nb_ion_adiabspecies);
    host_t<FieldSp<int>> charges(dom_allsp);
    for (IndexSp isp : dom_kinsp) {
        charges(isp) = kinetic_charges(isp);
    }
    if (nb_elec_adiabspecies + nb_ion_adiabspecies > 0) {
        charges(dom_kinsp.back() + 1) = nb_ion_adiabspecies - nb_elec_adiabspecies;
    }

    // Initialization of the distribution function
    ddc::init_discrete_space<IDimSp>(std::move(charges), std::move(masses));
    DFieldSpVxVy allfequilibrium(meshSpVxVy);
    MaxwellianEquilibrium const init_fequilibrium(
            std::move(density_eq),
            std::move(temperature_eq),
            std::move(mean_velocity_eq));
    init_fequilibrium(allfequilibrium);
    DFieldSpXYVxVy allfdistribu(meshSpXYVxVy);
    SingleModePerturbInitialization const
            init(allfequilibrium,
                 init_perturb_mode.span_cview(),
                 init_perturb_amplitude.span_cview());
    init(allfdistribu);
    auto allfequilibrium_host = ddc::create_mirror_view_and_copy(allfequilibrium.span_view());

    // --> Algorithm info
    double const deltat = PCpp_double(conf_voicexx, ".Algorithm.deltat");
    int const nbiter = static_cast<int>(PCpp_int(conf_voicexx, ".Algorithm.nbiter"));

    // --> Output info
    double const time_diag = PCpp_double(conf_voicexx, ".Output.time_diag");
    int const nbstep_diag = int(time_diag / deltat);

    // Create spline evaluator
    ddc::PeriodicExtrapolationRule<RDimX> bv_x_min;
    ddc::PeriodicExtrapolationRule<RDimX> bv_x_max;
    SplineXEvaluator const spline_x_evaluator(bv_x_min, bv_x_max);

    PreallocatableSplineInterpolator const spline_x_interpolator(builder_x, spline_x_evaluator);

    ddc::PeriodicExtrapolationRule<RDimY> bv_y_min;
    ddc::PeriodicExtrapolationRule<RDimY> bv_y_max;
    SplineYEvaluator const spline_y_evaluator(bv_y_min, bv_y_max);

    PreallocatableSplineInterpolator const spline_y_interpolator(builder_y, spline_y_evaluator);

    ddc::ConstantExtrapolationRule<RDimVx> bv_vx_min(vx_min);
    ddc::ConstantExtrapolationRule<RDimVx> bv_vx_max(vx_max);
    SplineVxEvaluator const spline_vx_evaluator(bv_vx_min, bv_vx_max);

    PreallocatableSplineInterpolator const spline_vx_interpolator(builder_vx, spline_vx_evaluator);

    ddc::ConstantExtrapolationRule<RDimVy> bv_vy_min(vy_min);
    ddc::ConstantExtrapolationRule<RDimVy> bv_vy_max(vy_max);
    SplineVyEvaluator const spline_vy_evaluator(bv_vy_min, bv_vy_max);

    PreallocatableSplineInterpolator const spline_vy_interpolator(builder_vy, spline_vy_evaluator);

    // Create advection operator
    BslAdvectionSpatial<GeometryXYVxVy, IDimX> const advection_x(spline_x_interpolator);
    BslAdvectionSpatial<GeometryXYVxVy, IDimY> const advection_y(spline_y_interpolator);
    BslAdvectionVelocity<GeometryXYVxVy, IDimVx> const advection_vx(spline_vx_interpolator);
    BslAdvectionVelocity<GeometryXYVxVy, IDimVy> const advection_vy(spline_vy_interpolator);

    SplitVlasovSolver const vlasov(advection_x, advection_y, advection_vx, advection_vy);

    ddc::init_discrete_space<IDimFx>(
            ddc::init_fourier_space<IDimFx>(ddc::select<IDimX>(meshSpXYVxVy)));
    ddc::init_discrete_space<IDimFy>(
            ddc::init_fourier_space<IDimFy>(ddc::select<IDimY>(meshSpXYVxVy)));

    host_t<DFieldVxVy> const quadrature_coeffs_host = neumann_spline_quadrature_coefficients(
            interpolation_domain_vxvy,
            builder_vx_1d,
            builder_vy_1d);
    auto quadrature_coeffs = ddc::create_mirror_view_and_copy(
            Kokkos::DefaultExecutionSpace(),
            quadrature_coeffs_host.span_view());
    ChargeDensityCalculator const rhs(quadrature_coeffs);
    FftPoissonSolver const poisson(rhs);

    // Create predcorr operator
    PredCorr const predcorr(vlasov, poisson);

    // Creating of mesh for output saving
    host_t<FieldX<CoordX>> meshX_coord(interpolation_domain_x);
    ddc::for_each(interpolation_domain_x, [&](IndexX const ix) {
        meshX_coord(ix) = ddc::coordinate(ix);
    });

    host_t<FieldY<CoordY>> meshY_coord(interpolation_domain_y);
    ddc::for_each(interpolation_domain_y, [&](IndexY const iy) {
        meshY_coord(iy) = ddc::coordinate(iy);
    });

    host_t<FieldVx<CoordVx>> meshVx_coord(interpolation_domain_vx);
    for (IndexVx const ivx : interpolation_domain_vx) {
        meshVx_coord(ivx) = ddc::coordinate(ivx);
    }

    host_t<FieldVy<CoordVy>> meshVy_coord(interpolation_domain_vy);
    for (IndexVy const ivy : interpolation_domain_vy) {
        meshVy_coord(ivy) = ddc::coordinate(ivy);
    }

    // Starting the code
    ddc::expose_to_pdi("Nx_spline_cells", x_ncells.value());
    ddc::expose_to_pdi("Ny_spline_cells", y_ncells.value());
    ddc::expose_to_pdi("Nvx_spline_cells", vx_ncells.value());
    ddc::expose_to_pdi("Nvy_spline_cells", vy_ncells.value());
    ddc::expose_to_pdi("MeshX", meshX_coord);
    ddc::expose_to_pdi("MeshY", meshY_coord);
    ddc::expose_to_pdi("MeshVx", meshVx_coord);
    ddc::expose_to_pdi("MeshVy", meshVy_coord);
    ddc::expose_to_pdi("nbstep_diag", nbstep_diag);
    ddc::expose_to_pdi("Nkinspecies", nb_kinspecies.value());
    ddc::expose_to_pdi("fdistribu_charges", ddc::discrete_space<IDimSp>().charges()[dom_kinsp]);
    ddc::expose_to_pdi("fdistribu_masses", ddc::discrete_space<IDimSp>().masses()[dom_kinsp]);
    ddc::PdiEvent("initial_state").with("fdistribu_eq", allfequilibrium_host);

    steady_clock::time_point const start = steady_clock::now();

    predcorr(allfdistribu.span_view(), deltat, nbiter);

    steady_clock::time_point const end = steady_clock::now();

    double const simulation_time = std::chrono::duration<double>(end - start).count();
    std::cout << "Simulation time: " << simulation_time << "s\n";

    PC_tree_destroy(&conf_pdi);

    PDI_finalize();

    PC_tree_destroy(&conf_voicexx);

    return EXIT_SUCCESS;
}
