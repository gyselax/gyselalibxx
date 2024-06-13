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
#include "fft_poisson_solver.hpp"
#include "geometry.hpp"
#include "input.hpp"
#include "maxwellianequilibrium.hpp"
#include "neumann_spline_quadrature.hpp"
#include "output.hpp"
#include "paraconfpp.hpp"
#include "params.yaml.hpp"
#include "pdi_out.yml.hpp"
#include "predcorr.hpp"
#include "qnsolver.hpp"
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
    IDomainX const mesh_x = init_spline_dependent_domain<
            IDimX,
            BSplinesX,
            SplineInterpPointsX>(conf_voicexx, "x");
    IDomainY const mesh_y = init_spline_dependent_domain<
            IDimY,
            BSplinesY,
            SplineInterpPointsY>(conf_voicexx, "y");
    IDomainVx const mesh_vx = init_spline_dependent_domain<
            IDimVx,
            BSplinesVx,
            SplineInterpPointsVx>(conf_voicexx, "vx");
    IDomainVy const mesh_vy = init_spline_dependent_domain<
            IDimVy,
            BSplinesVy,
            SplineInterpPointsVy>(conf_voicexx, "vy");
    IDomainXY const mesh_xy(mesh_x, mesh_y);
    IDomainVxVy mesh_vxvy(mesh_vx, mesh_vy);
    IDomainXYVxVy const meshXYVxVy(mesh_x, mesh_y, mesh_vx, mesh_vy);

    IVectSp const nb_kinspecies(PCpp_len(conf_voicexx, ".SpeciesInfo"));
    IDomainSp const dom_kinsp(IndexSp(0), nb_kinspecies);

    IDomainSpVxVy const meshSpVxVy(dom_kinsp, mesh_vx, mesh_vy);
    IDomainSpXYVxVy const meshSpXYVxVy(dom_kinsp, meshXYVxVy);

    SplineXBuilder const builder_x(meshXYVxVy);
    SplineYBuilder const builder_y(meshXYVxVy);
    SplineVxBuilder const builder_vx(meshXYVxVy);
    SplineVyBuilder const builder_vy(meshXYVxVy);
    SplineVxBuilder_1d const builder_vx_1d(mesh_vx);
    SplineVyBuilder_1d const builder_vy_1d(mesh_vy);

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

    ddc::ConstantExtrapolationRule<RDimVx> bv_vx_min(ddc::coordinate(mesh_vx.front()));
    ddc::ConstantExtrapolationRule<RDimVx> bv_vx_max(ddc::coordinate(mesh_vx.back()));
    SplineVxEvaluator const spline_vx_evaluator(bv_vx_min, bv_vx_max);

    PreallocatableSplineInterpolator const spline_vx_interpolator(builder_vx, spline_vx_evaluator);

    ddc::ConstantExtrapolationRule<RDimVy> bv_vy_min(ddc::coordinate(mesh_vy.front()));
    ddc::ConstantExtrapolationRule<RDimVy> bv_vy_max(ddc::coordinate(mesh_vy.back()));
    SplineVyEvaluator const spline_vy_evaluator(bv_vy_min, bv_vy_max);

    PreallocatableSplineInterpolator const spline_vy_interpolator(builder_vy, spline_vy_evaluator);

    // Create advection operator
    BslAdvectionSpatial<GeometryXYVxVy, IDimX> const advection_x(spline_x_interpolator);
    BslAdvectionSpatial<GeometryXYVxVy, IDimY> const advection_y(spline_y_interpolator);
    BslAdvectionVelocity<GeometryXYVxVy, IDimVx> const advection_vx(spline_vx_interpolator);
    BslAdvectionVelocity<GeometryXYVxVy, IDimVy> const advection_vy(spline_vy_interpolator);

    SplitVlasovSolver const vlasov(advection_x, advection_y, advection_vx, advection_vy);

    host_t<DFieldVxVy> const quadrature_coeffs_host
            = neumann_spline_quadrature_coefficients(mesh_vxvy, builder_vx_1d, builder_vy_1d);
    auto quadrature_coeffs = ddc::create_mirror_view_and_copy(
            Kokkos::DefaultExecutionSpace(),
            quadrature_coeffs_host.span_view());
    FFTPoissonSolver<IDomainXY, IDomainXY, Kokkos::DefaultExecutionSpace> fft_poisson_solver(
            mesh_xy);
    ChargeDensityCalculator const rhs(quadrature_coeffs);
    QNSolver const poisson(fft_poisson_solver, rhs);

    // Create predcorr operator
    PredCorr const predcorr(vlasov, poisson);

    // Starting the code
    ddc::expose_to_pdi("Nx_spline_cells", ddc::discrete_space<BSplinesX>().ncells());
    ddc::expose_to_pdi("Ny_spline_cells", ddc::discrete_space<BSplinesY>().ncells());
    ddc::expose_to_pdi("Nvx_spline_cells", ddc::discrete_space<BSplinesVx>().ncells());
    ddc::expose_to_pdi("Nvy_spline_cells", ddc::discrete_space<BSplinesVy>().ncells());
    expose_mesh_to_pdi("MeshX", mesh_x);
    expose_mesh_to_pdi("MeshY", mesh_y);
    expose_mesh_to_pdi("MeshVx", mesh_vx);
    expose_mesh_to_pdi("MeshVy", mesh_vy);
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
