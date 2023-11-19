// SPDX-License-Identifier: MIT

#include <chrono>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string_view>

#include <ddc/ddc.hpp>

#include <sll/constant_extrapolation_boundary_value.hpp>
#include <sll/spline_evaluator.hpp>

#include <paraconf.h>
#include <pdi.h>

#include "bsl_advection_vx.hpp"
#include "bsl_advection_x.hpp"
#ifdef PERIODIC_RDIMX
#include "femperiodicpoissonsolver.hpp"
#else
#include "femnonperiodicpoissonsolver.hpp"
#endif
#include "geometry.hpp"
#include "maxwellianequilibrium.hpp"
#include "paraconfpp.hpp"
#include "params.yaml.hpp"
#include "pdi_out.yml.hpp"
#include "predcorr.hpp"
#include "restartinitialization.hpp"
#include "singlemodeperturbinitialization.hpp"
#include "species_info.hpp"
#include "spline_interpolator.hpp"
#include "splitvlasovsolver.hpp"

using std::cerr;
using std::endl;
using std::chrono::steady_clock;
namespace fs = std::filesystem;

using PreallocatableSplineInterpolatorX
        = PreallocatableSplineInterpolator<IDimX, BSplinesX, SplineXBoundary, SplineXBoundary>;
using PreallocatableSplineInterpolatorVx = PreallocatableSplineInterpolator<
        IDimVx,
        BSplinesVx,
        BoundCond::HERMITE,
        BoundCond::HERMITE>;
using BslAdvectionX = BslAdvectionSpatial<GeometryXVx, IDimX>;
using BslAdvectionVx = BslAdvectionVelocity<GeometryXVx, IDimVx>;

int main(int argc, char** argv)
{
    ddc::ScopeGuard scope(argc, argv);

    long int iter_start(0);
    PC_tree_t conf_voicexx;
    if (argc == 2) {
        conf_voicexx = PC_parse_path(fs::path(argv[1]).c_str());
    } else if (argc == 3) {
        if (argv[1] == std::string_view("--dump-config")) {
            std::fstream file(argv[2], std::fstream::out);
            file << params_yaml;
            return EXIT_SUCCESS;
        }
    } else if (argc == 4) {
        if (argv[1] == std::string_view("--iter-restart")) {
            iter_start = std::strtol(argv[2], NULL, 10);
            conf_voicexx = PC_parse_path(fs::path(argv[3]).c_str());
        }
    } else {
        cerr << "usage: " << argv[0] << " [--dump-config] <config_file.yml>" << endl;
        cerr << "or to perform a restart" << argv[0] << " [--iter-restart] <iter> <config_file.yml>"
             << endl;
        return EXIT_FAILURE;
    }
    PC_tree_t conf_pdi = PC_parse_string(PDI_CFG);
    PC_errhandler(PC_NULL_HANDLER);
    PDI_init(conf_pdi);

    // Reading config
    // --> Mesh info
    CoordX const x_min(PCpp_double(conf_voicexx, ".Mesh.x_min"));
    CoordX const x_max(PCpp_double(conf_voicexx, ".Mesh.x_max"));
    IVectX const x_size(PCpp_int(conf_voicexx, ".Mesh.x_size"));
    CoordVx const vx_min(PCpp_double(conf_voicexx, ".Mesh.vx_min"));
    CoordVx const vx_max(PCpp_double(conf_voicexx, ".Mesh.vx_max"));
    IVectVx const vx_size(PCpp_int(conf_voicexx, ".Mesh.vx_size"));

    // Creating mesh & supports
    ddc::init_discrete_space<BSplinesX>(x_min, x_max, x_size);

    ddc::init_discrete_space<BSplinesVx>(vx_min, vx_max, vx_size);

    ddc::init_discrete_space<IDimX>(SplineInterpPointsX::get_sampling());
    ddc::init_discrete_space<IDimVx>(SplineInterpPointsVx::get_sampling());
    ddc::DiscreteDomain<IDimX> interpolation_domain_x(SplineInterpPointsX::get_domain());
    ddc::DiscreteDomain<IDimVx> interpolation_domain_vx(SplineInterpPointsVx::get_domain());

    SplineXBuilder const builder_x(interpolation_domain_x);

    SplineVxBuilder const builder_vx(interpolation_domain_vx);

    IVectSp const nb_kinspecies(PCpp_len(conf_voicexx, ".SpeciesInfo"));
    IDomainSp const dom_kinsp(IndexSp(0), nb_kinspecies);

    IDomainSpXVx const meshSpXVx(
            dom_kinsp,
            builder_x.interpolation_domain(),
            builder_vx.interpolation_domain());
    IDomainSpVx const meshSpVx(dom_kinsp, builder_vx.interpolation_domain());

    FieldSp<int> kinetic_charges(dom_kinsp);
    DFieldSp masses(dom_kinsp);
    DFieldSp density_eq(dom_kinsp);
    DFieldSp temperature_eq(dom_kinsp);
    DFieldSp mean_velocity_eq(dom_kinsp);
    DFieldSp init_perturb_amplitude(dom_kinsp);
    FieldSp<int> init_perturb_mode(dom_kinsp);
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
        init_perturb_mode(isp) = static_cast<int>(PCpp_int(conf_isp, ".perturb_mode"));
    }

    // Create the domain of all species including kinetic species + adiabatic species (if existing)
    IDomainSp const
            dom_allsp(IndexSp(0), nb_kinspecies + nb_elec_adiabspecies + nb_ion_adiabspecies);
    FieldSp<int> charges(dom_allsp);
    for (IndexSp isp : dom_kinsp) {
        charges(isp) = kinetic_charges(isp);
    }
    if (nb_elec_adiabspecies + nb_ion_adiabspecies > 0) {
        charges(dom_kinsp.back() + 1) = nb_ion_adiabspecies - nb_elec_adiabspecies;
    }

    // Initialization of the distribution function
    ddc::init_discrete_space<IDimSp>(
            std::move(charges),
            std::move(masses),
            std::move(init_perturb_amplitude),
            std::move(init_perturb_mode));
    device_t<DFieldSpVx> allfequilibrium_device(meshSpVx);
    MaxwellianEquilibrium const init_fequilibrium(
            std::move(density_eq),
            std::move(temperature_eq),
            std::move(mean_velocity_eq));
    init_fequilibrium(allfequilibrium_device);

    ddc::expose_to_pdi("iter_start", iter_start);

    device_t<DFieldSpXVx> allfdistribu_device(meshSpXVx);
    double time_start(0);
    if (iter_start == 0) {
        SingleModePerturbInitialization const
                init(allfequilibrium_device,
                     ddc::host_discrete_space<IDimSp>().perturb_modes(),
                     ddc::host_discrete_space<IDimSp>().perturb_amplitudes());
        init(allfdistribu_device);
    } else {
        RestartInitialization const restart(iter_start, time_start);
        restart(allfdistribu_device);
    }
    auto allfequilibrium = ddc::create_mirror_view_and_copy(allfequilibrium_device.span_view());
    auto allfdistribu = ddc::create_mirror_view_and_copy(allfdistribu_device.span_view());

    // --> Algorithm info
    double const deltat = PCpp_double(conf_voicexx, ".Algorithm.deltat");
    int const nbiter = static_cast<int>(PCpp_int(conf_voicexx, ".Algorithm.nbiter"));

    // --> Output info
    double const time_diag = PCpp_double(conf_voicexx, ".Output.time_diag");
    int const nbstep_diag = int(time_diag / deltat);

    ConstantExtrapolationBoundaryValue<BSplinesX> bv_x_min(x_min);
    ConstantExtrapolationBoundaryValue<BSplinesX> bv_x_max(x_max);

    // Creating operators
    SplineEvaluator<BSplinesX> const spline_x_evaluator(bv_x_min, bv_x_max);

    PreallocatableSplineInterpolatorX const spline_x_interpolator(builder_x, spline_x_evaluator);

    ConstantExtrapolationBoundaryValue<BSplinesVx> bv_v_min(vx_min);
    ConstantExtrapolationBoundaryValue<BSplinesVx> bv_v_max(vx_max);

    SplineEvaluator<BSplinesVx> const spline_vx_evaluator(bv_v_min, bv_v_max);

    PreallocatableSplineInterpolatorVx const
            spline_vx_interpolator(builder_vx, spline_vx_evaluator);

    BslAdvectionX const advection_x(spline_x_interpolator);

    BslAdvectionVx const advection_vx(spline_vx_interpolator);

    SplitVlasovSolver const vlasov(advection_x, advection_vx);

#ifdef PERIODIC_RDIMX
    using FemPoissonSolverX = FemPeriodicPoissonSolver;
#else
    using FemPoissonSolverX = FemNonPeriodicPoissonSolver;
#endif

    FemPoissonSolverX const poisson(builder_x, spline_x_evaluator, builder_vx, spline_vx_evaluator);

    PredCorr const predcorr(vlasov, poisson);

    // Creating of mesh for output saving
    IDomainX const gridx = ddc::select<IDimX>(meshSpXVx);
    FieldX<CoordX> meshX_coord(gridx);
    for (IndexX const ix : gridx) {
        meshX_coord(ix) = ddc::coordinate(ix);
    }

    IDomainVx const gridvx = ddc::select<IDimVx>(meshSpXVx);
    FieldVx<CoordVx> meshVx_coord(gridvx);
    for (IndexVx const ivx : gridvx) {
        meshVx_coord(ivx) = ddc::coordinate(ivx);
    }

    // Starting the code
    ddc::expose_to_pdi("Nx", x_size.value());
    ddc::expose_to_pdi("Nvx", vx_size.value());
    ddc::expose_to_pdi("MeshX", meshX_coord);
    ddc::expose_to_pdi("MeshVx", meshVx_coord);
    ddc::expose_to_pdi("nbstep_diag", nbstep_diag);
    ddc::expose_to_pdi("Nkinspecies", nb_kinspecies.value());
    ddc::expose_to_pdi("fdistribu_charges", ddc::discrete_space<IDimSp>().charges()[dom_kinsp]);
    ddc::expose_to_pdi("fdistribu_masses", ddc::discrete_space<IDimSp>().masses()[dom_kinsp]);
    ddc::PdiEvent("initial_state").with("fdistribu_eq", allfequilibrium);

    steady_clock::time_point const start = steady_clock::now();

    predcorr(allfdistribu_device, time_start, deltat, nbiter);

    steady_clock::time_point const end = steady_clock::now();

    double const simulation_time = std::chrono::duration<double>(end - start).count();
    std::cout << "Simulation time: " << simulation_time << "s\n";

    PC_tree_destroy(&conf_pdi);

    PDI_finalize();

    PC_tree_destroy(&conf_voicexx);

    return EXIT_SUCCESS;
}
