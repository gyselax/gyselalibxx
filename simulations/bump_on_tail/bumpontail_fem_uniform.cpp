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
#include "bumpontailequilibrium.hpp"
#include "femnonperiodicpoissonsolver.hpp"
#include "femperiodicpoissonsolver.hpp"
#include "geometry.hpp"
#include "paraconfpp.hpp"
#include "params.yaml.hpp"
#include "pdi_out.yml.hpp"
#include "predcorr.hpp"
#include "singlemodeperturbinitialization.hpp"
#include "species_info.hpp"
#include "spline_interpolator_vx.hpp"
#include "spline_interpolator_x.hpp"
#include "splitvlasovsolver.hpp"

using std::cerr;
using std::endl;
using std::chrono::steady_clock;
namespace fs = std::filesystem;

int main(int argc, char** argv)
{
    ScopeGuard scope(argc, argv);

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
    PC_errhandler(PC_NULL_HANDLER);

    // Reading config
    // --> Mesh info
    CoordX const x_min(PCpp_double(conf_voicexx, ".Mesh.x_min"));
    CoordX const x_max(PCpp_double(conf_voicexx, ".Mesh.x_max"));
    IVectX const x_size(PCpp_int(conf_voicexx, ".Mesh.x_size"));
    CoordVx const vx_min(PCpp_double(conf_voicexx, ".Mesh.vx_min"));
    CoordVx const vx_max(PCpp_double(conf_voicexx, ".Mesh.vx_max"));
    IVectVx const vx_size(PCpp_int(conf_voicexx, ".Mesh.vx_size"));

    // Creating mesh & supports
    init_discrete_space<BSplinesX>(x_min, x_max, x_size);

    init_discrete_space<BSplinesVx>(vx_min, vx_max, vx_size);

    SplineXBuilder const builder_x;

    SplineVxBuilder const builder_vx;

    IVectSp const nb_kinspecies(PCpp_len(conf_voicexx, ".SpeciesInfo"));
    IDomainSp const dom_kinsp(IndexSp(0), nb_kinspecies);

    IDomainSpXVx const meshSpXVx(
            dom_kinsp,
            builder_x.interpolation_domain(),
            builder_vx.interpolation_domain());
    IDomainSpVx const meshSpVx(dom_kinsp, builder_vx.interpolation_domain());

    FieldSp<int> kinetic_charges(dom_kinsp);
    DFieldSp masses(dom_kinsp);
    DFieldSp epsilon_bot(dom_kinsp);
    DFieldSp temperature_bot(dom_kinsp);
    DFieldSp mean_velocity_bot(dom_kinsp);
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
        epsilon_bot(isp) = PCpp_double(conf_isp, ".epsilon_bot");
        temperature_bot(isp) = PCpp_double(conf_isp, ".temperature_bot");
        mean_velocity_bot(isp) = PCpp_double(conf_isp, ".mean_velocity_bot");
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
    init_discrete_space<IDimSp>(
            std::move(charges),
            std::move(masses),
            std::move(init_perturb_amplitude),
            std::move(init_perturb_mode));
    DFieldSpVx allfequilibrium(meshSpVx);
    BumpontailEquilibrium const init_fequilibrium(
            std::move(epsilon_bot),
            std::move(temperature_bot),
            std::move(mean_velocity_bot));
    init_fequilibrium(allfequilibrium);
    DFieldSpXVx allfdistribu(meshSpXVx);
    SingleModePerturbInitialization const
            init(allfequilibrium,
                 discrete_space<IDimSp>().perturb_modes(),
                 discrete_space<IDimSp>().perturb_amplitudes());
    init(allfdistribu);

    // --> Algorithm info
    double const deltat = PCpp_double(conf_voicexx, ".Algorithm.deltat");
    int const nbiter = static_cast<int>(PCpp_int(conf_voicexx, ".Algorithm.nbiter"));

    // --> Output info
    double const time_diag = PCpp_double(conf_voicexx, ".Output.time_diag");
    int const nbstep_diag = int(time_diag / deltat);

    PC_tree_t conf_pdi = PC_parse_string(PDI_CFG);

    PDI_init(conf_pdi);

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

    using FemPoissonSolverX = std::
            conditional_t<RDimX::PERIODIC, FemPeriodicPoissonSolver, FemNonPeriodicPoissonSolver>;
    FemPoissonSolverX const poisson(builder_x, spline_x_evaluator, builder_vx, spline_vx_evaluator);

    PredCorr const predcorr(vlasov, poisson, deltat);

    // Creating of mesh for output saving
    IDomainX const gridx = select<IDimX>(meshSpXVx);
    FieldX<CoordX> meshX_coord(gridx);
    for (IndexX const ix : gridx) {
        meshX_coord(ix) = coordinate(ix);
    }

    IDomainVx const gridvx = select<IDimVx>(meshSpXVx);
    FieldVx<CoordVx> meshVx_coord(gridvx);
    for (IndexVx const ivx : gridvx) {
        meshVx_coord(ivx) = coordinate(ivx);
    }

    // Starting the code
    expose_to_pdi("Nx", x_size.value());
    expose_to_pdi("Nvx", vx_size.value());
    expose_to_pdi("MeshX", meshX_coord);
    expose_to_pdi("MeshVx", meshVx_coord);
    expose_to_pdi("nbstep_diag", nbstep_diag);
    expose_to_pdi("Nkinspecies", nb_kinspecies.value());
    expose_to_pdi("fdistribu_charges", discrete_space<IDimSp>().charges()[dom_kinsp]);
    expose_to_pdi("fdistribu_masses", discrete_space<IDimSp>().masses()[dom_kinsp]);
    PdiEvent("initial_state").with("fdistribu_eq", allfequilibrium);

    steady_clock::time_point const start = steady_clock::now();

    predcorr(allfdistribu, nbiter);

    steady_clock::time_point const end = steady_clock::now();

    double const simulation_time = std::chrono::duration<double>(end - start).count();
    std::cout << "Simulation time: " << simulation_time << "s\n";

    PC_tree_destroy(&conf_pdi);

    PDI_finalize();

    PC_tree_destroy(&conf_voicexx);

    return EXIT_SUCCESS;
}
