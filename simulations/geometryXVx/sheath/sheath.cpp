// SPDX-License-Identifier: MIT

#include <chrono>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#include <ddc/ddc.hpp>

#include <sll/constant_extrapolation_boundary_value.hpp>
#include <sll/spline_evaluator.hpp>

#include <collisions_inter.hpp>
#include <paraconf.h>
#include <pdi.h>

#include "bsl_advection_vx.hpp"
#include "bsl_advection_x.hpp"
#include "collisions_intra.hpp"
#ifdef PERIODIC_RDIMX
#include "femperiodicpoissonsolver.hpp"
#else
#include "femnonperiodicpoissonsolver.hpp"
#endif
#include "geometry.hpp"
#include "irighthandside.hpp"
#include "kinetic_source.hpp"
#include "krook_source_adaptive.hpp"
#include "krook_source_constant.hpp"
#include "maxwellianequilibrium.hpp"
#include "paraconfpp.hpp"
#include "pdi_out.yml.hpp"
#include "predcorr.hpp"
#include "sheath.yaml.hpp"
#include "singlemodeperturbinitialization.hpp"
#include "species_info.hpp"
#include "spline_interpolator.hpp"
#include "splitrighthandsidesolver.hpp"
#include "splitvlasovsolver.hpp"

using std::cerr;
using std::endl;
using std::chrono::steady_clock;
namespace fs = std::filesystem;

int main(int argc, char** argv)
{
    ddc::ScopeGuard scope(argc, argv);

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
    ddc::init_discrete_space<BSplinesX>(x_min, x_max, x_size);

    ddc::init_discrete_space<BSplinesVx>(vx_min, vx_max, vx_size);

    ddc::init_discrete_space<IDimX>(InterpPointsX::get_sampling());
    ddc::init_discrete_space<IDimVx>(InterpPointsVx::get_sampling());
    ddc::DiscreteDomain<IDimX> interpolation_domain_x(InterpPointsX::get_domain());
    ddc::DiscreteDomain<IDimVx> interpolation_domain_vx(InterpPointsVx::get_domain());

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
    DFieldSpVx allfequilibrium(meshSpVx);
    MaxwellianEquilibrium const init_fequilibrium(
            std::move(density_eq),
            std::move(temperature_eq),
            std::move(mean_velocity_eq));
    init_fequilibrium(allfequilibrium);
    DFieldSpXVx allfdistribu(meshSpXVx);
    SingleModePerturbInitialization const
            init(allfequilibrium,
                 ddc::discrete_space<IDimSp>().perturb_modes(),
                 ddc::discrete_space<IDimSp>().perturb_amplitudes());
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

    PreallocatableSplineInterpolator const spline_x_interpolator(builder_x, spline_x_evaluator);

    ConstantExtrapolationBoundaryValue<BSplinesVx> bv_v_min(vx_min);
    ConstantExtrapolationBoundaryValue<BSplinesVx> bv_v_max(vx_max);

    SplineEvaluator<BSplinesVx> const spline_vx_evaluator(bv_v_min, bv_v_max);

    PreallocatableSplineInterpolator const spline_vx_interpolator(builder_vx, spline_vx_evaluator);

    BslAdvectionSpatial<GeometryXVx, IDimX> const advection_x(spline_x_interpolator);

    BslAdvectionVelocity<GeometryXVx, IDimVx> const advection_vx(spline_vx_interpolator);

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

    // list of rhs operators
    std::vector<std::reference_wrapper<IRightHandSide const>> rhs_operators;
    std::vector<KrookSourceConstant> krook_source_constant_vector;
    std::vector<KrookSourceAdaptive> krook_source_adaptive_vector;
    // Krook operators initialization
    int const nb_rhsKrook(PCpp_len(conf_voicexx, ".Krook"));
    for (int ik = 0; ik < nb_rhsKrook; ++ik) {
        // --> Krook info
        PC_tree_t const conf_krook = PCpp_get(conf_voicexx, ".Krook[%d]", ik);

        static std::map<std::string, RhsType>
                str2rhstype {{"source", RhsType::Source}, {"sink", RhsType::Sink}};
        RhsType type = str2rhstype[PCpp_string(conf_krook, ".type")];
        std::string const krook_name = PCpp_string(conf_krook, ".name");
        if (krook_name == "constant") {
            krook_source_constant_vector.emplace_back(
                    gridx,
                    gridvx,
                    type,
                    PCpp_double(conf_krook, ".extent"),
                    PCpp_double(conf_krook, ".stiffness"),
                    PCpp_double(conf_krook, ".amplitude"),
                    PCpp_double(conf_krook, ".density"),
                    PCpp_double(conf_krook, ".temperature"));
            rhs_operators.emplace_back(krook_source_constant_vector.back());

        } else if (krook_name == "adaptive") {
            static std::map<std::string, RhsSolver> str2solver {{"rk2", RhsSolver::Rk2}};
            RhsSolver solver = str2solver[PCpp_string(conf_krook, ".solver")];

            krook_source_adaptive_vector.emplace_back(
                    gridx,
                    gridvx,
                    type,
                    solver,
                    PCpp_double(conf_krook, ".extent"),
                    PCpp_double(conf_krook, ".stiffness"),
                    PCpp_double(conf_krook, ".amplitude"),
                    PCpp_double(conf_krook, ".density"),
                    PCpp_double(conf_krook, ".temperature"));
            rhs_operators.emplace_back(krook_source_adaptive_vector.back());
        } else {
            throw std::invalid_argument(
                    "Invalid krook name, allowed values are: 'constant', or 'adaptive'.");
        }
    }

    // Kinetic source
    KineticSource const rhs_kinetic_source(
            gridx,
            gridvx,
            PCpp_double(conf_voicexx, ".KineticSource.extent"),
            PCpp_double(conf_voicexx, ".KineticSource.stiffness"),
            PCpp_double(conf_voicexx, ".KineticSource.amplitude"),
            PCpp_double(conf_voicexx, ".KineticSource.density"),
            PCpp_double(conf_voicexx, ".KineticSource.energy"),
            PCpp_double(conf_voicexx, ".KineticSource.temperature"));
    rhs_operators.emplace_back(rhs_kinetic_source);


    CollisionsIntra const
            collisions_intra(meshSpXVx, PCpp_double(conf_voicexx, ".CollisionsInfo.nustar0"));
    rhs_operators.emplace_back(collisions_intra);

    std::optional<CollisionsInter> collisions_inter;
    if (PCpp_bool(conf_voicexx, ".CollisionsInfo.enable_inter")) {
        collisions_inter.emplace(meshSpXVx, PCpp_double(conf_voicexx, ".CollisionsInfo.nustar0"));
        rhs_operators.emplace_back(*collisions_inter);
    }

    SplitVlasovSolver const vlasov(advection_x, advection_vx);
    SplitRightHandSideSolver const boltzmann(vlasov, rhs_operators);

#ifdef PERIODIC_RDIMX
    using FemPoissonSolverX = FemPeriodicPoissonSolver;
#else
    using FemPoissonSolverX = FemNonPeriodicPoissonSolver;
#endif
    FemPoissonSolverX const poisson(builder_x, spline_x_evaluator, builder_vx, spline_vx_evaluator);

    PredCorr const predcorr(boltzmann, poisson);

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

    predcorr(allfdistribu, deltat, nbiter);

    steady_clock::time_point const end = steady_clock::now();

    double const simulation_time = std::chrono::duration<double>(end - start).count();
    std::cout << "Simulation time: " << simulation_time << "s\n";

    PC_tree_destroy(&conf_pdi);

    PDI_finalize();

    PC_tree_destroy(&conf_voicexx);

    return EXIT_SUCCESS;
}
