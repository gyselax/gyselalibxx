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
#include <ddc/kernels/splines.hpp>

#include <collisions_inter.hpp>
#include <paraconf.h>
#include <pdi.h>

#include "bsl_advection_vx_batched.hpp"
#include "bsl_advection_x_batched.hpp"
#include "collisions_intra.hpp"
#ifdef PERIODIC_RDIMX
#include "femperiodicpoissonsolver.hpp"
#else
#include "femnonperiodicpoissonsolver.hpp"
#endif
#include "Lagrange_interpolator_batched.hpp"
#include "chargedensitycalculator.hpp"
#include "fftpoissonsolver.hpp"
#include "geometry.hpp"
#include "irighthandside.hpp"
#include "kinetic_source.hpp"
#include "krook_source_adaptive.hpp"
#include "krook_source_constant.hpp"
#include "maxwellianequilibrium.hpp"
#include "neumann_spline_quadrature.hpp"
#include "paraconfpp.hpp"
#include "pdi_out.yml.hpp"
#include "predcorr.hpp"
#include "restartinitialization.hpp"
#include "sheath.yaml.hpp"
#include "singlemodeperturbinitialization.hpp"
#include "species_info.hpp"
#include "spline_interpolator_batched.hpp"
#include "splitrighthandsidesolver.hpp"
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
    CoordX const x_min(PCpp_double(conf_voicexx, ".SplineMesh.x_min"));
    CoordX const x_max(PCpp_double(conf_voicexx, ".SplineMesh.x_max"));
    IVectX const x_ncells(PCpp_int(conf_voicexx, ".SplineMesh.x_ncells"));
    CoordVx const vx_min(PCpp_double(conf_voicexx, ".SplineMesh.vx_min"));
    CoordVx const vx_max(PCpp_double(conf_voicexx, ".SplineMesh.vx_max"));
    IVectVx const vx_ncells(PCpp_int(conf_voicexx, ".SplineMesh.vx_ncells"));

    // Creating mesh & supports
    ddc::init_discrete_space<BSplinesX>(x_min, x_max, x_ncells);

    ddc::init_discrete_space<BSplinesVx>(vx_min, vx_max, vx_ncells);

    ddc::init_discrete_space<IDimX>(SplineInterpPointsX::get_sampling());
    ddc::init_discrete_space<IDimVx>(SplineInterpPointsVx::get_sampling());
    IDomainX interpolation_domain_x(SplineInterpPointsX::get_domain());
    IDomainVx interpolation_domain_vx(SplineInterpPointsVx::get_domain());
    IDomainXVx meshXVx(interpolation_domain_x, interpolation_domain_vx);

    IVectSp const nb_kinspecies(PCpp_len(conf_voicexx, ".SpeciesInfo"));
    IDomainSp const dom_kinsp(IndexSp(0), nb_kinspecies);

    IDomainSpXVx const meshSpXVx(dom_kinsp, meshXVx);
    IDomainSpVx const meshSpVx(dom_kinsp, interpolation_domain_vx);

    SplineXBuilder const builder_x(meshXVx);
#ifndef PERIODIC_RDIMX
    SplineXBuilder_1d const builder_x_poisson(interpolation_domain_x);
#endif
    SplineVxBuilder const builder_vx(meshXVx);
    SplineVxBuilder_1d const builder_vx_poisson(interpolation_domain_vx);

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
        init_perturb_mode(isp) = static_cast<int>(PCpp_int(conf_isp, ".perturb_mode"));
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
    DFieldSpVx allfequilibrium(meshSpVx);
    MaxwellianEquilibrium const init_fequilibrium(
            std::move(density_eq),
            std::move(temperature_eq),
            std::move(mean_velocity_eq));
    init_fequilibrium(allfequilibrium);

    ddc::expose_to_pdi("iter_start", iter_start);

    DFieldSpXVx allfdistribu(meshSpXVx);
    double time_start(0);
    if (iter_start == 0) {
        SingleModePerturbInitialization const
                init(allfequilibrium,
                     init_perturb_mode.span_cview(),
                     init_perturb_amplitude.span_cview());
        init(allfdistribu);
    } else {
        RestartInitialization const restart(iter_start, time_start);
        restart(allfdistribu);
    }
    auto allfequilibrium_host = ddc::create_mirror_view_and_copy(allfequilibrium.span_view());

    // --> Algorithm info
    double const deltat = PCpp_double(conf_voicexx, ".Algorithm.deltat");
    int const nbiter = static_cast<int>(PCpp_int(conf_voicexx, ".Algorithm.nbiter"));

    // --> Output info
    double const time_diag = PCpp_double(conf_voicexx, ".Output.time_diag");
    int const nbstep_diag = int(time_diag / deltat);

#ifdef PERIODIC_RDIMX
    ddc::PeriodicExtrapolationRule<RDimX> bv_x_min;
    ddc::PeriodicExtrapolationRule<RDimX> bv_x_max;
#else
    ddc::ConstantExtrapolationRule<RDimX> bv_x_min(x_min);
    ddc::ConstantExtrapolationRule<RDimX> bv_x_max(x_max);
#endif

    // Creating batched operators
    SplineXEvaluator const spline_x_evaluator(bv_x_min, bv_x_max);
#ifndef PERIODIC_RDIMX
    SplineXEvaluator_1d const spline_x_evaluator_poisson(bv_x_min, bv_x_max);
#endif
    PreallocatableSplineInterpolatorBatched const
            spline_x_interpolator(builder_x, spline_x_evaluator);

    IVectVx static constexpr gwvx {0};
    LagrangeInterpolatorBatched<IDimVx, BCond::DIRICHLET, BCond::DIRICHLET, IDimX, IDimVx> const
            lagrange_vx_non_preallocatable_interpolator(3, gwvx);
    PreallocatableLagrangeInterpolatorBatched<
            IDimVx,
            BCond::DIRICHLET,
            BCond::DIRICHLET,
            IDimX,
            IDimVx> const lagrange_vx_interpolator(lagrange_vx_non_preallocatable_interpolator);

    BslAdvectionSpatialBatched<GeometryXVx, IDimX> const advection_x(spline_x_interpolator);
    BslAdvectionVelocityBatched<GeometryXVx, IDimVx> const advection_vx(lagrange_vx_interpolator);

    // Creating of mesh for output saving
    IDomainX const gridx = ddc::select<IDimX>(meshSpXVx);
    host_t<FieldX<CoordX>> meshX_coord(gridx);
    for (IndexX const ix : gridx) {
        meshX_coord(ix) = ddc::coordinate(ix);
    }

    IDomainVx const gridvx = ddc::select<IDimVx>(meshSpXVx);
    host_t<FieldVx<CoordVx>> meshVx_coord(gridvx);
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
            krook_source_adaptive_vector.emplace_back(
                    gridx,
                    gridvx,
                    type,
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

    host_t<DFieldVx> const quadrature_coeffs_host
            = neumann_spline_quadrature_coefficients(gridvx, builder_vx_poisson);

    auto const quadrature_coeffs = ddc::create_mirror_view_and_copy(
            Kokkos::DefaultExecutionSpace(),
            quadrature_coeffs_host.span_view());
    ChargeDensityCalculator rhs(quadrature_coeffs);
#ifdef PERIODIC_RDIMX
    ddc::init_discrete_space<IDimFx>(ddc::init_fourier_space<RDimX>(ddc::select<IDimX>(meshSpXVx)));
    FftPoissonSolver const poisson(rhs);
#else
    FemNonPeriodicPoissonSolver const poisson(builder_x_poisson, spline_x_evaluator_poisson, rhs);
#endif


    PredCorr const predcorr(boltzmann, poisson);

    // Starting the code
    ddc::expose_to_pdi("Nx_spline_cells", x_ncells.value());
    ddc::expose_to_pdi("Nvx_spline_cells", vx_ncells.value());
    ddc::expose_to_pdi("MeshX", meshX_coord);
    ddc::expose_to_pdi("MeshVx", meshVx_coord);
    ddc::expose_to_pdi("Lx", ddcHelper::total_interval_length(gridx));
    ddc::expose_to_pdi("nbstep_diag", nbstep_diag);
    ddc::expose_to_pdi("Nkinspecies", nb_kinspecies.value());
    ddc::expose_to_pdi("fdistribu_charges", ddc::discrete_space<IDimSp>().charges()[dom_kinsp]);
    ddc::expose_to_pdi("fdistribu_masses", ddc::discrete_space<IDimSp>().masses()[dom_kinsp]);
    ddc::PdiEvent("initial_state").with("fdistribu_eq", allfequilibrium_host);

    steady_clock::time_point const start = steady_clock::now();

    predcorr(allfdistribu, time_start, deltat, nbiter);

    steady_clock::time_point const end = steady_clock::now();

    double const simulation_time = std::chrono::duration<double>(end - start).count();
    std::cout << "Simulation time: " << simulation_time << "s\n";

    PC_tree_destroy(&conf_pdi);

    PDI_finalize();

    PC_tree_destroy(&conf_voicexx);

    return EXIT_SUCCESS;
}
