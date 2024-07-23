// SPDX-License-Identifier: MIT

#include <chrono>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>

#include <ddc/ddc.hpp>

#include <paraconf.h>
#include <pdi.h>

#include "CollisionSpVparMu.hpp"
#include "geometry.hpp"
#include "input.hpp"
#include "maxwellianequilibrium.hpp"
#include "noperturbinitialization.hpp"
#include "output.hpp"
#include "paraconfpp.hpp"
#include "params.yaml.hpp"
#include "pdi_out.yml.hpp"
#include "simpson_quadrature.hpp"
#include "species_info.hpp"
#include "species_init.hpp"

using std::cerr;
using std::cout;
using std::endl;
using std::chrono::steady_clock;

int main(int argc, char** argv)
{
    // Environments variables for profiling
    setenv("KOKKOS_TOOLS_LIBS", KP_KERNEL_TIMER_PATH, false);
    setenv("KOKKOS_TOOLS_TIMER_JSON", "true", false);

    Kokkos::ScopeGuard kokkos_scope(argc, argv);
    ddc::ScopeGuard ddc_scope(argc, argv);

    PC_tree_t conf_collision = parse_executable_arguments(argc, argv, params_yaml);
    PC_tree_t conf_pdi = PC_parse_string(PDI_CFG);
    PC_errhandler(PC_NULL_HANDLER);
    PDI_init(conf_pdi);


    // --------- INITIALISATION ---------
    // ---> Reading of the mesh configuration from input YAML file
    // -----> Reading of mesh info
    IdxRangeVpar const idxrange_vpar = init_spline_dependent_domain<
            GridVpar,
            BSplinesVpar,
            SplineInterpPointsVpar>(conf_collision, "vpar");
    IdxRangeMu const idxrange_mu = init_spline_dependent_domain<
            GridMu,
            BSplinesMu,
            SplineInterpPointsMu>(conf_collision, "mu");
    // -----> Reading of species info
    IdxRangeSp const idxrange_kinsp = init_species(conf_collision);
    IdxRangeSpVparMu const idxrange_spvparmu(idxrange_kinsp, idxrange_vpar, idxrange_mu);

    // ---> Initialisation of the Maxwellian equilibrium distribution
    DFieldSpVparMu allfequilibrium(idxrange_spvparmu);
    MaxwellianEquilibrium const init_fequilibrium
            = MaxwellianEquilibrium::init_from_input(idxrange_kinsp, conf_collision);
    init_fequilibrium(allfequilibrium);

    // ---> Initialisation of the distribution function as a pertubed Maxwellian
    DFieldSpVparMu allfdistribu(idxrange_spvparmu);
    NoPerturbInitialization const init(allfequilibrium);
    init(allfdistribu);

    // ---> Expose unchanged data (related to mesh and species) to PDI
    ddc::expose_to_pdi("Nvpar_spline_cells", idxrange_vpar.size());
    ddc::expose_to_pdi("Nmu_spline_cells", idxrange_mu.size());
    expose_mesh_to_pdi("grid_vpar", idxrange_vpar);
    expose_mesh_to_pdi("grid_mu", idxrange_mu);
    ddc::expose_to_pdi("Nkinspecies", idxrange_kinsp.size());
    ddc::expose_to_pdi(
            "fdistribu_charges",
            ddc::discrete_space<IDimSp>().charges()[idxrange_kinsp]);
    ddc::expose_to_pdi("fdistribu_masses", ddc::discrete_space<IDimSp>().masses()[idxrange_kinsp]);


    // --------- OPERATOR INITIALISATION ---------
    // ---> Initialisation of the Collision operator
    // nustar fixed to 1. => normalisation of time is equal to 1./nustar
    double const nustar = 1.0;
    // rg fixed to 1. => because already included in nustar
    double const rg = 1.0;
    // safety factor fixed to 1. => because already included in nustar
    double const safety_factor = 1.0;
    // B_norm is fixed to 1. => normalisation at the localisation of collisions
    double const B_norm = 1.0;
    // TODO: Simplify the construction of coeff_intdmu and coeff_indmu as soon as the possibililty to define the quadrature coefficients directly on GPU is available
    host_t<DFieldVpar> const coeff_intdvpar_host
            = simpson_quadrature_coefficients_1d(allfdistribu.domain<GridVpar>());
    host_t<DFieldMu> const coeff_intdmu_host
            = simpson_quadrature_coefficients_1d(allfdistribu.domain<GridMu>());
    auto coeff_intdvpar = ddc::create_mirror_view_and_copy(
            Kokkos::DefaultExecutionSpace(),
            coeff_intdvpar_host.span_cview());
    auto coeff_intdmu = ddc::create_mirror_view_and_copy(
            Kokkos::DefaultExecutionSpace(),
            coeff_intdmu_host.span_cview());
    CollisionSpVparMu collision_operator(
            idxrange_spvparmu,
            coeff_intdmu.span_cview(),
            coeff_intdvpar.span_cview(),
            nustar,
            rg,
            safety_factor,
            B_norm);


    // --------- TIME ITERATION ---------
    // ---> Reading of algorithm info from input YAML file
    double const deltat = PCpp_double(conf_collision, ".Algorithm.deltat");
    if (deltat > 1.e-2) {
        throw std::runtime_error(
                "Deltat must be inferior to 1.e-2 for numerical stability of the scheme.");
    }
    int const nbiter = static_cast<int>(PCpp_int(conf_collision, ".Algorithm.nbiter"));

    // ---> Reading of output info from input YAML file
    double const time_diag = PCpp_double(conf_collision, ".Output.time_diag");
    int const nbstep_diag = int(time_diag / deltat);
    ddc::expose_to_pdi("nbstep_diag", nbstep_diag);

    // ---> Algorithm time loop
    long int const iter_start(0);
    double const time_start(0.0);
    ddc::expose_to_pdi("iter_start", iter_start);

    steady_clock::time_point const start = steady_clock::now();

    int iter = 0;
    for (; iter < nbiter; ++iter) {
        double const time_iter = time_start + iter * deltat;
        cout << "iter = " << iter << " ; time_iter = " << time_iter << endl;

        // Apply collision operator
        collision_operator(allfdistribu, deltat);
        auto allfdistribu_host = ddc::create_mirror_view_and_copy(allfdistribu.span_view());

        // Write distribution function
        ddc::PdiEvent("write_fdistribu")
                .with("iter", iter)
                .and_with("time_saved", time_iter)
                .and_with("fdistribu", allfdistribu_host);
    }

    steady_clock::time_point const end = steady_clock::now();
    double const simulation_time = std::chrono::duration<double>(end - start).count();
    std::cout << "Simulation time: " << simulation_time << "s\n";


    // --------- FINALISATION ---------
    PDI_finalize();
    PC_tree_destroy(&conf_pdi);
    PC_tree_destroy(&conf_collision);

    return EXIT_SUCCESS;
}
