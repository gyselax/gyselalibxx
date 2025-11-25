// SPDX-License-Identifier: MIT
#include <cstdint>
#include <iostream>
#include <string>
#include <mpi.h>

#include <ddc/ddc.hpp>

#include <paraconf.h>
#include <pdi.h>

#include "ddc_alias_inline_functions.hpp"
#include "geometry.hpp"
#include "input.hpp"
#include "output.hpp"
#include "paraconfpp.hpp"
#include "pdi_helper.hpp"
#include "pdi_default.yml.hpp"
#include "species_init.hpp"

using std::cout;
using std::endl;

int main(int argc, char** argv)
{
    Kokkos::ScopeGuard scope(argc, argv);
    ddc::ScopeGuard ddc_scope(argc, argv);
    MPI_Init(&argc, &argv);
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    //----------------------------------------------
    // Parse configuration files
    //----------------------------------------------
    PC_tree_t conf_gyselax;
    // Parse YAML parameters file (first argument) or use empty config
    if (argc > 1) {
        conf_gyselax = PC_parse_path(argv[1]);
    } else {
        // Use default empty configuration if no file provided
        conf_gyselax = PC_parse_string("");
    }
    PC_errhandler(PC_NULL_HANDLER);

    // Parse PDI YAML file (second argument) or use default
    PC_tree_t conf_pdi;
    if (argc > 2) {
        conf_pdi = PC_parse_path(argv[2]);
    } else {
        // Use default PDI configuration from header file
        conf_pdi = PC_parse_string(PDI_CFG);
    }
        
    PDI_init(conf_pdi);
    
    std::string const gysela_io_filename(PCpp_string(conf_gyselax, ".FileNames.gysela_io_filename"));
    int64_t const gysela_io_filename_size = gysela_io_filename.size();
    

    if (rank == 0) {
        cout << "==========================================" << endl;
        cout << "           GYSELA MINI APP               " << endl;
        cout << "==========================================" << endl;
    }

    // Print the configuration
    if (rank == 0) {
        cout << "Configuration parameters:" << endl;
        cout << "  gysela_io_filename: " << gysela_io_filename << endl;
        cout << "  Mesh parameters:" << endl;
        cout << "  npts_tor1: " << PCpp_int(conf_gyselax, ".Mesh.npts_tor1") << endl;
        cout << "  npts_tor2: " << PCpp_int(conf_gyselax, ".Mesh.npts_tor2") << endl;
        cout << "  npts_tor3: " << PCpp_int(conf_gyselax, ".Mesh.npts_tor3") << endl;
        cout << "  npts_vpar: " << PCpp_int(conf_gyselax, ".Mesh.npts_vpar") << endl;
        cout << "  npts_mu: " << PCpp_int(conf_gyselax, ".Mesh.npts_mu") << endl;
        cout << "  min_tor1: " << PCpp_double(conf_gyselax, ".Mesh.min_tor1") << endl;
        cout << "  max_tor1: " << PCpp_double(conf_gyselax, ".Mesh.max_tor1") << endl;
        cout << "  min_tor2: " << PCpp_double(conf_gyselax, ".Mesh.min_tor2") << endl;
        cout << "  max_tor2: " << PCpp_double(conf_gyselax, ".Mesh.max_tor2") << endl;
        cout << "  min_tor3: " << PCpp_double(conf_gyselax, ".Mesh.min_tor3") << endl;
        cout << "  max_tor3: " << PCpp_double(conf_gyselax, ".Mesh.max_tor3") << endl;
        cout << "  min_vpar: " << PCpp_double(conf_gyselax, ".Mesh.min_vpar") << endl;
        cout << "  max_vpar: " << PCpp_double(conf_gyselax, ".Mesh.max_vpar") << endl;
        cout << "  min_mu: " << PCpp_double(conf_gyselax, ".Mesh.min_mu") << endl;
        cout << "  max_mu: " << PCpp_double(conf_gyselax, ".Mesh.max_mu") << endl;
    }

        
    if (rank == 0) {
            cout << "Initializing 5D particle distribution function..." << endl;
    }
        
    //----------------------------------------------
    // Read species information from PDI
    //----------------------------------------------
    int const nb_species = PCpp_len(conf_gyselax, ".SpeciesInfo");
    IdxRangeSp idxrange_glob_kinsp(IdxSp(0), IdxStepSp(nb_species));
    
    host_t<DFieldMemSp> kinetic_charges(idxrange_glob_kinsp);
    host_t<DFieldMemSp> kinetic_masses(idxrange_glob_kinsp);
    
    for (int i = 0; i < nb_species; ++i) {
        std::string const base = ".SpeciesInfo[" + std::to_string(i) + "]";
        kinetic_charges(IdxSp(i)) = PCpp_double(conf_gyselax, (base + ".charge").c_str());
        kinetic_masses(IdxSp(i)) = PCpp_double(conf_gyselax, (base + ".mass").c_str());
    }
    
    ddc::init_discrete_space<Species>(
            std::move(kinetic_charges),
            std::move(kinetic_masses));  
    if (rank == 0) {
        cout << "Number of kinetic species: " << idxrange_glob_kinsp.size() << endl;
    }
        
    // //----------------------------------------------
    // // Initialization of the IdxRange for all dimensions
    // //----------------------------------------------
    // // -- Read breakpoints and grid via PDI and deduce the IdxRange of each direction
    // IdxRangeTor1 const idxrange_glob_tor1 = init_spline_dependent_idx_range<
    //         GridTor1,
    //         BSplinesTor1,
    //         SplineInterpPointsR>(conf_gyselax, "tor1");
    // IdxRangeTor2 const idxrange_glob_tor2 = init_spline_dependent_idx_range<
    //         GridTor2,
    //         BSplinesTor2,
    //         SplineInterpPointsTor2>(conf_gyselax, "tor2");
    // IdxRangeTor3 const idxrange_glob_tor3 = init_spline_dependent_idx_range<
    //         GridTor3,
    //         BSplinesTor3,
    //         SplineInterpPointsTor3>(conf_gyselax, "tor3");
    // IdxRangeVpar const idxrange_glob_vpar = init_spline_dependent_idx_range<
    //         GridVpar,
    //         BSplinesVpar,
    //         SplineInterpPointsVpar>(conf_gyselax, "vpar");
    // IdxRangeMu const idxrange_glob_mu = init_spline_dependent_idx_range<
    //         GridMu,
    //         BSplinesMu,
    //         SplineInterpPointsMu>(conf_gyselax, "mu");

    // if (rank == 0) {
    //     cout << "Grid sizes:" << endl;
    //     cout << "  tor1: " << idxrange_glob_tor1.size() << endl;
    //     cout << "  tor2: " << idxrange_glob_tor2.size() << endl;
    //     cout << "  tor3: " << idxrange_glob_tor3.size() << endl;
    //     cout << "  vpar: " << idxrange_glob_vpar.size() << endl;
    //     cout << "  mu: " << idxrange_glob_mu.size() << endl;
    // }

    // // -- Deduction of the useful global IdxRange for 5D distribution function (species is an index)
    // IdxRangeSpTor3DV2D const idxrange_glob_sptor3Dv2D(
    //         idxrange_glob_kinsp,
    //         idxrange_glob_tor1,
    //         idxrange_glob_tor2,
    //         idxrange_glob_tor3,
    //         idxrange_glob_vpar,
    //         idxrange_glob_mu);

    // if (rank == 0) {
    //     cout << "5D distribution function size (with species index): " << idxrange_glob_sptor3Dv2D.size() << endl;
    // }

    // //----------------------------------------------
    // // Mesh saving (grid and breakpoints in all directions)
    // //----------------------------------------------
    // expose_mesh_to_pdi(
    //         "breakpoints_tor1",
    //         ddc::discrete_space<BSplinesTor1>().break_point_domain());
    // expose_mesh_to_pdi("grid_tor1", idxrange_glob_tor1);
    // expose_mesh_to_pdi(
    //         "breakpoints_tor2",
    //         ddc::discrete_space<BSplinesTor2>().break_point_domain());
    // expose_mesh_to_pdi("grid_tor2", idxrange_glob_tor2);
    // expose_mesh_to_pdi(
    //         "breakpoints_vpar",
    //         ddc::discrete_space<BSplinesVpar>().break_point_domain());
    // expose_mesh_to_pdi("grid_vpar", idxrange_glob_vpar);
    // expose_mesh_to_pdi(
    //         "breakpoints_tor3",
    //         ddc::discrete_space<BSplinesTor3>().break_point_domain());
    // expose_mesh_to_pdi("grid_tor3", idxrange_glob_tor3);
    // expose_mesh_to_pdi(
    //         "breakpoints_vpar",
    //         ddc::discrete_space<BSplinesVpar>().break_point_domain());
    // expose_mesh_to_pdi("grid_vpar", idxrange_glob_vpar);
    // expose_mesh_to_pdi(
    //         "breakpoints_mu",
    //         ddc::discrete_space<BSplinesMu>().break_point_domain());
    // expose_mesh_to_pdi("grid_mu", idxrange_glob_mu);

    // //----------------------------------------------
    // // Read the 5D distribution function (tor1, tor2, tor3, vpar, mu) for each species
    // //----------------------------------------------
    // DFieldMemSpTor3DV2D_host fdistribu_sptor3Dv2D_host(idxrange_glob_sptor3Dv2D);
    
    // // Expose index range for parallel I/O
    // PDI_expose_idx_range(idxrange_glob_sptor3Dv2D, "local_fdistribu");
    
    // // Read distribution function via PDI
    // ddc::PdiEvent("read_fdistribu")
    //         .with("fdistribu_sptor3Dv2D", fdistribu_sptor3Dv2D_host);

    // if (rank == 0) {
    //     cout << "Distribution function read successfully." << endl;
    //     cout << "Total size: " << idxrange_glob_sptor3Dv2D.size() << " elements" << endl;
    // }

    // // Copy to device memory (optional, for further processing)
    // auto fdistribu_sptor3Dv2D_alloc = ddc::create_mirror_view_and_copy(
    //         Kokkos::DefaultExecutionSpace(),
    //         get_field(fdistribu_sptor3Dv2D_host));
    // DFieldSpTor3DV2D fdistribu_sptor3Dv2D(get_field(fdistribu_sptor3Dv2D_alloc));

    // if (rank == 0) {
    //     cout << "5D distribution function initialized successfully." << endl;
    // }

    PDI_finalize();
    MPI_Finalize();

    PC_tree_destroy(&conf_pdi);
    PC_tree_destroy(&conf_gyselax);

    return EXIT_SUCCESS;
}

