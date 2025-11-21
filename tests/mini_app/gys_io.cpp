// SPDX-License-Identifier: MIT
#include <iostream>
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
//     PC_tree_t conf_gyselax;
    
//     // Parse YAML parameters file (first argument) or use empty config
//     if (argc > 1) {
//         conf_gyselax = PC_parse_path(argv[1]);
//     } else {
//         // Use default empty configuration if no file provided
//         conf_gyselax = PC_parse_string("");
//     }
//     PC_errhandler(PC_NULL_HANDLER);

//     // Parse PDI YAML file (second argument) or use default
//     PC_tree_t conf_pdi;
//     if (argc > 2) {
//         conf_pdi = PC_parse_path(argv[2]);
//     } else {
//         // Use default minimal PDI configuration
//         constexpr char const* const default_pdi_cfg = R"PDI_CFG(
// metadata:
//   gysela_io_filename_size: size_t
//   gysela_io_filename: {type: array, subtype: char, size: "$gysela_io_filename_size"}

// data:

// plugins:
//   decl_hdf5:
//     - file: '${gysela_io_filename}'
//       on_event: [read_species, read_tor1, read_tor2, read_vpar, read_mu, read_fdistribu]
//       read: [species, charges, masses, breakpoints_tor1, grid_tor1, breakpoints_tor2, grid_tor2, breakpoints_vpar, grid_vpar, breakpoints_mu, grid_mu, fdistribu_sptor2Dv2D]
// )PDI_CFG";
//         conf_pdi = PC_parse_string(default_pdi_cfg);
//     }

//     PDI_init(conf_pdi);

//     if (rank == 0) {
//         cout << "Initializing 5D particle distribution function..." << endl;
//     }

//     //----------------------------------------------
//     // Read species information from PDI
//     //----------------------------------------------
//     IdxRangeSp idxrange_glob_kinsp = init_kinetic_species();
//     if (rank == 0) {
//         cout << "Number of kinetic species: " << idxrange_glob_kinsp.size() << endl;
//     }

//     //----------------------------------------------
//     // Initialization of the IdxRange for all dimensions
//     //----------------------------------------------
//     // -- Read breakpoints and grid via PDI and deduce the IdxRange of each direction
//     IdxRangeTor1 const idxrange_glob_tor1 = init_spline_dependent_idx_range<
//             GridTor1,
//             BSplinesTor1,
//             SplineInterpPointsR>(conf_gyselax, "tor1");
//     IdxRangeTor2 const idxrange_glob_tor2 = init_spline_dependent_idx_range<
//             GridTor2,
//             BSplinesTor2,
//             SplineInterpPointsTor2>(conf_gyselax, "tor2");
//     IdxRangeTor3 const idxrange_glob_tor3 = init_spline_dependent_idx_range<
//             GridTor3,
//             BSplinesTor3,
//             SplineInterpPointsTor3>(conf_gyselax, "tor3");
//     IdxRangeVpar const idxrange_glob_vpar = init_spline_dependent_idx_range<
//             GridVpar,
//             BSplinesVpar,
//             SplineInterpPointsVpar>(conf_gyselax, "vpar");
//     IdxRangeMu const idxrange_glob_mu = init_spline_dependent_idx_range<
//             GridMu,
//             BSplinesMu,
//             SplineInterpPointsMu>(conf_gyselax, "mu");

//     if (rank == 0) {
//         cout << "Grid sizes:" << endl;
//         cout << "  tor1: " << idxrange_glob_tor1.size() << endl;
//         cout << "  tor2: " << idxrange_glob_tor2.size() << endl;
//         cout << "  tor3: " << idxrange_glob_tor3.size() << endl;
//         cout << "  vpar: " << idxrange_glob_vpar.size() << endl;
//         cout << "  mu: " << idxrange_glob_mu.size() << endl;
//     }

//     // -- Deduction of the useful global IdxRange for 5D distribution function (species is an index)
//     IdxRangeSpTor3DV2D const idxrange_glob_sptor3Dv2D(
//             idxrange_glob_kinsp,
//             idxrange_glob_tor1,
//             idxrange_glob_tor2,
//             idxrange_glob_tor3,
//             idxrange_glob_vpar,
//             idxrange_glob_mu);

//     if (rank == 0) {
//         cout << "5D distribution function size (with species index): " << idxrange_glob_sptor3Dv2D.size() << endl;
//     }

//     //----------------------------------------------
//     // Mesh saving (grid and breakpoints in all directions)
//     //----------------------------------------------
//     expose_mesh_to_pdi(
//             "breakpoints_tor1",
//             ddc::discrete_space<BSplinesTor1>().break_point_domain());
//     expose_mesh_to_pdi("grid_tor1", idxrange_glob_tor1);
//     expose_mesh_to_pdi(
//             "breakpoints_tor2",
//             ddc::discrete_space<BSplinesTor2>().break_point_domain());
//     expose_mesh_to_pdi("grid_tor2", idxrange_glob_tor2);
//     expose_mesh_to_pdi(
//             "breakpoints_vpar",
//             ddc::discrete_space<BSplinesVpar>().break_point_domain());
//     expose_mesh_to_pdi("grid_vpar", idxrange_glob_vpar);
//     expose_mesh_to_pdi(
//             "breakpoints_tor3",
//             ddc::discrete_space<BSplinesTor3>().break_point_domain());
//     expose_mesh_to_pdi("grid_tor3", idxrange_glob_tor3);
//     expose_mesh_to_pdi(
//             "breakpoints_vpar",
//             ddc::discrete_space<BSplinesVpar>().break_point_domain());
//     expose_mesh_to_pdi("grid_vpar", idxrange_glob_vpar);
//     expose_mesh_to_pdi(
//             "breakpoints_mu",
//             ddc::discrete_space<BSplinesMu>().break_point_domain());
//     expose_mesh_to_pdi("grid_mu", idxrange_glob_mu);

//     //----------------------------------------------
//     // Read the 5D distribution function (tor1, tor2, tor3, vpar, mu) for each species
//     //----------------------------------------------
//     DFieldMemSpTor3DV2D_host fdistribu_sptor3Dv2D_host(idxrange_glob_sptor3Dv2D);
    
//     // Expose index range for parallel I/O
//     PDI_expose_idx_range(idxrange_glob_sptor3Dv2D, "local_fdistribu");
    
//     // Read distribution function via PDI
//     ddc::PdiEvent("read_fdistribu")
//             .with("fdistribu_sptor3Dv2D", fdistribu_sptor3Dv2D_host);

//     if (rank == 0) {
//         cout << "Distribution function read successfully." << endl;
//         cout << "Total size: " << idxrange_glob_sptor3Dv2D.size() << " elements" << endl;
//     }

//     // Copy to device memory (optional, for further processing)
//     auto fdistribu_sptor3Dv2D_alloc = ddc::create_mirror_view_and_copy(
//             Kokkos::DefaultExecutionSpace(),
//             get_field(fdistribu_sptor3Dv2D_host));
//     DFieldSpTor3DV2D fdistribu_sptor3Dv2D(get_field(fdistribu_sptor3Dv2D_alloc));

//     if (rank == 0) {
//         cout << "5D distribution function initialized successfully." << endl;
//     }

    PDI_finalize();
    MPI_Finalize();

    // PC_tree_destroy(&conf_pdi);
    // PC_tree_destroy(&conf_gyselax);

    return EXIT_SUCCESS;
}

