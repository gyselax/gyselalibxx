// SPDX-License-Identifier: MIT
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <string>
#include <mpi.h>

#include <ddc/ddc.hpp>

#include <paraconf.h>
#include <pdi.h>

#include "ddc_alias_inline_functions.hpp"
#include "geometry.hpp"
#include "input.hpp"
#include "mesh_builder.hpp"
#include "output.hpp"
#include "paraconfpp.hpp"
#include "pdi_helper.hpp"
#include "pdi_default.yml.hpp"
#include "species_init.hpp"
// #include "maxwellianequilibrium.hpp"

using std::cout;
using std::endl;
using std::string;

namespace {

struct ConfigHandles
{
    PC_tree_t conf_gyselax;
    PC_tree_t conf_pdi;
};

ConfigHandles parse_config_files(int argc, char** argv)
{
    ConfigHandles configs{};
    if (argc > 1) {
        configs.conf_gyselax = PC_parse_path(argv[1]);
    } else {
        configs.conf_gyselax = PC_parse_string("");
    }
    PC_errhandler(PC_NULL_HANDLER);

    if (argc > 2) {
        configs.conf_pdi = PC_parse_path(argv[2]);
    } else {
        configs.conf_pdi = PC_parse_string(PDI_CFG);
    }
    return configs;
}

void print_banner(int rank)
{
    if (rank != 0) {
        return;
    }
    cout << "==========================================" << endl;
    cout << "           GYSELA MINI APP               " << endl;
    cout << "==========================================" << endl;
}

void yaml_params_to_log(int rank, PC_tree_t conf_gyselax)
{
    if (rank != 0) {
        return;
    }
    cout << "Configuration parameters:" << endl;
    string const gysela_io_filename(PCpp_string(conf_gyselax, ".FileNames.gysela_io_filename"));
    cout << "  gysela_io_filename: " << gysela_io_filename << endl;
    cout << "  Mesh parameters:" << endl;
    for (string const& dim : {"Tor1", "Tor2", "Tor3", "Vpar", "Mu"}) {
        for (string const& key : {"ncells", "min", "max"}) {
            string const key_key = ".SplineMesh." + dim + "_" + key;
            cout << dim  << "_" << key << ":" << PCpp_double(conf_gyselax, key_key.c_str()) << endl;
        }
    }
}

IdxRangeSp init_species_from_yaml(PC_tree_t conf_gyselax)
{
    int const nb_species = PCpp_len(conf_gyselax, ".SpeciesInfo");
    IdxRangeSp idx_range_sp(IdxSp(0), IdxStepSp(nb_species));

    host_t<DFieldMemSp> kinetic_charges(idx_range_sp);
    host_t<DFieldMemSp> kinetic_masses(idx_range_sp);
    for (int i = 0; i < nb_species; ++i) {
        std::string const base = ".SpeciesInfo[" + std::to_string(i) + "]";
        kinetic_charges(IdxSp(i)) = PCpp_double(conf_gyselax, (base + ".charge").c_str());
        kinetic_masses(IdxSp(i)) = PCpp_double(conf_gyselax, (base + ".mass").c_str());
    }

    ddc::init_discrete_space<Species>(std::move(kinetic_charges), std::move(kinetic_masses));
    return idx_range_sp;
}

} // namespace

int main(int argc, char** argv)
{
    Kokkos::ScopeGuard scope(argc, argv);
    ddc::ScopeGuard ddc_scope(argc, argv);
    MPI_Init(&argc, &argv);
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    ConfigHandles configs = parse_config_files(argc, argv);
    PDI_init(configs.conf_pdi);

    print_banner(rank);
    yaml_params_to_log(rank, configs.conf_gyselax);

    if (rank == 0) {
        cout << "Initializing 5D particle distribution function..." << endl;
    }

    IdxRangeSp const idx_range_sp = init_species_from_yaml(configs.conf_gyselax);
    if (rank == 0) {
        cout << "Number of kinetic species: " << idx_range_sp.size() << endl;
    }

    //---------------------------------------------------------
    // Initialisation of the IdxRange using the input yaml file
    //----------------------------------------------
    IdxRange<GridTor1> const idx_range_tor1 = init_spline_dependent_idx_range<GridTor1,BSplinesTor1,SplineInterpPointsTor1>(configs.conf_gyselax, "Tor1");
    IdxRange<GridTor2> const idx_range_tor2 = init_spline_dependent_idx_range<GridTor2,BSplinesTor2,SplineInterpPointsTor2>(configs.conf_gyselax, "Tor2");
    IdxRange<GridTor3> const idx_range_tor3 = init_spline_dependent_idx_range<GridTor3,BSplinesTor3,SplineInterpPointsTor3>(configs.conf_gyselax, "Tor3");
    IdxRange<GridVpar> const idx_range_vpar = init_spline_dependent_idx_range<GridVpar,BSplinesVpar,SplineInterpPointsVpar>(configs.conf_gyselax, "Vpar");
    IdxRange<GridMu> const idx_range_mu = init_spline_dependent_idx_range<GridMu,BSplinesMu,SplineInterpPointsMu>(configs.conf_gyselax, "Mu");
    
    IdxRangeSpGrid const meshGridSp(idx_range_sp, idx_range_tor1, idx_range_tor2, idx_range_tor3, idx_range_vpar, idx_range_mu);
    IdxRangeSpVparMu const meshGridSpVparMu(idx_range_sp, idx_range_vpar, idx_range_mu);
    //---------------------------------------------------------
    // Print the grid sizes
    //---------------------------------------------------------
    if (rank == 0) {
        cout << "Grid sizes:" << endl;
        cout << "  tor1: " << idx_range_tor1.size() << endl;
        cout << "  tor2: " << idx_range_tor2.size() << endl;
        cout << "  tor3: " << idx_range_tor3.size() << endl;
        cout << "  vpar: " << idx_range_vpar.size() << endl;
        cout << "  mu: " << idx_range_mu.size() << endl;
    }

    //---------------------------------------------------------
    // Initialisation of the distribution function
    //---------------------------------------------------------
    DFieldMemSpVparMu allfequilibrium(meshGridSpVparMu);
    // std::unique_ptr<IEquilibrium> const init_fequilibrium
    //         = equilibrium::init_from_input(idx_range_sp, configs.conf_gyselax);
    // (*init_fequilibrium)(get_field(allfequilibrium));

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
    //         idx_range_sp,
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

    PC_tree_destroy(&configs.conf_pdi);
    PC_tree_destroy(&configs.conf_gyselax);

    return EXIT_SUCCESS;
}

