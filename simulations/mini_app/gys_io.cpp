// SPDX-License-Identifier: MIT
#include <algorithm>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
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
#include "mpitransposealltoall.hpp"
#include "transpose.hpp"
// #include "maxwellianequilibrium.hpp"

using std::cout;
using std::endl;
using std::string;
using std::chrono::steady_clock;

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

struct MeshInitializationResult
{
    IdxRangeSp idx_range_sp;
    IdxRangeSpGrid mesh_sp;
    IdxRangeSpVparMu mesh_sp_vparmu;
    IdxRange<GridTor1> idx_range_tor1;
    IdxRange<GridTor2> idx_range_tor2;
    IdxRange<GridTor3> idx_range_tor3;
    IdxRange<GridVpar> idx_range_vpar;
    IdxRange<GridMu> idx_range_mu;
};

MeshInitializationResult initialize_mesh(int rank, PC_tree_t conf_gyselax)
{
    IdxRangeSp const idx_range_sp = init_species_from_yaml(conf_gyselax);
    if (rank == 0) {
        cout << "Number of kinetic species: " << idx_range_sp.size() << endl;
    }

    IdxRange<GridTor1> const idx_range_tor1
            = init_spline_dependent_idx_range<GridTor1, BSplinesTor1, SplineInterpPointsTor1>(
                    conf_gyselax,
                    "Tor1");
    IdxRange<GridTor2> const idx_range_tor2
            = init_spline_dependent_idx_range<GridTor2, BSplinesTor2, SplineInterpPointsTor2>(
                    conf_gyselax,
                    "Tor2");
    IdxRange<GridTor3> const idx_range_tor3
            = init_spline_dependent_idx_range<GridTor3, BSplinesTor3, SplineInterpPointsTor3>(
                    conf_gyselax,
                    "Tor3");
    IdxRange<GridVpar> const idx_range_vpar
            = init_spline_dependent_idx_range<GridVpar, BSplinesVpar, SplineInterpPointsVpar>(
                    conf_gyselax,
                    "Vpar");
    IdxRange<GridMu> const idx_range_mu
            = init_spline_dependent_idx_range<GridMu, BSplinesMu, SplineInterpPointsMu>(
                    conf_gyselax,
                    "Mu");

    if (rank == 0) {
        cout << "Grid sizes:" << endl;
        cout << "  tor1: " << idx_range_tor1.size() << endl;
        cout << "  tor2: " << idx_range_tor2.size() << endl;
        cout << "  tor3: " << idx_range_tor3.size() << endl;
        cout << "  vpar: " << idx_range_vpar.size() << endl;
        cout << "  mu: " << idx_range_mu.size() << endl;
    }

    MeshInitializationResult result{
            idx_range_sp,
            IdxRangeSpGrid(idx_range_sp, idx_range_tor1, idx_range_tor2, idx_range_tor3, idx_range_vpar, idx_range_mu),
            IdxRangeSpVparMu(idx_range_sp, idx_range_vpar, idx_range_mu),
            idx_range_tor1,
            idx_range_tor2,
            idx_range_tor3,
            idx_range_vpar,
            idx_range_mu};
    return result;
}

void init_distribution_fun(
        DFieldMemSpGrid& allfdistribu,
        IdxRangeSpVparMu const& meshGridSpVparMu,
        IdxRangeSpGrid const& meshGridSp)
{
    double const mean_velocity = 1.0;
    double const temperature = 1.0;
    DFieldMemSpVparMu allfequilibrium(meshGridSpVparMu);
    auto allfequilibrium_field = get_field(allfequilibrium);
    auto allfdistribu_field = get_field(allfdistribu);

    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            meshGridSpVparMu,
            KOKKOS_LAMBDA(IdxSpVparMu const ispvparmu) {
                CoordVpar const vpar = ddc::coordinate(ddc::select<GridVpar>(ispvparmu));
                double const vpar_value = static_cast<double>(vpar);
                allfequilibrium_field(ispvparmu)
                        = Kokkos::exp(-(vpar_value - mean_velocity) * (vpar_value - mean_velocity)
                                      / (2. * temperature));
            });

    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            meshGridSp,
            KOKKOS_LAMBDA(IdxSpGrid const ispgrid) {
                IdxSpVparMu const ispvparmu(
                        ddc::select<Species>(ispgrid),
                        ddc::select<GridVpar>(ispgrid),
                        ddc::select<GridMu>(ispgrid));
                allfdistribu_field(ispgrid) = allfequilibrium_field(ispvparmu);
            });
}

void write_fdistribu(
        int rank,
        MeshInitializationResult const& mesh,
        host_t<DFieldMemSpGrid> const& allfdistribu_host)
{
    if (rank == 0) {
        cout << "Writing 5D distribution function and coordinates to file." << endl;
    }

    // Expose index range for parallel I/O
    PDI_expose_idx_range(mesh.mesh_sp, "local_fdistribu");

    // Expose species extents
    std::size_t species_extent = mesh.idx_range_sp.size();
    std::array<std::size_t, 1> species_extents_arr = {species_extent};
    PDI_expose("species_extents", species_extents_arr.data(), PDI_OUT);

    // Expose coordinate extents
    std::array<std::size_t, 1> tor1_extents_arr = {mesh.idx_range_tor1.size()};
    std::array<std::size_t, 1> tor2_extents_arr = {mesh.idx_range_tor2.size()};
    std::array<std::size_t, 1> tor3_extents_arr = {mesh.idx_range_tor3.size()};
    std::array<std::size_t, 1> vpar_extents_arr = {mesh.idx_range_vpar.size()};
    std::array<std::size_t, 1> mu_extents_arr = {mesh.idx_range_mu.size()};
    PDI_expose("tor1_extents", tor1_extents_arr.data(), PDI_OUT);
    PDI_expose("tor2_extents", tor2_extents_arr.data(), PDI_OUT);
    PDI_expose("tor3_extents", tor3_extents_arr.data(), PDI_OUT);
    PDI_expose("vpar_extents", vpar_extents_arr.data(), PDI_OUT);
    PDI_expose("mu_extents", mu_extents_arr.data(), PDI_OUT);

    // Expose coordinates to PDI
    expose_mesh_to_pdi("tor1", mesh.idx_range_tor1);
    expose_mesh_to_pdi("tor2", mesh.idx_range_tor2);
    expose_mesh_to_pdi("tor3", mesh.idx_range_tor3);
    expose_mesh_to_pdi("vpar", mesh.idx_range_vpar);
    expose_mesh_to_pdi("mu", mesh.idx_range_mu);

    // Expose distribution function to PDI and trigger write event
    ddc::PdiEvent("write_fdistribu")
            .with("fdistribu_sptor3Dv2D", allfdistribu_host);

    if (rank == 0) {
        cout << "5D distribution function and coordinates written successfully." << endl;
    }
}

void write_cpu_time_stats(
        int rank,
        double const* durations,
        std::vector<std::string> const& names,
        std::size_t num_entries)
{
    if (rank != 0) {
        return;
    }

    // Validate input
    std::size_t num_names = names.size();
    if (num_names != num_entries) {
        throw std::runtime_error("Number of names must match number of durations");
    }

    // Convert names to 2D char array [num_entries][max_str_len]
    std::size_t max_str_len = 0;
    for (auto const& name : names) {
        max_str_len = std::max(max_str_len, name.size());
    }
    max_str_len += 1; // null terminator
    
    std::vector<char> timing_names_2d(num_entries * max_str_len, '\0');
    for (std::size_t i = 0; i < num_names; ++i) {
        std::copy(names[i].begin(), names[i].end(), 
                  timing_names_2d.begin() + i * max_str_len);
    }

    // Expose timing data (sizes as scalars, arrays via PDI_expose)
    ddc::expose_to_pdi("timing_values_size", num_entries);
    ddc::expose_to_pdi("timing_names_count", num_entries);
    ddc::expose_to_pdi("timing_names_strlen", max_str_len);
    PDI_expose("timing_values", durations, PDI_OUT);
    PDI_expose("timing_names", timing_names_2d.data(), PDI_OUT);
    ddc::PdiEvent("write_timing");
}

} // namespace

int main(int argc, char** argv)
{

    /**************************************
    *   Main program: mini_app            *
    *  
    ****************************************/
    Kokkos::ScopeGuard scope(argc, argv);
    ddc::ScopeGuard ddc_scope(argc, argv);
    MPI_Init(&argc, &argv);
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    steady_clock::time_point time_points[5];

    print_banner(rank);
    //---------------------------------------------------------
    // Read and initialize the configuration
    //---------------------------------------------------------
    if (rank == 0) {
        cout << "Initializing 5D particle distribution function." << endl;
    }
    ConfigHandles configs = parse_config_files(argc, argv);
    PDI_init(configs.conf_pdi);
    //---------------------------------------------------------
    // Initialisation of the mesh (sp, space, phase-space)
    //---------------------------------------------------------
    MeshInitializationResult const mesh = initialize_mesh(rank, configs.conf_gyselax);
    IdxRangeSp const idx_range_sp = mesh.idx_range_sp;
    IdxRangeSpGrid const meshGridSp = mesh.mesh_sp;
    IdxRangeSpVparMu const meshGridSpVparMu = mesh.mesh_sp_vparmu;
    //---------------------------------------------------------
    // Initialisation of the distribution function
    //---------------------------------------------------------
    time_points[0] = steady_clock::now();
    DFieldMemSpGrid allfdistribu(meshGridSp);
    init_distribution_fun(allfdistribu, meshGridSpVparMu, meshGridSp);
    time_points[1] = steady_clock::now();
    
    //---------------------------------------------------------
    // Read application version from YAML config
    //---------------------------------------------------------
    string const version = PCpp_string(configs.conf_gyselax, ".Application.version");
    if (rank == 0) {
        cout << "Performing Version: " << version << "" << endl;
    }

    if (version == "mpi_transpose") {
        MPITransposeAllToAll<Tor3DSplit, V2DSplit> transpose(meshGridSp, MPI_COMM_WORLD);
    }
    time_points[2] = steady_clock::now();
    //---------------------------------------------------------
    // Write 5D distribution function and coordinates to file using PDI
    //---------------------------------------------------------
    // Create host version of distribution function for I/O (needed for PDI)
    host_t<DFieldMemSpGrid> allfdistribu_host(mesh.mesh_sp);
    ddc::parallel_deepcopy(allfdistribu_host, allfdistribu);
    time_points[3] = steady_clock::now();
    write_fdistribu(rank, mesh, allfdistribu_host);
    time_points[4] = steady_clock::now();
    //---------------------------------------------------------
    // Finalize PDI and MPI
    //---------------------------------------------------------
    if (rank == 0) {
        double durations[5];
        for (int i = 0; i < 4; i++) {
            durations[i] = std::chrono::duration<double>(time_points[i+1] - time_points[i]).count();
        }
        durations[4] = std::chrono::duration<double>(time_points[4] - time_points[0]).count();
        cout << "Time initialization: " << durations[0] << "s" << endl;
        cout << "Time transpose: " << durations[1] << "s" << endl;
        cout << "Time gpu2cpu: " << durations[2] << "s" << endl;
        cout << "Time write: " << durations[3] << "s" << endl;
        cout << "Time total: " << durations[4] << "s" << endl;

        // Use the new function to write timing stats as a table
        std::vector<std::string> timing_names = {
            "initialization",
            "transpose",
            "gpu2cpu",
            "write",
            "total"
        };
        write_cpu_time_stats(rank, durations, timing_names, timing_names.size());
    }
    PDI_finalize();
    MPI_Finalize();

    PC_tree_destroy(&configs.conf_pdi);
    PC_tree_destroy(&configs.conf_gyselax);

    return EXIT_SUCCESS;
}

