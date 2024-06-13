// SPDX-License-Identifier: MIT

#include <chrono>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>

#include <ddc/ddc.hpp>

#include <paraconf.h>
#include <pdi.h>

#include "CollisionSpVparMu.hpp"
#include "geometry.hpp"
#include "paraconfpp.hpp"
#include "pdi_out.yml.hpp"
#include "simpson_quadrature.hpp"
#include "testcollisions.yaml.hpp"

using std::cerr;
using std::cout;
using std::endl;
using std::chrono::steady_clock;
namespace fs = std::filesystem;

int main(int argc, char** argv)
{
    Kokkos::ScopeGuard scope(argc, argv);
    ddc::ScopeGuard ddc_scope(argc, argv);
    CollisionsGuard a_collision_guard {8, 0};

    long int iter_start(0);
    PC_tree_t conf_gyselax;
    if (argc == 2) {
        conf_gyselax = PC_parse_path(fs::path(argv[1]).c_str());
    } else if (argc == 3) {
        if (argv[1] == std::string_view("--dump-config")) {
            std::fstream file(argv[2], std::fstream::out);
            file << params_yaml;
            return EXIT_SUCCESS;
        }
    } else if (argc == 4) {
        if (argv[1] == std::string_view("--iter-restart")) {
            iter_start = std::strtol(argv[2], NULL, 10);
            conf_gyselax = PC_parse_path(fs::path(argv[3]).c_str());
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

    ddc::expose_to_pdi("iter_start", iter_start);

    // Input and output file names info
    std::string readRestartFileName(PCpp_string(conf_gyselax, ".InputFileNames.read_restart"));
    cout << "Input read restart: " << readRestartFileName << endl;
    std::string writeRestartFileName(PCpp_string(conf_gyselax, ".InputFileNames.write_restart"));
    cout << "Input write restart: " << writeRestartFileName << endl;
    int64_t read_restart_filename_size = readRestartFileName.size();
    int64_t write_restart_filename_size = writeRestartFileName.size();
    PDI_multi_expose(
            "restartFile",
            "read_restart_filename_size",
            &read_restart_filename_size,
            PDI_OUT,
            "read_restart_filename",
            readRestartFileName.c_str(),
            PDI_OUT,
            "write_restart_filename_size",
            &write_restart_filename_size,
            PDI_OUT,
            "write_restart_filename",
            writeRestartFileName.c_str(),
            PDI_OUT,
            NULL);

    std::vector<size_t> grid_tor1_extents(1);
    std::vector<size_t> grid_tor2_extents(1);
    std::vector<size_t> grid_tor3_extents(1);
    std::vector<size_t> grid_vpar_extents(1);
    std::vector<size_t> grid_mu_extents(1);
    std::vector<size_t> species_extents(1);
    std::vector<size_t> charges_extents(1);
    std::vector<size_t> masses_extents(1);
    PDI_multi_expose(
            "read_grid_extents",
            "grid_tor1_extents",
            grid_tor1_extents.data(),
            PDI_INOUT,
            "grid_tor2_extents",
            grid_tor2_extents.data(),
            PDI_INOUT,
            "grid_tor3_extents",
            grid_tor3_extents.data(),
            PDI_INOUT,
            "grid_vpar_extents",
            grid_vpar_extents.data(),
            PDI_INOUT,
            "grid_mu_extents",
            grid_mu_extents.data(),
            PDI_INOUT,
            "species_extents",
            species_extents.data(),
            PDI_INOUT,
            "charges_extents",
            charges_extents.data(),
            PDI_INOUT,
            "masses_extents",
            masses_extents.data(),
            PDI_INOUT,
            NULL);
    std::vector<double> grid_tor1(grid_tor1_extents[0]);
    std::vector<double> grid_tor2(grid_tor2_extents[0]);
    std::vector<double> grid_tor3(grid_tor3_extents[0]);
    std::vector<double> grid_vpar(grid_vpar_extents[0]);
    std::vector<double> grid_mu(grid_mu_extents[0]);
    std::vector<int> species(species_extents[0]);
    std::vector<double> charges(charges_extents[0]);
    std::vector<double> masses(masses_extents[0]);
    PDI_multi_expose(
            "read_grid",
            "grid_tor1",
            grid_tor1.data(),
            PDI_INOUT,
            "grid_tor2",
            grid_tor2.data(),
            PDI_INOUT,
            "grid_tor3",
            grid_tor3.data(),
            PDI_INOUT,
            "grid_vpar",
            grid_vpar.data(),
            PDI_INOUT,
            "grid_mu",
            grid_mu.data(),
            PDI_INOUT,
            "species",
            species.data(),
            PDI_INOUT,
            "charges",
            charges.data(),
            PDI_INOUT,
            "masses",
            masses.data(),
            PDI_INOUT,
            NULL);
    ddc::init_discrete_space<GridTor1>(grid_tor1);
    DDomTor1 dom_tor1(IdxTor1(0), DVecTor1(grid_tor1.size()));
    ddc::init_discrete_space<GridTor2>(grid_tor2);
    DDomTor2 dom_tor2(IdxTor2(0), DVecTor2(grid_tor2.size()));
    ddc::init_discrete_space<GridTor3>(grid_tor3);
    DDomTor3 dom_tor3(IdxTor3(0), DVecTor3(grid_tor3.size()));
    ddc::init_discrete_space<GridVpar>(grid_vpar);
    DDomVpar dom_vpar(IdxVpar(0), DVecVpar(grid_vpar.size()));
    ddc::init_discrete_space<GridMu>(grid_mu);
    DDomMu dom_mu(IdxMu(0), DVecMu(grid_mu.size()));
    DVecSp const kinspecies(charges.size());
    DDomSp const dom_kinsp(IdxSp(0), kinspecies);

    DFieldTor1 field_grid_tor1(dom_tor1);
    ddc::parallel_deepcopy(field_grid_tor1, DViewTor1(grid_tor1.data(), dom_tor1));
    auto field_grid_tor1_host = ddc::create_mirror_view_and_copy(field_grid_tor1.span_view());
    DFieldTor2 field_grid_tor2(dom_tor2);
    ddc::parallel_deepcopy(field_grid_tor2, DViewTor2(grid_tor2.data(), dom_tor2));
    auto field_grid_tor2_host = ddc::create_mirror_view_and_copy(field_grid_tor2.span_view());
    DFieldTor3 field_grid_tor3(dom_tor3);
    ddc::parallel_deepcopy(field_grid_tor3, DViewTor3(grid_tor3.data(), dom_tor3));
    auto field_grid_tor3_host = ddc::create_mirror_view_and_copy(field_grid_tor3.span_view());
    DFieldVpar field_grid_vpar(dom_vpar);
    ddc::parallel_deepcopy(field_grid_vpar, DViewVpar(grid_vpar.data(), dom_vpar));
    auto field_grid_vpar_host = ddc::create_mirror_view_and_copy(field_grid_vpar.span_view());
    DFieldMu field_grid_mu(dom_mu);
    ddc::parallel_deepcopy(field_grid_mu, DViewMu(grid_mu.data(), dom_mu));
    auto field_grid_mu_host = ddc::create_mirror_view_and_copy(field_grid_mu.span_view());
    FieldSp<int> field_species(dom_kinsp);
    ddc::parallel_deepcopy(field_species, ViewSp<int>(species.data(), dom_kinsp));
    auto field_species_host = ddc::create_mirror_view_and_copy(field_species.span_view());
    DFieldSp field_charges(dom_kinsp);
    ddc::parallel_deepcopy(field_charges, DViewSp(charges.data(), dom_kinsp));
    auto field_charges_host = ddc::create_mirror_view_and_copy(field_charges.span_view());
    DFieldSp field_masses(dom_kinsp);
    ddc::parallel_deepcopy(field_masses, DViewSp(masses.data(), dom_kinsp));
    auto field_masses_host = ddc::create_mirror_view_and_copy(field_masses.span_view());

    // Algorithm Info
    double const deltat = PCpp_double(conf_gyselax, ".Algorithm.deltat");
    int const nbiter = static_cast<int>(PCpp_int(conf_gyselax, ".Algorithm.nbiter"));

    // Output info
    double const time_diag = PCpp_double(conf_gyselax, ".Output.time_diag");
    int const nbstep_diag = int(time_diag / deltat);

    cout << "nbiter = " << nbiter << " nbstep_diag = " << nbstep_diag << endl;

    // Poloidal cross-section of the 3 moments: density, temperature and Upar
    DDomSpTorCS const dom_sp_torCS(dom_kinsp, dom_tor2, dom_tor1);
    DFieldSpTorCS_host density_torCS_host(dom_sp_torCS);
    DFieldSpTorCS_host temperature_torCS_host(dom_sp_torCS);
    DFieldSpTorCS_host Upar_torCS_host(dom_sp_torCS);
    ddc::PdiEvent("read_profiles")
            .with("densityTorCS", density_torCS_host)
            .with("temperatureTorCS", temperature_torCS_host)
            .with("UparTorCS", Upar_torCS_host);

    // fdistribu
    DDomSpTor3DV2D const
            dom_sp_tor3D_v2D(dom_kinsp, dom_tor3, dom_tor2, dom_tor1, dom_vpar, dom_mu);
    DFieldSpTor3DV2D_host allfdistribu_host(dom_sp_tor3D_v2D);
    double time_saved;
    ddc::PdiEvent("read_fdistribu")
            .with("time_saved", time_saved)
            .and_with("fdistribu", allfdistribu_host);
    cout << "Reading of time " << time_saved << endl;
    auto allfdistribu_alloc = ddc::create_mirror_view_and_copy(
            Kokkos::DefaultExecutionSpace(),
            allfdistribu_host.span_view());
    auto allfdistribu = allfdistribu_alloc.span_view();

    // Collision operator initialisation
    DFieldTor1 nustar0_r(dom_tor1);
    ddc::parallel_fill(nustar0_r, 0.0); //ATTENTION: Must be changed
    DDomTorCS dom_tor1_tor2(dom_tor2, dom_tor1);
    DFieldTorCS B_norm(dom_tor1_tor2);
    ddc::parallel_fill(B_norm, 1.0);
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
            dom_sp_tor3D_v2D,
            coeff_intdmu.span_cview(),
            coeff_intdvpar.span_cview(),
            nustar0_r.span_view(),
            B_norm.span_view());

    steady_clock::time_point const start = steady_clock::now();

    collision_operator(allfdistribu, deltat);

    long int iter_saved(iter_start + 1);
    int iter = 1;
    time_saved = time_saved + iter * deltat;
    cout << "iter_saved = " << iter_saved << " ; time_saved = " << time_saved << endl;
    ddc::parallel_deepcopy(allfdistribu_host, allfdistribu);
    ddc::PdiEvent("write_restart")
            .with("iter_saved", iter_saved)
            .with("time_saved", time_saved)
            .with("grid_tor1", field_grid_tor1_host)
            .with("grid_tor2", field_grid_tor2_host)
            .with("grid_tor3", field_grid_tor3_host)
            .with("grid_vpar", field_grid_vpar_host)
            .with("grid_mu", field_grid_mu_host)
            .with("species", field_species_host)
            .with("masses", field_masses_host)
            .with("charges", field_charges_host)
            .with("densityTorCS", density_torCS_host)
            .with("temperatureTorCS", temperature_torCS_host)
            .with("UparTorCS", Upar_torCS_host)
            .with("fdistribu", allfdistribu_host);

    steady_clock::time_point const end = steady_clock::now();

    double const simulation_time = std::chrono::duration<double>(end - start).count();
    std::cout << "Simulation time: " << simulation_time << "s\n";

    PDI_finalize();

    PC_tree_destroy(&conf_pdi);

    PC_tree_destroy(&conf_gyselax);

    return EXIT_SUCCESS;
}
