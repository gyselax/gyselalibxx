// SPDX-License-Identifier: MIT
#include <cmath>
#include <iostream>

#include <ddc/ddc.hpp>

#include <paraconf.h>
#include <pdi.h>

#include "ddc_alias_inline_functions.hpp"
#include "geometry.hpp"
#include "input.hpp"
#include "mpitransposealltoall.hpp"
#include "non_uniform_interpolation_points.hpp"
#include "output.hpp"
#include "paraconfpp.hpp"
#include "pdi_out.yml.hpp"
#include "test_collTokamAxi.yaml.hpp"

using std::cerr;
using std::cout;
using std::endl;
using std::chrono::steady_clock;

int main(int argc, char** argv)
{
    long int iter_start;
    PC_tree_t conf_gyselax;
    parse_executable_arguments(conf_gyselax, iter_start, argc, argv, params_yaml);
    PC_tree_t conf_pdi = PC_parse_string(PDI_CFG);
    PC_errhandler(PC_NULL_HANDLER);

    Kokkos::ScopeGuard scope(argc, argv);
    ddc::ScopeGuard ddc_scope(argc, argv);
    MPI_Init(&argc, &argv);
    PDI_init(conf_pdi);

    ddc::expose_to_pdi("iter_start", iter_start);

    //----------------------------------------------
    // Extract information from the YAML file
    //----------------------------------------------
    //  - Input and output file names info
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

    // - Algorithm Info
    double const deltat = PCpp_double(conf_gyselax, ".Algorithm.deltat");
    //int const nbiter = static_cast<int>(PCpp_int(conf_gyselax, ".Algorithm.nbiter"));

    // - Output info
    //double const time_diag = PCpp_double(conf_gyselax, ".Output.time_diag");
    //int const nbstep_diag = int(time_diag / deltat);

    //----------------------------------------------
    // Read information from the HDF5 file read_restart_filename
    //----------------------------------------------
    steady_clock::time_point const start_read = steady_clock::now();
    // - Read species info
    std::vector<int> species;
    std::vector<double> charges;
    std::vector<double> masses;
    PDI_get_arrays("read_species", "species", species, "charges", charges, "masses", masses);
    IdxStepSp const kinspecies(charges.size());
    IdxRangeSp const idxrange_kinsp(IdxSp(0), kinspecies);

    // - Read mesh
    // -- Read breakpoints and grid
    IdxRangeR const idxrange_r = init_spline_dependent_idx_range<
            GridR,
            BSplinesR,
            SplineInterpPointsR>(conf_gyselax, "tor1");
    IdxRangeTheta const idxrange_theta = init_spline_dependent_idx_range<
            GridTheta,
            BSplinesTheta,
            SplineInterpPointsTheta>(conf_gyselax, "tor2");
    IdxRangeVpar const idxrange_vpar = init_spline_dependent_idx_range<
            GridVpar,
            BSplinesVpar,
            SplineInterpPointsVpar>(conf_gyselax, "vpar");
    IdxRangeMu const idxrange_mu = init_spline_dependent_idx_range<
            GridMu,
            BSplinesMu,
            SplineInterpPointsMu>(conf_gyselax, "mu");
    expose_mesh_to_pdi("breakpoints_tor1", ddc::discrete_space<BSplinesR>().break_point_domain());
    expose_mesh_to_pdi("grid_tor1", idxrange_r);
    expose_mesh_to_pdi(
            "breakpoints_tor2",
            ddc::discrete_space<BSplinesTheta>().break_point_domain());
    expose_mesh_to_pdi("grid_tor2", idxrange_theta);
    expose_mesh_to_pdi(
            "breakpoints_vpar",
            ddc::discrete_space<BSplinesVpar>().break_point_domain());
    expose_mesh_to_pdi("grid_vpar", idxrange_vpar);
    expose_mesh_to_pdi("breakpoints_mu", ddc::discrete_space<BSplinesMu>().break_point_domain());
    expose_mesh_to_pdi("grid_mu", idxrange_mu);

    // - Read magnetic configuration
    // -- Read R_matrix and Z_matrix
    // -- Read normB_matrix
    IdxRangeTor2D const idxrange_tor2D(idxrange_r, idxrange_theta);
    DFieldMemTor2D_host R_matrix_host(idxrange_tor2D);
    DFieldMemTor2D_host Z_matrix_host(idxrange_tor2D);
    DFieldMemTor2D_host normB_matrix_host(idxrange_tor2D);
    ddc::PdiEvent("read_magnetic_config")
            .with("R_matrix", R_matrix_host)
            .with("Z_matrix", Z_matrix_host)
            .with("normB_matrix", normB_matrix_host);

    // - Read poloidal cross-section of the 3 moments: density, temperature and Upar
    IdxRangeSpTor2D const idxrange_sptor2D(idxrange_kinsp, idxrange_r, idxrange_theta);
    DFieldMemSpTor2D_host density_torCS_host(idxrange_sptor2D);
    DFieldMemSpTor2D_host temperature_torCS_host(idxrange_sptor2D);
    DFieldMemSpTor2D_host Upar_torCS_host(idxrange_sptor2D);
    ddc::PdiEvent("read_profiles")
            .with("densityTorCS", density_torCS_host)
            .with("temperatureTorCS", temperature_torCS_host)
            .with("UparTorCS", Upar_torCS_host);

    // - Read the distribution function
    IdxRangeSpTor2DV2D const idxrange_sptor2Dv2D(
            idxrange_kinsp,
            idxrange_theta,
            idxrange_r,
            idxrange_vpar,
            idxrange_mu);
    MPITransposeAllToAll<Tor2DDistributed, V2DDistributed>
            transpose(idxrange_sptor2Dv2D, MPI_COMM_WORLD);
    DFieldMemSpTor2DV2D_host allfdistribu_host(transpose.get_local_idx_range<Tor2DDistributed>());
    PDI_expose_idx_range(get_idx_range(allfdistribu_host), "local_fdistribu");
    double time_saved;
    ddc::PdiEvent("read_fdistribu")
            .with("time_saved", time_saved)
            .and_with("fdistribu", allfdistribu_host);

    steady_clock::time_point const end_read = steady_clock::now();
    double const time_read = std::chrono::duration<double>(end_read - start_read).count();

    //----------------------------------------------
    // Algorithm
    //----------------------------------------------
    steady_clock::time_point const start_algo = steady_clock::now();

    // allfdistribu on the Tor2DDistributed layout
    auto allfdistribu_vpar_mu_alloc = ddc::create_mirror_view_and_copy(
            Kokkos::DefaultExecutionSpace(),
            get_field(allfdistribu_host));
    DFieldSpTor2DV2D allfdistribu_vpar_mu = get_field(allfdistribu_vpar_mu_alloc);

    // unintialised allfdistribu on the V2DDistributed layout
    DFieldMemSpV2DTor2D allfdistribu_r_theta_alloc(transpose.get_local_idx_range<V2DDistributed>());
    DFieldSpV2DTor2D allfdistribu_r_theta = get_field(allfdistribu_r_theta_alloc);

    // [TODO] Apply collision operator

    // Get the values of allfdistribu on the V2DDistributed layout
    transpose(
            Kokkos::DefaultExecutionSpace(),
            allfdistribu_r_theta,
            get_const_field(allfdistribu_vpar_mu));

    // [TODO] Compute the three fluid moments

    // Get the values of allfdistribu on the Tor2DDistributed layout
    transpose(
            Kokkos::DefaultExecutionSpace(),
            allfdistribu_vpar_mu,
            get_const_field(allfdistribu_r_theta));

    long int iter_saved(iter_start + 1);
    int iter = 1;
    time_saved = time_saved + iter * deltat;

    steady_clock::time_point const end_algo = steady_clock::now();
    double const time_algo = std::chrono::duration<double>(end_algo - start_algo).count();

    //----------------------------------------------
    // Save the HDF5 restart file write_restart_filename
    //----------------------------------------------
    steady_clock::time_point const start_write = steady_clock::now();
    ddc::parallel_deepcopy(allfdistribu_host, allfdistribu_vpar_mu);
    ddc::PdiEvent("write_restart")
            .with("iter_saved", iter_saved)
            .with("time_saved", time_saved)
            .with("densityTorCS", density_torCS_host)
            .with("temperatureTorCS", temperature_torCS_host)
            .with("UparTorCS", Upar_torCS_host)
            .with("fdistribu", allfdistribu_host);
    steady_clock::time_point const end_write = steady_clock::now();
    double const time_write = std::chrono::duration<double>(end_write - start_write).count();

    std::cout << " - Reading time : " << time_read << "s\n";
    std::cout << " - Algorithm time : " << time_algo << "s\n";
    std::cout << " - Writing time : " << time_write << "s\n";

    PDI_finalize();

    MPI_Finalize();

    PC_tree_destroy(&conf_pdi);

    PC_tree_destroy(&conf_gyselax);

    return EXIT_SUCCESS;
}
