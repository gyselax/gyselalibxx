// SPDX-License-Identifier: MIT
#include <cmath>
#include <iostream>

#include <ddc/ddc.hpp>

#include <paraconf.h>
#include <pdi.h>

#include "CollisionSpVparMu.hpp"
#include "bsl_advection_1d.hpp"
#include "collisioninfo_radial.hpp"
#include "ddc_alias_inline_functions.hpp"
#include "euler.hpp"
#include "geometry.hpp"
#include "input.hpp"
#include "mpitransposealltoall.hpp"
#include "neumann_spline_quadrature.hpp"
#include "non_uniform_interpolation_points.hpp"
#include "output.hpp"
#include "paraconfpp.hpp"
#include "pdi_out.yml.hpp"
#include "simpson_quadrature.hpp"
#include "species_info.hpp"
#include "species_init.hpp"
#include "spline_interpolator.hpp"
#include "test_collTokamAxi.yaml.hpp"

using std::cerr;
using std::cout;
using std::endl;
using std::chrono::steady_clock;

// Define aliases for operators
using AdvectionVpar = BslAdvection1D<
        GridVpar,
        IdxRangeVpar,
        IdxRangeSpTor2DV2D,
        SplineVparBuilder_1d,
        SplineVparEvaluator_1d>;

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

    // - Collisions Info
    double const nustar0_rpeak = PCpp_double(conf_gyselax, ".Collisions.nustar0_rpeak");
    std::int8_t const collisions_interspecies = PCpp_bool(conf_gyselax, ".Collisions.interspecies");
    double const rhopeak = PCpp_double(conf_gyselax, ".Collisions.rhopeak");

    // - Algorithm Info
    double const deltat = PCpp_double(conf_gyselax, ".Algorithm.deltat");
    int const nbiter = static_cast<int>(PCpp_int(conf_gyselax, ".Algorithm.nbiter"));

    // - Output info
    double const time_diag = PCpp_double(conf_gyselax, ".Output.time_diag");
    int const nbstep_diag = int(time_diag / deltat);
    ddc::expose_to_pdi("nbstep_diag", nbstep_diag);

    //----------------------------------------------
    // Read information from the HDF5 file read_restart_filename
    //----------------------------------------------
    steady_clock::time_point const start_read = steady_clock::now();
    // - Read species info
    IdxRangeSp idxrange_kinsp = init_kinetic_species();

    // - Read mesh
    // -- Read breakpoints and grid
    IdxRangeR const global_idxrange_r = init_spline_dependent_idx_range<
            GridR,
            BSplinesR,
            SplineInterpPointsR>(conf_gyselax, "tor1");
    IdxRangeTheta const global_idxrange_theta = init_spline_dependent_idx_range<
            GridTheta,
            BSplinesTheta,
            SplineInterpPointsTheta>(conf_gyselax, "tor2");
    IdxRangeVpar const global_idxrange_vpar = init_spline_dependent_idx_range<
            GridVpar,
            BSplinesVpar,
            SplineInterpPointsVpar>(conf_gyselax, "vpar");
    IdxRangeMu const global_idxrange_mu = init_spline_dependent_idx_range<
            GridMu,
            BSplinesMu,
            SplineInterpPointsMu>(conf_gyselax, "mu");
    expose_mesh_to_pdi("breakpoints_tor1", ddc::discrete_space<BSplinesR>().break_point_domain());
    expose_mesh_to_pdi("grid_tor1", global_idxrange_r);
    expose_mesh_to_pdi(
            "breakpoints_tor2",
            ddc::discrete_space<BSplinesTheta>().break_point_domain());
    expose_mesh_to_pdi("grid_tor2", global_idxrange_theta);
    expose_mesh_to_pdi(
            "breakpoints_vpar",
            ddc::discrete_space<BSplinesVpar>().break_point_domain());
    expose_mesh_to_pdi("grid_vpar", global_idxrange_vpar);
    expose_mesh_to_pdi("breakpoints_mu", ddc::discrete_space<BSplinesMu>().break_point_domain());
    expose_mesh_to_pdi("grid_mu", global_idxrange_mu);

    // - Read magnetic configuration
    // -- Read R_matrix and Z_matrix
    // -- Read normB_matrix
    IdxRangeTor2D const global_idxrange_tor2D(global_idxrange_r, global_idxrange_theta);
    DFieldMemTor2D_host R_matrix_host(global_idxrange_tor2D);
    DFieldMemTor2D_host Z_matrix_host(global_idxrange_tor2D);
    DFieldMemTor2D_host normB_matrix_host(global_idxrange_tor2D);
    DFieldMemR_host safety_factor_host(global_idxrange_r);
    ddc::PdiEvent("read_magnetic_config")
            .with("R_matrix", R_matrix_host)
            .with("Z_matrix", Z_matrix_host)
            .with("normB_matrix", normB_matrix_host)
            .with("safety_factor", safety_factor_host);

    // - Read poloidal cross-section of the 3 moments: density, temperature and Upar
    IdxRangeSpTor2D const
            global_idxrange_sptor2D(idxrange_kinsp, global_idxrange_r, global_idxrange_theta);
    DFieldMemSpTor2D_host density_torCS_host(global_idxrange_sptor2D);
    DFieldMemSpTor2D_host temperature_torCS_host(global_idxrange_sptor2D);
    DFieldMemSpTor2D_host Upar_torCS_host(global_idxrange_sptor2D);
    ddc::PdiEvent("read_profiles")
            .with("densityTorCS", density_torCS_host)
            .with("temperatureTorCS", temperature_torCS_host)
            .with("UparTorCS", Upar_torCS_host);

    // - Read the distribution function
    IdxRangeSpTor2DV2D const global_idxrange_sptor2Dv2D(
            idxrange_kinsp,
            global_idxrange_theta,
            global_idxrange_r,
            global_idxrange_vpar,
            global_idxrange_mu);
    MPITransposeAllToAll<Tor2DDistributed, V2DDistributed>
            transpose(global_idxrange_sptor2Dv2D, MPI_COMM_WORLD);
    IdxRangeSpTor2DV2D local_idxrange_sptor2Dv2D(transpose.get_local_idx_range<Tor2DDistributed>());
    IdxRangeSpV2DTor2D local_idxrange_spv2Dtor2D(transpose.get_local_idx_range<V2DDistributed>());
    DFieldMemSpTor2DV2D_host allfdistribu_host(local_idxrange_sptor2Dv2D);
    PDI_expose_idx_range(local_idxrange_sptor2Dv2D, "local_fdistribu");
    double time_saved;
    ddc::PdiEvent("read_fdistribu")
            .with("time_saved", time_saved)
            .and_with("fdistribu", allfdistribu_host);

    steady_clock::time_point const end_read = steady_clock::now();
    double const time_read = std::chrono::duration<double>(end_read - start_read).count();

    //------------------------------------------------------
    // Initialisation of the distribution function layouts
    //------------------------------------------------------
    steady_clock::time_point const start_init = steady_clock::now();
    // allfdistribu on the Tor2DDistributed layout
    auto allfdistribu_vpar_mu_alloc = ddc::create_mirror_view_and_copy(
            Kokkos::DefaultExecutionSpace(),
            get_field(allfdistribu_host));
    DFieldSpTor2DV2D allfdistribu_vpar_mu = get_field(allfdistribu_vpar_mu_alloc);

    // uninitialised allfdistribu on the V2DDistributed layout
    DFieldMemSpV2DTor2D allfdistribu_r_theta_alloc(local_idxrange_spv2Dtor2D);
    DFieldSpV2DTor2D allfdistribu_r_theta = get_field(allfdistribu_r_theta_alloc);

    //----------------------------------------------
    // Initialisation of the splines operators
    //----------------------------------------------
    // Create extrapolation rules for the spline evaluators
    ddc::NullExtrapolationRule extrapol_vpar_min;
    ddc::NullExtrapolationRule extrapol_vpar_max;

    // Create a builder and evaluator for splines in the vpar direction batched over all other directions
    SplineVparBuilder spline_builder_vpar(local_idxrange_sptor2Dv2D);
    SplineVparEvaluator spline_eval_vpar(extrapol_vpar_min, extrapol_vpar_max);

    // Create a builder and evaluator for splines in the vpar direction without a batch direction
    // This is used for the null advection field which only depends on the vpar direction
    SplineVparBuilder_1d spline_builder_adv_field_vpar(global_idxrange_vpar);
    SplineVparEvaluator_1d spline_eval_adv_field_vpar(extrapol_vpar_min, extrapol_vpar_max);

    // Create a spline interpolator in vpar direction
    PreallocatableSplineInterpolator const
            spline_interpolator_vpar(spline_builder_vpar, spline_eval_vpar);

    //---------------------------------------------
    // Initialisation of the advection operators
    //----------------------------------------------
    // Create a timestepper to calculate the foot of the characteristics in the advection
    Euler<FieldMemVpar<CoordVpar>, DFieldMemVpar> characteristic_timestepper(global_idxrange_vpar);

    // uninitialised advection field A(vpar)
    DFieldMemVpar null_advection_field_alloc(global_idxrange_vpar);
    DFieldVpar null_advection_field(get_field(null_advection_field_alloc));
    DFieldMem<IdxRange<ddc::Deriv<Vpar>>> null_advection_field_deriv_vpar_min_alloc(
            spline_builder_adv_field_vpar.batched_derivs_xmin_domain());
    DFieldMem<IdxRange<ddc::Deriv<Vpar>>> null_advection_field_deriv_vpar_max_alloc(
            spline_builder_adv_field_vpar.batched_derivs_xmax_domain());
    DField<IdxRange<ddc::Deriv<Vpar>>> null_advection_field_deriv_vpar_min(
            null_advection_field_deriv_vpar_min_alloc);
    DField<IdxRange<ddc::Deriv<Vpar>>> null_advection_field_deriv_vpar_max(
            null_advection_field_deriv_vpar_max_alloc);
    // initialise the advection field A(vpar) = 0.0
    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            global_idxrange_vpar,
            KOKKOS_LAMBDA(IdxVpar idx) { null_advection_field(idx) = 0.0; });
    // initialise the derivative of the advection field dA(vpar) = 0.0 at the boundaries
    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            get_idx_range(null_advection_field_deriv_vpar_min),
            KOKKOS_LAMBDA(Idx<ddc::Deriv<Vpar>> idx) {
                null_advection_field_deriv_vpar_min(idx) = 0.0;
            });
    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            get_idx_range(null_advection_field_deriv_vpar_max),
            KOKKOS_LAMBDA(Idx<ddc::Deriv<Vpar>> idx) {
                null_advection_field_deriv_vpar_max(idx) = 0.0;
            });

    // Create a vpar advection operator
    AdvectionVpar vpar_advection(
            spline_interpolator_vpar,
            spline_builder_adv_field_vpar,
            spline_eval_adv_field_vpar,
            characteristic_timestepper);

    //----------------------------------------------
    // Initialisation of the collision operator
    //----------------------------------------------
    // Create and initialise a collision operator
    int const Nr = global_idxrange_r.size();
    IdxR idx_r(int(Nr * rhopeak));
    double const rpeak = ddc::coordinate(idx_r);
    double const q_rpeak = safety_factor_host(idx_r);
    double const R0 = R_matrix_host(global_idxrange_r.front(), global_idxrange_theta.front());
    IdxRangeR local_idxrange_r(local_idxrange_sptor2Dv2D);
    CollisionInfoRadial<GridR> collision_info(
            rpeak,
            q_rpeak,
            R0,
            nustar0_rpeak,
            collisions_interspecies,
            local_idxrange_r);

    DFieldMemVpar const coeff_intdvpar(
            neumann_spline_quadrature_coefficients<
                    Kokkos::DefaultExecutionSpace>(global_idxrange_vpar, spline_builder_vpar));
    DFieldMemMu const coeff_intdmu(
            simpson_trapezoid_quadrature_coefficients_1d<
                    Kokkos::DefaultExecutionSpace>(global_idxrange_mu, Extremity::BACK));
    // --> Initialisation of the magnetic configuration: B norm and Bstar_s
    IdxRangeTor2D local_idxrange_tor2D(global_idxrange_tor2D);
    DFieldMemTor2D B_norm_alloc(local_idxrange_tor2D);
    DFieldTor2D B_norm(get_field(B_norm_alloc));
    ddc::parallel_deepcopy(B_norm, normB_matrix_host[local_idxrange_tor2D]);
    // [TODO]: Bstar_s is forced equal to B but in the future
    //         Bstar_s should be an input variable of the constructor
    IdxRangeSpTor2DVpar idxrange_Bstar(local_idxrange_sptor2Dv2D);
    DFieldMemSpTor2DVpar Bstar_s_alloc(idxrange_Bstar);
    DFieldSpTor2DVpar Bstar_s(get_field(Bstar_s_alloc));
    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            idxrange_Bstar,
            KOKKOS_LAMBDA(IdxSpTor2DVpar const idx_Bstar) {
                IdxTor2D idx_thetar(idx_Bstar);
                Bstar_s(idx_Bstar) = B_norm(idx_thetar);
            });
    CollisionSpVparMu<
            CollisionInfoRadial<GridR>,
            IdxRangeSpTor2DV2D,
            GridVpar,
            GridMu,
            DConstFieldTor2D>
            collision_operator(
                    collision_info,
                    local_idxrange_sptor2Dv2D,
                    get_const_field(coeff_intdmu),
                    get_const_field(coeff_intdvpar),
                    get_const_field(B_norm),
                    get_const_field(Bstar_s));

    steady_clock::time_point const end_init = steady_clock::now();
    double const time_init = std::chrono::duration<double>(end_init - start_init).count();

    //----------------------------------------------
    // Algorithm
    //----------------------------------------------
    // ---> Algorithm time loop
    long int iter_globsimu(0); // [TODO] Read this value from HDF5 file
    double time_globsimu(0.0); // [TODO] Read this value from HDF5 file
    ddc::expose_to_pdi("iter_start", iter_globsimu);

    steady_clock::time_point const start_algo = steady_clock::now();

    int iter_run = 0;
    for (; iter_run < nbiter; ++iter_run) {
        cout << "iter_run = " << iter_run << " ; time_simu = " << time_globsimu << endl;

        // Advection in vpar direction on dt/0.5
        // [TODO] Replace the null advection by the good one
        vpar_advection(
                allfdistribu_vpar_mu,
                null_advection_field,
                0.5 * deltat,
                null_advection_field_deriv_vpar_min,
                null_advection_field_deriv_vpar_max);

        // Get the values of allfdistribu on the V2DDistributed layout
        transpose(
                Kokkos::DefaultExecutionSpace(),
                allfdistribu_r_theta,
                get_const_field(allfdistribu_vpar_mu));

        // [TODO] Advection in (r,theta) direction on dt

        // Get the values of allfdistribu on the Tor2DDistributed layout
        transpose(
                Kokkos::DefaultExecutionSpace(),
                allfdistribu_vpar_mu,
                get_const_field(allfdistribu_r_theta));

        // Collision operator on dt
        collision_operator(allfdistribu_vpar_mu, deltat);

        // Advection in vpar direction on dt/0.5
        // [TODO] Replace the null advection by the good one
        vpar_advection(
                allfdistribu_vpar_mu,
                null_advection_field,
                0.5 * deltat,
                null_advection_field_deriv_vpar_min,
                null_advection_field_deriv_vpar_max);

        // [TODO] Compute the three fluid moments

        time_globsimu = time_globsimu + deltat;
        iter_globsimu = iter_globsimu + 1;
    }

    steady_clock::time_point const end_algo = steady_clock::now();
    double const time_algo = std::chrono::duration<double>(end_algo - start_algo).count();

    //----------------------------------------------
    // Save the HDF5 restart file write_restart_filename
    //----------------------------------------------
    steady_clock::time_point const start_write = steady_clock::now();
    ddc::parallel_deepcopy(allfdistribu_host, allfdistribu_vpar_mu);
    ddc::PdiEvent("write_restart")
            .with("iter_saved", iter_globsimu)
            .with("time_saved", time_globsimu)
            .with("densityTorCS", density_torCS_host)
            .with("temperatureTorCS", temperature_torCS_host)
            .with("UparTorCS", Upar_torCS_host)
            .with("fdistribu", allfdistribu_host);
    steady_clock::time_point const end_write = steady_clock::now();
    double const time_write = std::chrono::duration<double>(end_write - start_write).count();

    std::cout << " - Reading time : " << time_read << "s\n";
    std::cout << " - Initialisation time : " << time_init << "s\n";
    std::cout << " - Algorithm time : " << time_algo << "s\n";
    std::cout << " - Writing time : " << time_write << "s\n";

    PDI_finalize();

    MPI_Finalize();

    PC_tree_destroy(&conf_pdi);

    PC_tree_destroy(&conf_gyselax);

    return EXIT_SUCCESS;
}
