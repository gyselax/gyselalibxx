// SPDX-License-Identifier: MIT
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string_view>

#include <ddc/ddc.hpp>
#include <ddc/pdi.hpp>

#include <paraconf.h>
#include <pdi.h>

#include "bsl_advection_vx.hpp"
#include "bsl_advection_3d_rot_exact_splitting.hpp"
#include "bsl_advection_x.hpp"
#include "ddc_alias_inline_functions.hpp"
#include "fft_hybrid_solver_1d.hpp"
#include "geometry.hpp"
#include "input.hpp"
#include "maxwellianequilibrium.hpp"
#include "mpimomentscalculator.hpp"
#include "mpisplithybridvlasovsolver.hpp"
//#include "mpichargedensitycalculator.hpp"
#include "mpisplitvlasovsolver.hpp"
#include "mpitransposealltoall.hpp"
#include "neumann_spline_quadrature.hpp"
#include "output.hpp"
#include "paraconfpp.hpp"
#include "params.yaml.hpp"
//#include "qnsolver.hpp"
#include "pdi_out.yml.hpp"
#include "hybridsplitting.hpp"
#include "hybridfieldsolver.hpp"
#include "singlemodeperturbinitialisation.hpp"
#include "hybridmodel_field_initialisation.hpp"
#include "species_info.hpp"
#include "species_init.hpp"
#include "spline_interpolator.hpp"

using std::cerr;
using std::endl;
using std::chrono::steady_clock;
namespace fs = std::filesystem;

int main(int argc, char** argv)
{
    PC_tree_t conf_gyselalibxx = parse_executable_arguments(argc, argv, params_yaml);
    PC_tree_t conf_pdi = PC_parse_string(PDI_CFG);
    PC_errhandler(PC_NULL_HANDLER);
    MPI_Init(&argc, &argv);
    PDI_init(conf_pdi);

    Kokkos::ScopeGuard kokkos_scope(argc, argv);
    ddc::ScopeGuard ddc_scope(argc, argv);

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // Reading config
    // --> Mesh info
    IdxRangeX const idxrange_x = init_spline_dependent_idx_range<
            GridX,
            BSplinesX,
            SplineInterpPointsX>(conf_gyselalibxx, "x");

    IdxRangeVx const idxrange_vx = init_spline_dependent_idx_range<
            GridVx,
            BSplinesVx,
            SplineInterpPointsVx>(conf_gyselalibxx, "vx");
    IdxRangeVy const idxrange_vy = init_spline_dependent_idx_range<
            GridVy,
            BSplinesVy,
            SplineInterpPointsVy>(conf_gyselalibxx, "vy");
    IdxRangeVz const idxrange_vz = init_spline_dependent_idx_range<
            GridVz,
            BSplinesVz,
            SplineInterpPointsVz>(conf_gyselalibxx, "vz");
 
    IdxRangeVxVyVz idxrange_vxvyvz(idxrange_vx, idxrange_vy, idxrange_vz);
    IdxRangeXVxVyVz const idxrange_xvxvyvz(idxrange_x, idxrange_vx, idxrange_vy, idxrange_vz);

    IdxRangeSp const idx_range_kinsp = init_species(conf_gyselalibxx);

    // set the index range of charge and mean velocity dnesity of each kind of ions.
    IdxRangeSpX const idxrange_spx(idx_range_kinsp, idxrange_x);

    IdxRangeSpXVxVyVz const idxrange_glob_spxvxvyvz(idx_range_kinsp, idxrange_xvxvyvz);

    MPITransposeAllToAll<X1DSplit, V3DSplit> transpose(idxrange_glob_spxvxvyvz, MPI_COMM_WORLD);

    IdxRangeSpXVxVyVz idxrange_spxvxvyvz_x1Dsplit(transpose.get_local_idx_range<X1DSplit>());
    IdxRangeSpVxVyVzX idxrange_spvxvyvzx_v3Dsplit(transpose.get_local_idx_range<V3DSplit>());

    IdxRangeVxVyVz idxrange_vxvyvz_v3Dsplit(idxrange_spvxvyvzx_v3Dsplit);
    IdxRangeVz idxrange_vz_v3Dsplit(idxrange_spvxvyvzx_v3Dsplit);
    IdxRangeVxVyVz idxrange_vxvyvz_x1Dsplit(idxrange_spxvxvyvz_x1Dsplit);

    IdxRangeVxVyVzX idxrange_vxvyvzx_v3Dsplit(idxrange_spvxvyvzx_v3Dsplit);
    IdxRangeXVxVyVz idxrange_xvxvyvz_x1Dsplit(idxrange_spxvxvyvz_x1Dsplit);
    SplineXBuilder const builder_x(idxrange_vxvyvzx_v3Dsplit);
    SplineVxBuilder const builder_vx(idxrange_xvxvyvz_x1Dsplit);
    SplineVyBuilder const builder_vy(idxrange_xvxvyvz_x1Dsplit);
    SplineVzBuilder const builder_vz(idxrange_xvxvyvz_x1Dsplit);

    IdxRangeSpVxVyVz idxrange_spvxvyvz_local(idxrange_spxvxvyvz_x1Dsplit);
    // Initialisation of the distribution function
    DFieldMemSpVxVyVz allfequilibrium(idxrange_spvxvyvz_local);
    MaxwellianEquilibrium const init_fequilibrium
            = MaxwellianEquilibrium::init_from_input(idx_range_kinsp, conf_gyselalibxx);
    init_fequilibrium(get_field(allfequilibrium));
    DFieldMemSpXVxVyVz allfdistribu_x1D_split(idxrange_spxvxvyvz_x1Dsplit);
    DFieldMemSpVxVyVzX allfdistribu_v3D_split(idxrange_spvxvyvzx_v3Dsplit);
    SingleModePerturbInitialisation const init = SingleModePerturbInitialisation::
            init_from_input(get_const_field(allfequilibrium), idx_range_kinsp, conf_gyselalibxx);
    init(get_field(allfdistribu_x1D_split));

    // Define the date structure of the mean velocities and charge density for each kind of ions.
    DFieldMemSpX mean_current_x_each(idxrange_spx);
    DFieldMemSpX mean_current_y_each(idxrange_spx);
    DFieldMemSpX mean_current_z_each(idxrange_spx);
    DFieldMemSpX rho_each(idxrange_spx);
    // Define the date structure of the mean velocity and charge densities for all kinds of ions
    DFieldMemX mean_current_x(idxrange_x);
    DFieldMemX mean_current_y(idxrange_x);
    DFieldMemX mean_current_z(idxrange_x);
    DFieldMemX rho(idxrange_x);
    // Define the kinetic density 
    DFieldMemX kinetic(idxrange_x);
    // Define the date structure of the magnetic field and the pressure
    DFieldMemX magnetic_field_x(idxrange_x);
    DFieldMemX magnetic_field_y(idxrange_x);
    DFieldMemX magnetic_field_z(idxrange_x);
    DFieldMemX pressure(idxrange_x);
    int magnetic_mode = 1;
    double magnetic_amplitude = 0.1;
    int pressure_mode = 1;
    double pressure_amplitude = 0.1;
    Hybridmodel_field_initialisation const init_fields(magnetic_mode, magnetic_amplitude, pressure_mode, pressure_amplitude);
    init_fields(get_field(magnetic_field_x), get_field(magnetic_field_y),get_field(magnetic_field_z), get_field(pressure));
    // momentum 
    DFieldMemX momentum_x(idxrange_x);
    DFieldMemX momentum_y(idxrange_x);
    DFieldMemX momentum_z(idxrange_x);

    // --> Algorithm info
    double const deltat = PCpp_double(conf_gyselalibxx, ".Algorithm.deltat");
    int const nbiter = static_cast<int>(PCpp_int(conf_gyselalibxx, ".Algorithm.nbiter"));

    // --> Output info
    double const time_diag = PCpp_double(conf_gyselalibxx, ".Output.time_diag");
    int const nbstep_diag = int(time_diag / deltat);

    // Create spline evaluator
    ddc::PeriodicExtrapolationRule<X> bv_x_min;
    ddc::PeriodicExtrapolationRule<X> bv_x_max;
    SplineXEvaluator const spline_x_evaluator(bv_x_min, bv_x_max);

    PreallocatableSplineInterpolator const spline_x_interpolator(builder_x, spline_x_evaluator, idxrange_vxvyvzx_v3Dsplit);

    ddc::ConstantExtrapolationRule<Vx> bv_vx_min(ddc::coordinate(idxrange_vx.front()));
    ddc::ConstantExtrapolationRule<Vx> bv_vx_max(ddc::coordinate(idxrange_vx.back()));
    SplineVxEvaluator const spline_vx_evaluator(bv_vx_min, bv_vx_max);

    PreallocatableSplineInterpolator const spline_vx_interpolator(builder_vx, spline_vx_evaluator, idxrange_xvxvyvz_x1Dsplit);

    ddc::ConstantExtrapolationRule<Vy> bv_vy_min(ddc::coordinate(idxrange_vy.front()));
    ddc::ConstantExtrapolationRule<Vy> bv_vy_max(ddc::coordinate(idxrange_vy.back()));
    SplineVyEvaluator const spline_vy_evaluator(bv_vy_min, bv_vy_max);

    PreallocatableSplineInterpolator const spline_vy_interpolator(builder_vy, spline_vy_evaluator, idxrange_xvxvyvz_x1Dsplit);

    ddc::ConstantExtrapolationRule<Vz> bv_vz_min(ddc::coordinate(idxrange_vz.front()));
    ddc::ConstantExtrapolationRule<Vz> bv_vz_max(ddc::coordinate(idxrange_vz.back()));
    SplineVzEvaluator const spline_vz_evaluator(bv_vz_min, bv_vz_max);

    PreallocatableSplineInterpolator const spline_vz_interpolator(builder_vz, spline_vz_evaluator, idxrange_xvxvyvz_x1Dsplit);

    // ============== spline evaluator purely on 3D velocity ==================
    SplineVxEvaluator const spline_exact_vx_evaluator(bv_vx_min, bv_vx_max);
    SplineVxBuilder const builder_exact_vx(idxrange_vxvyvz_x1Dsplit);
    PreallocatableSplineInterpolator const spline_exact_vx_interpolator(builder_exact_vx, spline_exact_vx_evaluator, idxrange_vxvyvz_x1Dsplit);

    SplineVyEvaluator const spline_exact_vy_evaluator(bv_vy_min, bv_vy_max);
    SplineVyBuilder const builder_exact_vy(idxrange_vxvyvz_x1Dsplit);
    PreallocatableSplineInterpolator const spline_exact_vy_interpolator(builder_exact_vy, spline_exact_vy_evaluator, idxrange_vxvyvz_x1Dsplit);

    SplineVzEvaluator const spline_exact_vz_evaluator(bv_vz_min, bv_vz_max);
    SplineVzBuilder const builder_exact_vz(idxrange_vxvyvz_x1Dsplit);
    PreallocatableSplineInterpolator const spline_exact_vz_interpolator(builder_exact_vz, spline_exact_vz_evaluator, idxrange_vxvyvz_x1Dsplit);


    // Create advection operator
    BslAdvectionSpatial<GeometryVxVyVzX, GridX> const advection_x(spline_x_interpolator);
    
    BslAdvectionVelocity<GeometryXVxVyVz, GridVx> const advection_vx(spline_vx_interpolator);
    BslAdvectionVelocity<GeometryXVxVyVz, GridVy> const advection_vy(spline_vy_interpolator);
    BslAdvectionVelocity<GeometryXVxVyVz, GridVz> const advection_vz(spline_vz_interpolator);

    //BslAdvectionVelocityRot3DVx<GeometryXVxVyVz, GridVx, GridVy, GridVz> const advec_3d_rot_vx(spline_vx_interpolator);
    //BslAdvectionVelocityRot3DVy<GeometryXVxVyVz, GridVx, GridVy, GridVz> const advec_3d_rot_vy(spline_vy_interpolator);
    //BslAdvectionVelocityRot3DVz<GeometryXVxVyVz, GridVx, GridVy, GridVz> const advec_3d_rot_vz(spline_vz_interpolator);

    BslAdvectionVelocityRot3DExact<GeometryXVxVyVz, GridX, GridVx, GridVy, GridVz> const advec_3d_rot(spline_exact_vx_interpolator, spline_exact_vy_interpolator, spline_exact_vz_interpolator);

    MpiSplitHybridVlasovSolver const
            vlasov(advection_x, advection_vx, advection_vy, advection_vz, advec_3d_rot, transpose);
    
    // 3D velocity quafrature coefficients 
    DFieldMemVxVyVz const quadrature_coeffs(
            neumann_spline_quadrature_coefficients<
                    Kokkos::DefaultExecutionSpace>(idxrange_vxvyvz, builder_vx, builder_vy, builder_vz));
    DFieldMemVxVyVz local_quadrature_coeffs(idxrange_vxvyvz_v3Dsplit);
    ddc::parallel_deepcopy(
            get_field(local_quadrature_coeffs),
            quadrature_coeffs[idxrange_vxvyvz_v3Dsplit]);

    // 1D (Vz) velocity quafrature coefficients 
    DFieldMemVz const quadrature_coeffs_Vz(
            neumann_spline_quadrature_coefficients<
                    Kokkos::DefaultExecutionSpace>(idxrange_vz, builder_vz));
    DFieldMemVz local_quadrature_coeffs_Vz(idxrange_vz_v3Dsplit);
    ddc::parallel_deepcopy(
            get_field(local_quadrature_coeffs_Vz),
            quadrature_coeffs_Vz[idxrange_vz_v3Dsplit]);
    
    FFTHybridSolver1D<IdxRangeX> fft_hybrid_solver(idxrange_x);
    MomentsCalculator const rhs_local(get_const_field(local_quadrature_coeffs), get_const_field(local_quadrature_coeffs_Vz));
    MpiMomentsCalculator const rhs(MPI_COMM_WORLD, rhs_local);
    HybridFieldSolver const hybrid_field(fft_hybrid_solver, rhs);
    
    // Create time splitting operator
    Hybridsplitting const hybridsplitting(vlasov, hybrid_field, rhs);
    
    // Starting the code
    ddc::expose_to_pdi("Nx_spline_cells", ddc::discrete_space<BSplinesX>().ncells());
    
    ddc::expose_to_pdi("Nvx_spline_cells", ddc::discrete_space<BSplinesVx>().ncells());
    ddc::expose_to_pdi("Nvy_spline_cells", ddc::discrete_space<BSplinesVy>().ncells());
    ddc::expose_to_pdi("Nvz_spline_cells", ddc::discrete_space<BSplinesVz>().ncells());
    
    expose_mesh_to_pdi("MeshX", idxrange_x);
    expose_mesh_to_pdi("MeshVx", idxrange_vx);
    expose_mesh_to_pdi("MeshVy", idxrange_vy);
    expose_mesh_to_pdi("MeshVz", idxrange_vz);

    ddc::expose_to_pdi("nbstep_diag", nbstep_diag);
    ddc::expose_to_pdi("Nkinspecies", idx_range_kinsp.size());
    ddc::expose_to_pdi(
            "fdistribu_charges",
            ddc::discrete_space<Species>().charges()[idx_range_kinsp]);
    ddc::expose_to_pdi(
            "fdistribu_masses",
            ddc::discrete_space<Species>().masses()[idx_range_kinsp]);
    if (rank == 0) {
        auto allfequilibrium_host = ddc::create_mirror_view_and_copy(get_field(allfequilibrium));
        ddc::PdiEvent("initial_state").with("fdistribu_eq", allfequilibrium_host);
    }

    steady_clock::time_point const start = steady_clock::now();

    transpose(
            Kokkos::DefaultExecutionSpace(),
            get_field(allfdistribu_v3D_split),
            get_const_field(allfdistribu_x1D_split));

    // Save the output index range
    IdxRangeSpXVxVyVz idxrange_spxyvxvy_v3Dsplit(idxrange_spvxvyvzx_v3Dsplit);
    PDI_expose_idx_range(idxrange_spxyvxvy_v3Dsplit, "local_fdistribu");

    hybridsplitting(get_field(allfdistribu_v3D_split), get_field(mean_current_x_each), get_field(mean_current_y_each), get_field(mean_current_z_each), 
                   get_field(mean_current_x), get_field(mean_current_y), get_field(mean_current_z), get_field(momentum_x), get_field(momentum_y), get_field(momentum_z), 
                   get_field(magnetic_field_x), get_field(magnetic_field_y), get_field(magnetic_field_z), get_field(pressure),
                   get_field(rho_each), get_field(rho), get_field(kinetic), deltat, nbiter);
    
    steady_clock::time_point const end = steady_clock::now();

    double const simulation_time = std::chrono::duration<double>(end - start).count();
    std::cout << "Simulation time: " << simulation_time << "s\n";

    PC_tree_destroy(&conf_pdi);

    PDI_finalize();

    MPI_Finalize();

    PC_tree_destroy(&conf_gyselalibxx);

    return EXIT_SUCCESS;
}
