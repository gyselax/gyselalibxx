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
#include "bsl_advection_2d_rot_vx.hpp"
#include "bsl_advection_2d_rot_vy.hpp"
#include "bsl_advection_x.hpp"
#include "momentscalculator.hpp"
#include "ddc_alias_inline_functions.hpp"
#include "fft_poisson_solver.hpp"
#include "fft_hybrid_solver.hpp"
#include "geometry.hpp"
#include "input.hpp"
#include "maxwellianequilibrium.hpp"
#include "mpimomentscalculator.hpp"
#include "mpisplithybridvlasovsolver.hpp"
#include "mpichargedensitycalculator.hpp"
#include "mpisplitvlasovsolver.hpp"
#include "mpitransposealltoall.hpp"
#include "neumann_spline_quadrature.hpp"
#include "output.hpp"
#include "paraconfpp.hpp"
#include "params.yaml.hpp"
#include "qnsolver.hpp"
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

    IdxRangeY const idxrange_y = init_spline_dependent_idx_range<
            GridY,
            BSplinesY,
            SplineInterpPointsY>(conf_gyselalibxx, "y");
    IdxRangeVx const idxrange_vx = init_spline_dependent_idx_range<
            GridVx,
            BSplinesVx,
            SplineInterpPointsVx>(conf_gyselalibxx, "vx");
    IdxRangeVy const idxrange_vy = init_spline_dependent_idx_range<
            GridVy,
            BSplinesVy,
            SplineInterpPointsVy>(conf_gyselalibxx, "vy");
    IdxRangeXY const idxrange_xy(idxrange_x, idxrange_y);
    IdxRangeVxVy idxrange_vxvy(idxrange_vx, idxrange_vy);
    IdxRangeXYVxVy const idxrange_xyvxvy(idxrange_x, idxrange_y, idxrange_vx, idxrange_vy);

    IdxRangeSp const idx_range_kinsp = init_species(conf_gyselalibxx);

    // set the index range of charge and mean velocity dnesity of each kind of ions.
    IdxRangeSpXY const idxrange_spxy(idx_range_kinsp, idxrange_x, idxrange_y);

    IdxRangeSpXYVxVy const idxrange_glob_spxyvxvy(idx_range_kinsp, idxrange_xyvxvy);

    MPITransposeAllToAll<X2DSplit, V2DSplit> transpose(idxrange_glob_spxyvxvy, MPI_COMM_WORLD);

    IdxRangeSpXYVxVy idxrange_spxyvxvy_x2Dsplit(transpose.get_local_idx_range<X2DSplit>());
    IdxRangeSpVxVyXY idxrange_spvxvyxy_v2Dsplit(transpose.get_local_idx_range<V2DSplit>());

    IdxRangeVxVy idxrange_vxvy_v2Dsplit(idxrange_spvxvyxy_v2Dsplit);

    IdxRangeVxVyXY idxrange_vxvyxy_v2Dsplit(idxrange_spvxvyxy_v2Dsplit);
    IdxRangeXYVxVy idxrange_xyvxvy_x2Dsplit(idxrange_spxyvxvy_x2Dsplit);
    SplineXBuilder const builder_x(idxrange_vxvyxy_v2Dsplit);
    SplineYBuilder const builder_y(idxrange_vxvyxy_v2Dsplit);
    SplineVxBuilder const builder_vx(idxrange_xyvxvy_x2Dsplit);
    SplineVyBuilder const builder_vy(idxrange_xyvxvy_x2Dsplit);

    IdxRangeSpVxVy idxrange_spvxvy_local(idxrange_spxyvxvy_x2Dsplit);
    // Initialisation of the distribution function
    DFieldMemSpVxVy allfequilibrium(idxrange_spvxvy_local);
    MaxwellianEquilibrium const init_fequilibrium
            = MaxwellianEquilibrium::init_from_input(idx_range_kinsp, conf_gyselalibxx);
    init_fequilibrium(get_field(allfequilibrium));
    DFieldMemSpXYVxVy allfdistribu_x2D_split(idxrange_spxyvxvy_x2Dsplit);
    DFieldMemSpVxVyXY allfdistribu_v2D_split(idxrange_spvxvyxy_v2Dsplit);
    SingleModePerturbInitialisation const init = SingleModePerturbInitialisation::
            init_from_input(get_const_field(allfequilibrium), idx_range_kinsp, conf_gyselalibxx);
    init(get_field(allfdistribu_x2D_split));

    // Define the date structure of the mean velocities and charge density for each kind of ions.
    DFieldMemSpXY mean_current_x_each(idxrange_spxy);
    DFieldMemSpXY mean_current_y_each(idxrange_spxy);
    DFieldMemSpXY rho_each(idxrange_spxy);
    // Define the date structure of the mean velocity and charge densities for all kinds of ions
    DFieldMemXY mean_current_x(idxrange_xy);
    DFieldMemXY mean_current_y(idxrange_xy);
    DFieldMemXY rho(idxrange_xy);
    // Define the kinetic density 
    DFieldMemXY kinetic(idxrange_xy);
    // Define the date structure of the magnetic field Bz and the pressure
    DFieldMemXY magnetic_field_z(idxrange_xy);
    DFieldMemXY pressure(idxrange_xy);
    int magnetic_mode = 1;
    double magnetic_amplitude = 0.1;
    int pressure_mode = 1;
    double pressure_amplitude = 0.1;
    Hybridmodel_field_initialisation const init_fields(magnetic_mode, magnetic_amplitude, pressure_mode, pressure_amplitude);
    init_fields(get_field(magnetic_field_z), get_field(pressure));
    // momentum 
    DFieldMemXY momentum_x(idxrange_xy);
    DFieldMemXY momentum_y(idxrange_xy);

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

    PreallocatableSplineInterpolator const spline_x_interpolator(builder_x, spline_x_evaluator, idxrange_vxvyxy_v2Dsplit);

    ddc::PeriodicExtrapolationRule<Y> bv_y_min;
    ddc::PeriodicExtrapolationRule<Y> bv_y_max;
    SplineYEvaluator const spline_y_evaluator(bv_y_min, bv_y_max);

    PreallocatableSplineInterpolator const spline_y_interpolator(builder_y, spline_y_evaluator, idxrange_vxvyxy_v2Dsplit);

    ddc::ConstantExtrapolationRule<Vx> bv_vx_min(ddc::coordinate(idxrange_vx.front()));
    ddc::ConstantExtrapolationRule<Vx> bv_vx_max(ddc::coordinate(idxrange_vx.back()));
    SplineVxEvaluator const spline_vx_evaluator(bv_vx_min, bv_vx_max);

    PreallocatableSplineInterpolator const spline_vx_interpolator(builder_vx, spline_vx_evaluator, idxrange_xyvxvy_x2Dsplit);

    ddc::ConstantExtrapolationRule<Vy> bv_vy_min(ddc::coordinate(idxrange_vy.front()));
    ddc::ConstantExtrapolationRule<Vy> bv_vy_max(ddc::coordinate(idxrange_vy.back()));
    SplineVyEvaluator const spline_vy_evaluator(bv_vy_min, bv_vy_max);

    PreallocatableSplineInterpolator const spline_vy_interpolator(builder_vy, spline_vy_evaluator, idxrange_xyvxvy_x2Dsplit);

    // Create advection operator
    BslAdvectionSpatial<GeometryVxVyXY, GridX> const advection_x(spline_x_interpolator);
    BslAdvectionSpatial<GeometryVxVyXY, GridY> const advection_y(spline_y_interpolator);
    BslAdvectionVelocity<GeometryXYVxVy, GridVx> const advection_vx(spline_vx_interpolator);
    BslAdvectionVelocity<GeometryXYVxVy, GridVy> const advection_vy(spline_vy_interpolator);
    BslAdvectionVelocityRot2DVx<GeometryXYVxVy, GridVx> const advection_2d_rot_vx(spline_vx_interpolator);
    BslAdvectionVelocityRot2DVy<GeometryXYVxVy, GridVy> const advection_2d_rot_vy(spline_vy_interpolator);

    MpiSplitHybridVlasovSolver const
            vlasov(advection_x, advection_y, advection_vx, advection_vy, advection_2d_rot_vx, advection_2d_rot_vy, transpose);

    DFieldMemVxVy const quadrature_coeffs(
            neumann_spline_quadrature_coefficients<
                    Kokkos::DefaultExecutionSpace>(idxrange_vxvy, builder_vx, builder_vy));
    DFieldMemVxVy local_quadrature_coeffs(idxrange_vxvy_v2Dsplit);
    ddc::parallel_deepcopy(
            get_field(local_quadrature_coeffs),
            quadrature_coeffs[idxrange_vxvy_v2Dsplit]);

    FFTHybridSolver<IdxRangeXY> fft_hybrid_solver(idxrange_xy);
    MomentsCalculator const rhs_local(get_const_field(local_quadrature_coeffs));
    MpiMomentsCalculator const rhs(MPI_COMM_WORLD, rhs_local);
    HybridFieldSolver const hybrid_field(fft_hybrid_solver, rhs);

    //FFTPoissonSolver<IdxRangeXY> fft_poisson_solver(idxrange_xy);
    //ChargeDensityCalculator const rhs_local_poisson(get_const_field(local_quadrature_coeffs));
    //MpiChargeDensityCalculator const rhs_poisson(MPI_COMM_WORLD, rhs_local_poisson);
    //QNSolver const poisson(fft_poisson_solver, rhs_poisson);

    // Create predcorr operator
    Hybridsplitting const hybridsplitting(vlasov, hybrid_field, rhs);

    // Starting the code
    ddc::expose_to_pdi("Nx_spline_cells", ddc::discrete_space<BSplinesX>().ncells());
    ddc::expose_to_pdi("Ny_spline_cells", ddc::discrete_space<BSplinesY>().ncells());
    ddc::expose_to_pdi("Nvx_spline_cells", ddc::discrete_space<BSplinesVx>().ncells());
    ddc::expose_to_pdi("Nvy_spline_cells", ddc::discrete_space<BSplinesVy>().ncells());
    expose_mesh_to_pdi("MeshX", idxrange_x);
    expose_mesh_to_pdi("MeshY", idxrange_y);
    expose_mesh_to_pdi("MeshVx", idxrange_vx);
    expose_mesh_to_pdi("MeshVy", idxrange_vy);
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
            get_field(allfdistribu_v2D_split),
            get_const_field(allfdistribu_x2D_split));

    // Save the output index range
    IdxRangeSpXYVxVy idxrange_spxyvxvy_v2Dsplit(idxrange_spvxvyxy_v2Dsplit);
    PDI_expose_idx_range(idxrange_spxyvxvy_v2Dsplit, "local_fdistribu");

    hybridsplitting(get_field(allfdistribu_v2D_split), get_field(mean_current_x_each), get_field(mean_current_y_each), 
                   get_field(mean_current_x), get_field(mean_current_y), get_field(momentum_x), get_field(momentum_y), 
                   get_field(magnetic_field_z), get_field(pressure),
                   get_field(rho_each), get_field(rho), get_field(kinetic), deltat, nbiter);

    //hybridsplitting(get_field(allfdistribu_v2D_split), deltat, nbiter);

    steady_clock::time_point const end = steady_clock::now();

    double const simulation_time = std::chrono::duration<double>(end - start).count();
    std::cout << "Simulation time: " << simulation_time << "s\n";

    PC_tree_destroy(&conf_pdi);

    PDI_finalize();

    MPI_Finalize();

    PC_tree_destroy(&conf_gyselalibxx);

    return EXIT_SUCCESS;
}
