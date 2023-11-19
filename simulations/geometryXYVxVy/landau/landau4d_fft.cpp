// SPDX-License-Identifier: MIT

#include <chrono>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string_view>

#include <ddc/ddc.hpp>

#include <sll/constant_extrapolation_boundary_value.hpp>
#include <sll/null_boundary_value.hpp>
#include <sll/spline_evaluator.hpp>

#include <geometry.hpp>
#include <paraconf.h>
#include <pdi.h>

#include "bsl_advection_vx.hpp"
#include "bsl_advection_x.hpp"
#include "fftpoissonsolver.hpp"
#include "maxwellianequilibrium.hpp"
#include "paraconfpp.hpp"
#include "params.yaml.hpp"
#include "pdi_out.yml.hpp"
#include "predcorr.hpp"
#include "singlemodeperturbinitialization.hpp"
//#include "species_info.hpp"
#include "spline_interpolator.hpp"
#include "splitvlasovsolver.hpp"

using std::cerr;
using std::endl;
using std::chrono::steady_clock;
namespace fs = std::filesystem;

using PreallocatableSplineInterpolatorX
        = PreallocatableSplineInterpolator<IDimX, BSplinesX, SplineXBoundary, SplineXBoundary>;
using PreallocatableSplineInterpolatorY
        = PreallocatableSplineInterpolator<IDimY, BSplinesY, SplineYBoundary, SplineYBoundary>;
using PreallocatableSplineInterpolatorVx = PreallocatableSplineInterpolator<
        IDimVx,
        BSplinesVx,
        BoundCond::HERMITE,
        BoundCond::HERMITE>;
using PreallocatableSplineInterpolatorVy = PreallocatableSplineInterpolator<
        IDimVy,
        BSplinesVy,
        BoundCond::HERMITE,
        BoundCond::HERMITE>;
using BslAdvectionX = BslAdvectionSpatial<GeometryXYVxVy, IDimX>;
using BslAdvectionY = BslAdvectionSpatial<GeometryXYVxVy, IDimY>;
using BslAdvectionVx = BslAdvectionVelocity<GeometryXYVxVy, IDimVx>;
using BslAdvectionVy = BslAdvectionVelocity<GeometryXYVxVy, IDimVy>;

int main(int argc, char** argv)
{
    ddc::ScopeGuard scope(argc, argv);

    PC_tree_t conf_voicexx;
    if (argc == 2) {
        conf_voicexx = PC_parse_path(fs::path(argv[1]).c_str());
    } else if (argc == 3) {
        if (argv[1] == std::string_view("--dump-config")) {
            std::fstream file(argv[2], std::fstream::out);
            file << params_yaml;
            return EXIT_SUCCESS;
        }
    } else {
        cerr << "usage: " << argv[0] << " [--dump-config] <config_file.yml>" << endl;
        return EXIT_FAILURE;
    }
    PC_tree_t conf_pdi = PC_parse_string(PDI_CFG);
    PC_errhandler(PC_NULL_HANDLER);
    PDI_init(conf_pdi);

    // Reading config
    // --> Mesh info
    CoordX const x_min(PCpp_double(conf_voicexx, ".Mesh.x_min"));
    CoordX const x_max(PCpp_double(conf_voicexx, ".Mesh.x_max"));
    IVectX const x_size(PCpp_int(conf_voicexx, ".Mesh.x_size"));

    CoordY const y_min(PCpp_double(conf_voicexx, ".Mesh.y_min"));
    CoordY const y_max(PCpp_double(conf_voicexx, ".Mesh.y_max"));
    IVectY const y_size(PCpp_int(conf_voicexx, ".Mesh.y_size"));

    CoordVx const vx_min(PCpp_double(conf_voicexx, ".Mesh.vx_min"));
    CoordVx const vx_max(PCpp_double(conf_voicexx, ".Mesh.vx_max"));
    IVectVx const vx_size(PCpp_int(conf_voicexx, ".Mesh.vx_size"));

    CoordVy const vy_min(PCpp_double(conf_voicexx, ".Mesh.vy_min"));
    CoordVy const vy_max(PCpp_double(conf_voicexx, ".Mesh.vy_max"));
    IVectVy const vy_size(PCpp_int(conf_voicexx, ".Mesh.vy_size"));

    // Creating mesh & supports
    ddc::init_discrete_space<BSplinesX>(x_min, x_max, x_size);
    ddc::init_discrete_space<IDimX>(SplineInterpPointsX::get_sampling());
    ddc::DiscreteDomain<IDimX> interpolation_domain_x(SplineInterpPointsX::get_domain());
    SplineXBuilder const builder_x(interpolation_domain_x);

    ddc::init_discrete_space<BSplinesY>(y_min, y_max, y_size);
    ddc::init_discrete_space<IDimY>(SplineInterpPointsY::get_sampling());
    ddc::DiscreteDomain<IDimY> interpolation_domain_y(SplineInterpPointsY::get_domain());
    SplineYBuilder const builder_y(interpolation_domain_y);

    ddc::DiscreteDomain<IDimX, IDimY>
            interpolation_domain_xy(interpolation_domain_x, interpolation_domain_y);
    SplineXYBuilder const builder_xy(interpolation_domain_xy);

    ddc::init_discrete_space<BSplinesVx>(vx_min, vx_max, vx_size);
    ddc::init_discrete_space<IDimVx>(SplineInterpPointsVx::get_sampling());
    ddc::DiscreteDomain<IDimVx> interpolation_domain_vx(SplineInterpPointsVx::get_domain());
    SplineVxBuilder const builder_vx(interpolation_domain_vx);

    ddc::init_discrete_space<BSplinesVy>(vy_min, vy_max, vy_size);
    ddc::init_discrete_space<IDimVy>(SplineInterpPointsVy::get_sampling());
    ddc::DiscreteDomain<IDimVy> interpolation_domain_vy(SplineInterpPointsVy::get_domain());
    SplineVyBuilder const builder_vy(interpolation_domain_vy);

    ddc::DiscreteDomain<IDimVx, IDimVy>
            interpolation_domain_vxvy(interpolation_domain_vx, interpolation_domain_vy);
    SplineVxVyBuilder const builder_vxvy(interpolation_domain_vxvy);

    IVectSp const nb_kinspecies(PCpp_len(conf_voicexx, ".SpeciesInfo"));
    IDomainSp const dom_kinsp(IndexSp(0), nb_kinspecies);

    IDomainSpXYVxVy const meshSpXYVxVy(
            dom_kinsp,
            builder_x.interpolation_domain(),
            builder_y.interpolation_domain(),
            builder_vx.interpolation_domain(),
            builder_vy.interpolation_domain());

    IDomainSpVxVy const meshSpVxVy(
            dom_kinsp,
            builder_vx.interpolation_domain(),
            builder_vy.interpolation_domain());

    FieldSp<int> kinetic_charges(dom_kinsp);
    DFieldSp masses(dom_kinsp);
    DFieldSp density_eq(dom_kinsp);
    DFieldSp temperature_eq(dom_kinsp);
    DFieldSp mean_velocity_eq(dom_kinsp);
    DFieldSp init_perturb_amplitude(dom_kinsp);
    FieldSp<int> init_perturb_mode(dom_kinsp);
    int nb_elec_adiabspecies = 1;
    int nb_ion_adiabspecies = 1;

    for (IndexSp const isp : dom_kinsp) {
        // --> SpeciesInfo info
        PC_tree_t const conf_isp = PCpp_get(conf_voicexx, ".SpeciesInfo[%d]", isp.uid());

        kinetic_charges(isp) = static_cast<int>(PCpp_int(conf_isp, ".charge"));
        if (kinetic_charges(isp) == -1) {
            nb_elec_adiabspecies = 0;
        } else {
            nb_ion_adiabspecies = 0;
        }

        masses(isp) = PCpp_double(conf_isp, ".mass");
        density_eq(isp) = PCpp_double(conf_isp, ".density_eq");
        temperature_eq(isp) = PCpp_double(conf_isp, ".temperature_eq");
        mean_velocity_eq(isp) = PCpp_double(conf_isp, ".mean_velocity_eq");
        init_perturb_amplitude(isp) = PCpp_double(conf_isp, ".perturb_amplitude");
        init_perturb_mode(isp) = PCpp_double(conf_isp, ".perturb_mode");
    }

    // Create the domain of all species including kinetic species + adiabatic species (if existing)
    IDomainSp const
            dom_allsp(IndexSp(0), nb_kinspecies + nb_elec_adiabspecies + nb_ion_adiabspecies);
    FieldSp<int> charges(dom_allsp);
    for (IndexSp isp : dom_kinsp) {
        charges(isp) = kinetic_charges(isp);
    }
    if (nb_elec_adiabspecies + nb_ion_adiabspecies > 0) {
        charges(dom_kinsp.back() + 1) = nb_ion_adiabspecies - nb_elec_adiabspecies;
    }

    // Initialization of the distribution function
    ddc::init_discrete_space<IDimSp>(
            std::move(charges),
            std::move(masses),
            std::move(init_perturb_amplitude),
            std::move(init_perturb_mode));
    DFieldSpVxVy allfequilibrium(meshSpVxVy);
    MaxwellianEquilibrium const init_fequilibrium(
            std::move(density_eq),
            std::move(temperature_eq),
            std::move(mean_velocity_eq));
    init_fequilibrium(allfequilibrium);
    DFieldSpXYVxVy allfdistribu(meshSpXYVxVy);
    SingleModePerturbInitialization const
            init(allfequilibrium,
                 ddc::host_discrete_space<IDimSp>().perturb_modes(),
                 ddc::host_discrete_space<IDimSp>().perturb_amplitudes());
    init(allfdistribu);

    // --> Algorithm info
    double const deltat = PCpp_double(conf_voicexx, ".Algorithm.deltat");
    int const nbiter = static_cast<int>(PCpp_int(conf_voicexx, ".Algorithm.nbiter"));

    // --> Output info
    double const time_diag = PCpp_double(conf_voicexx, ".Output.time_diag");
    int const nbstep_diag = int(time_diag / deltat);

    // Create spline evaluator
    ConstantExtrapolationBoundaryValue<BSplinesX> bv_x_min(x_min);
    ConstantExtrapolationBoundaryValue<BSplinesX> bv_x_max(x_max);
    SplineEvaluator<BSplinesX> const spline_x_evaluator(bv_x_min, bv_x_max);
    PreallocatableSplineInterpolatorX const spline_x_interpolator(builder_x, spline_x_evaluator);

    ConstantExtrapolationBoundaryValue<BSplinesY> bv_y_min(y_min);
    ConstantExtrapolationBoundaryValue<BSplinesY> bv_y_max(y_max);
    SplineEvaluator<BSplinesY> const spline_y_evaluator(bv_y_min, bv_y_max);
    PreallocatableSplineInterpolatorY const spline_y_interpolator(builder_y, spline_y_evaluator);

    ConstantExtrapolationBoundaryValue<BSplinesVx> bv_vx_min(vx_min);
    ConstantExtrapolationBoundaryValue<BSplinesVx> bv_vx_max(vx_max);
    SplineEvaluator<BSplinesVx> const spline_vx_evaluator(bv_vx_min, bv_vx_max);
    PreallocatableSplineInterpolatorVx const
            spline_vx_interpolator(builder_vx, spline_vx_evaluator);

    ConstantExtrapolationBoundaryValue<BSplinesVy> bv_vy_min(vy_min);
    ConstantExtrapolationBoundaryValue<BSplinesVy> bv_vy_max(vy_max);
    SplineEvaluator<BSplinesVy> const spline_vy_evaluator(bv_vy_min, bv_vy_max);
    PreallocatableSplineInterpolatorVy const
            spline_vy_interpolator(builder_vy, spline_vy_evaluator);

    SplineXYEvaluator const spline_xy_evaluator(
            g_null_boundary_2d<BSplinesX, BSplinesY>,
            g_null_boundary_2d<BSplinesX, BSplinesY>,
            g_null_boundary_2d<BSplinesX, BSplinesY>,
            g_null_boundary_2d<BSplinesX, BSplinesY>);

    SplineVxVyEvaluator const spline_vxvy_evaluator(
            g_null_boundary_2d<BSplinesVx, BSplinesVy>,
            g_null_boundary_2d<BSplinesVx, BSplinesVy>,
            g_null_boundary_2d<BSplinesVx, BSplinesVy>,
            g_null_boundary_2d<BSplinesVx, BSplinesVy>);

    // Create advection operator
    BslAdvectionX const advection_x(spline_x_interpolator);
    BslAdvectionY const advection_y(spline_y_interpolator);
    BslAdvectionVx const advection_vx(spline_vx_interpolator);
    BslAdvectionVy const advection_vy(spline_vy_interpolator);

    SplitVlasovSolver const vlasov(advection_x, advection_y, advection_vx, advection_vy);

    ddc::init_fourier_space<RDimX, RDimY>(ddc::select<IDimX, IDimY>(meshSpXYVxVy));

    FftPoissonSolver const
            poisson(builder_xy, spline_xy_evaluator, builder_vxvy, spline_vxvy_evaluator);

    // Create predcorr operator
    PredCorr const predcorr(vlasov, poisson);

    // Creating of mesh for output saving
    IDomainX const gridx = ddc::select<IDimX>(meshSpXYVxVy);
    FieldX<CoordX> meshX_coord(gridx);
    ddc::for_each(gridx, [&](IndexX const ix) { meshX_coord(ix) = ddc::coordinate(ix); });

    IDomainY const gridy = ddc::select<IDimY>(meshSpXYVxVy);
    FieldY<CoordY> meshY_coord(gridy);
    ddc::for_each(gridy, [&](IndexY const iy) { meshY_coord(iy) = ddc::coordinate(iy); });

    IDomainVx const gridvx = ddc::select<IDimVx>(meshSpVxVy);
    FieldVx<CoordVx> meshVx_coord(gridvx);
    for (IndexVx const ivx : gridvx) {
        meshVx_coord(ivx) = ddc::coordinate(ivx);
    }

    IDomainVy const gridvy = ddc::select<IDimVy>(meshSpVxVy);
    FieldVy<CoordVy> meshVy_coord(gridvy);
    for (IndexVy const ivy : gridvy) {
        meshVy_coord(ivy) = ddc::coordinate(ivy);
    }

    // Starting the code
    ddc::expose_to_pdi("Nx", x_size.value());
    ddc::expose_to_pdi("Ny", y_size.value());
    ddc::expose_to_pdi("Nvx", vx_size.value());
    ddc::expose_to_pdi("Nvy", vy_size.value());
    ddc::expose_to_pdi("MeshX", meshX_coord);
    ddc::expose_to_pdi("MeshY", meshY_coord);
    ddc::expose_to_pdi("MeshVx", meshVx_coord);
    ddc::expose_to_pdi("MeshVy", meshVy_coord);
    ddc::expose_to_pdi("nbstep_diag", nbstep_diag);
    ddc::expose_to_pdi("Nkinspecies", nb_kinspecies.value());
    ddc::expose_to_pdi("fdistribu_charges", ddc::discrete_space<IDimSp>().charges()[dom_kinsp]);
    ddc::expose_to_pdi("fdistribu_masses", ddc::discrete_space<IDimSp>().masses()[dom_kinsp]);
    ddc::PdiEvent("initial_state").with("fdistribu_eq", allfequilibrium);

    steady_clock::time_point const start = steady_clock::now();

    predcorr(allfdistribu, deltat, nbiter);

    steady_clock::time_point const end = steady_clock::now();

    double const simulation_time = std::chrono::duration<double>(end - start).count();
    std::cout << "Simulation time: " << simulation_time << "s\n";

    PC_tree_destroy(&conf_pdi);

    PDI_finalize();

    PC_tree_destroy(&conf_voicexx);

    return EXIT_SUCCESS;
}
