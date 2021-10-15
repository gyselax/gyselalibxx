#include <chrono>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <iostream>

#include <ddc/MCoord>
#include <ddc/ProductMDomain>
#include <ddc/pdi.hpp>

#include <sll/null_boundary_value.hpp>
#include <sll/spline_evaluator.hpp>

#include <paraconf.h>
#include <pdi.h>

#include "fftpoissonsolver.hpp"
#include "fftw.hpp"
#include "geometry.hpp"
#include "ifftw.hpp"
#include "pdi_out.yml.hpp"
#include "predcorr.hpp"
#include "singlemodeperturbinitialization.hpp"
#include "species_info.hpp"
#include "splineadvectionvx.hpp"
#include "splineadvectionx.hpp"
#include "splitvlasovsolver.hpp"

using std::cerr;
using std::endl;
using std::chrono::steady_clock;
namespace fs = std::filesystem;

int main(int argc, char** argv)
{
    PC_tree_t conf_voicexx;
    if (argc > 1) {
        conf_voicexx = PC_parse_path(fs::path(argv[1]).c_str());
    } else {
        cerr << "usage: " << argv[0] << " <config_file.yml>" << endl;
        return EXIT_FAILURE;
    }

    // Reading config
    // --> Mesh info
    RCoordX x_min = [&]() {
        double x_min;
        PC_double(PC_get(conf_voicexx, ".Mesh.x_min"), &x_min);
        return RCoordX(x_min);
    }();
    RCoordX x_max = [&]() {
        double x_max;
        PC_double(PC_get(conf_voicexx, ".Mesh.x_max"), &x_max);
        return RCoordX(x_max);
    }();
    MLengthX x_size = [&]() {
        long x_size;
        PC_int(PC_get(conf_voicexx, ".Mesh.x_size"), &x_size);
        return MLengthX(x_size);
    }();
    RCoordVx vx_min = [&]() {
        double vx_min;
        PC_double(PC_get(conf_voicexx, ".Mesh.vx_min"), &vx_min);
        return RCoordVx(vx_min);
    }();
    RCoordVx vx_max = [&]() {
        double vx_max;
        PC_double(PC_get(conf_voicexx, ".Mesh.vx_max"), &vx_max);
        return RCoordVx(vx_max);
    }();
    MLengthVx vx_size = [&]() {
        double vx_size;
        PC_double(PC_get(conf_voicexx, ".Mesh.vx_size"), &vx_size);
        return MLengthVx(vx_size);
    }();

    // Creating mesh & supports
    BSplinesX const bsplines_x(x_min, x_max, x_size);

    SplineXBuilder const builder_x(bsplines_x);

    BSplinesVx const bsplines_vx(vx_min, vx_max, vx_size);

    SplineVxBuilder const builder_vx(bsplines_vx);

    MeshSp species;

    MDomainSp dom_sp(species, MLengthSp(2));

    MDomainSpXVx const
            mesh(dom_sp, builder_x.interpolation_domain(), builder_vx.interpolation_domain());

    BlockSp<int> charges(dom_sp);
    DBlockSp masses(dom_sp);
    DBlockSp density_eq(dom_sp);
    DBlockSp temperature_eq(dom_sp);
    DBlockSp mean_velocity_eq(dom_sp);
    BlockSp<int> init_perturb_mode(dom_sp);
    DBlockSp init_perturb_amplitude(dom_sp);
    for (MCoordSp isp : dom_sp) {
        // --> SpeciesInfo info
        long charge;
        PC_int(PC_get(conf_voicexx, ".SpeciesInfo.charge[%d]", isp.value()), &charge);
        charges(isp) = charge;
        PC_double(PC_get(conf_voicexx, ".SpeciesInfo.mass[%d]", isp.value()), &masses(isp));
        PC_double(
                PC_get(conf_voicexx, ".SpeciesInfo.density_eq[%d]", isp.value()),
                &density_eq(isp));
        PC_double(
                PC_get(conf_voicexx, ".SpeciesInfo.temperature_eq[%d]", isp.value()),
                &temperature_eq(isp));
        PC_double(
                PC_get(conf_voicexx, ".SpeciesInfo.mean_velocity_eq[%d]", isp.value()),
                &mean_velocity_eq(isp));

        // --> Perturbation info
        PC_double(
                PC_get(conf_voicexx, ".Perturbation.amplitude[%d]", isp.value()),
                &init_perturb_amplitude(isp));
        long init_perturb_mode_sp;
        PC_int(PC_get(conf_voicexx, ".Perturbation.mode[%d]", isp.value()), &init_perturb_mode_sp);
        init_perturb_mode(isp) = init_perturb_mode_sp;
    }

    // Initialization of the distribution function
    SpeciesInformation species_info(
            std::move(charges),
            std::move(masses),
            std::move(density_eq),
            std::move(temperature_eq),
            std::move(mean_velocity_eq),
            mesh);

    DBlockSpXVx allfdistribu(mesh.restrict(MDomainSp(species, species_info.ielec(), MLengthSp(1))));

    SingleModePerturbInitialization init(species_info, init_perturb_mode, init_perturb_amplitude);
    init(allfdistribu);

    // --> Algorithm info
    double deltat;
    PC_double(PC_get(conf_voicexx, ".Algorithm.deltat"), &deltat);
    long nbiter;
    PC_int(PC_get(conf_voicexx, ".Algorithm.nbiter"), &nbiter);

    // --> Output info
    double time_diag;
    PC_double(PC_get(conf_voicexx, ".Output.time_diag"), &time_diag);
    int const nbstep_diag = int(time_diag / deltat);

    PC_tree_t conf_pdi = PC_parse_string(PDI_CFG);

    PDI_init(conf_pdi);

    // Creating operators

    SplineEvaluator<BSplinesVx>
            spline_vx_evaluator(bsplines_vx, NullBoundaryValue::value, NullBoundaryValue::value);

    SplineEvaluator<BSplinesX>
            spline_x_evaluator(bsplines_x, NullBoundaryValue::value, NullBoundaryValue::value);

    SplineAdvectionX const advection_x(species_info, builder_x, spline_x_evaluator);

    SplineAdvectionVx const advection_vx(
            species_info,
            builder_x,
            spline_x_evaluator,
            builder_vx,
            spline_vx_evaluator);

    SplitVlasovSolver const vlasov(advection_x, advection_vx);

    FftwFourierTransform<Dim::X> fft;

    FftwInverseFourierTransform<Dim::X> ifft;

    FftPoissonSolver poisson(species_info, fft, ifft, builder_vx, spline_vx_evaluator);

    PredCorr const predcorr(vlasov, poisson, deltat);

    // Creating of mesh for output saving
    MDomainX gridx = select<MeshX>(mesh);
    BlockX<RCoordX> meshX_coord(gridx);
    for (MCoordX ix : gridx) {
        meshX_coord(ix) = gridx.to_real(ix);
    }

    MDomainVx gridvx = select<MeshVx>(mesh);
    BlockVx<RCoordVx> meshVx_coord(gridvx);
    for (MCoordVx ivx : gridvx) {
        meshVx_coord(ivx) = gridvx.to_real(ivx);
    }

    // Starting the code
    expose_to_pdi("Nx", x_size);
    expose_to_pdi("Nvx", vx_size);
    expose_to_pdi("MeshX", meshX_coord);
    expose_to_pdi("MeshVx", meshVx_coord);
    expose_to_pdi("nbstep_diag", nbstep_diag);
    PdiEvent("initial_state");

    steady_clock::time_point start = steady_clock::now();

    predcorr(allfdistribu, nbiter);

    steady_clock::time_point end = steady_clock::now();

    double simulation_time = std::chrono::duration<double>(end - start).count();
    std::cout << "Simulation time: " << simulation_time << "s\n";

    PC_tree_destroy(&conf_pdi);

    PDI_finalize();

    PC_tree_destroy(&conf_voicexx);

    return EXIT_SUCCESS;
}
