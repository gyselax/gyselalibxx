#include <chrono>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <iostream>

#include <ddc/DiscreteCoordinate>
#include <ddc/DiscreteDomain>
#include <ddc/pdi.hpp>

#include <sll/null_boundary_value.hpp>
#include <sll/spline_evaluator.hpp>

#include <paraconf.h>
#include <pdi.h>

#include "bsl_advection_vx.hpp"
#include "bsl_advection_x.hpp"
#include "fftpoissonsolver.hpp"
#include "fftw.hpp"
#include "geometry.hpp"
#include "ifftw.hpp"
#include "pdi_out.yml.hpp"
#include "predcorr.hpp"
#include "singlemodeperturbinitialization.hpp"
#include "species_info.hpp"
#include "spline_interpolator_vx.hpp"
#include "spline_interpolator_x.hpp"
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
    CoordX const x_min = [&]() {
        double x_min;
        PC_double(PC_get(conf_voicexx, ".Mesh.x_min"), &x_min);
        return CoordX(x_min);
    }();
    CoordX const x_max = [&]() {
        double x_max;
        PC_double(PC_get(conf_voicexx, ".Mesh.x_max"), &x_max);
        return CoordX(x_max);
    }();
    IVectX const x_size = [&]() {
        long x_size;
        PC_int(PC_get(conf_voicexx, ".Mesh.x_size"), &x_size);
        return IVectX(x_size);
    }();
    CoordVx const vx_min = [&]() {
        double vx_min;
        PC_double(PC_get(conf_voicexx, ".Mesh.vx_min"), &vx_min);
        return CoordVx(vx_min);
    }();
    CoordVx const vx_max = [&]() {
        double vx_max;
        PC_double(PC_get(conf_voicexx, ".Mesh.vx_max"), &vx_max);
        return CoordVx(vx_max);
    }();
    IVectVx const vx_size = [&]() {
        double vx_size;
        PC_double(PC_get(conf_voicexx, ".Mesh.vx_size"), &vx_size);
        return IVectVx(vx_size);
    }();

    // Creating mesh & supports
    init_discretization<BSplinesX>(x_min, x_max, x_size);

    init_discretization<BSplinesVx>(vx_min, vx_max, vx_size);

    SplineXBuilder const builder_x;

    SplineVxBuilder const builder_vx;

    IDomainSp const dom_sp(IVectSp(2));

    IDomainSpXVx const
            mesh(dom_sp, builder_x.interpolation_domain(), builder_vx.interpolation_domain());

    FieldSp<int> charges(dom_sp);
    DFieldSp masses(dom_sp);
    DFieldSp density_eq(dom_sp);
    DFieldSp temperature_eq(dom_sp);
    DFieldSp mean_velocity_eq(dom_sp);
    FieldSp<int> init_perturb_mode(dom_sp);
    DFieldSp init_perturb_amplitude(dom_sp);
    for (IndexSp const isp : dom_sp) {
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
    SpeciesInformation const species_info(
            std::move(charges),
            std::move(masses),
            std::move(density_eq),
            std::move(temperature_eq),
            std::move(mean_velocity_eq),
            mesh);

    DFieldSpXVx allfdistribu(mesh.restrict(IDomainSp(species_info.ielec(), IVectSp(1))));

    SingleModePerturbInitialization const
            init(species_info, init_perturb_mode, init_perturb_amplitude);
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

    SplineEvaluator<BSplinesVx> const
            spline_vx_evaluator(NullBoundaryValue::value, NullBoundaryValue::value);

    SplineEvaluator<BSplinesX> const
            spline_x_evaluator(NullBoundaryValue::value, NullBoundaryValue::value);

    PreallocatableSplineInterpolatorX const spline_x_interpolator(builder_x, spline_x_evaluator);

    PreallocatableSplineInterpolatorVx const
            spline_vx_interpolator(builder_vx, spline_vx_evaluator);

    BslAdvectionX const advection_x(species_info, spline_x_interpolator);

    BslAdvectionVx const
            advection_vx(species_info, builder_x, spline_x_evaluator, spline_vx_interpolator);

    SplitVlasovSolver const vlasov(advection_x, advection_vx);

    FftwFourierTransform<RDimX> const fft;

    FftwInverseFourierTransform<RDimX> const ifft;

    FftPoissonSolver const poisson(species_info, fft, ifft, builder_vx, spline_vx_evaluator);

    PredCorr const predcorr(vlasov, poisson, deltat);

    // Creating of mesh for output saving
    IDomainX const gridx = select<IDimX>(mesh);
    FieldX<CoordX> meshX_coord(gridx);
    for (IndexX const ix : gridx) {
        meshX_coord(ix) = to_real(ix);
    }

    IDomainVx const gridvx = select<IDimVx>(mesh);
    FieldVx<CoordVx> meshVx_coord(gridvx);
    for (IndexVx const ivx : gridvx) {
        meshVx_coord(ivx) = to_real(ivx);
    }

    // Starting the code
    expose_to_pdi("Nx", x_size.value());
    expose_to_pdi("Nvx", vx_size.value());
    expose_to_pdi("MeshX", meshX_coord);
    expose_to_pdi("MeshVx", meshVx_coord);
    expose_to_pdi("nbstep_diag", nbstep_diag);
    PdiEvent("initial_state");

    steady_clock::time_point const start = steady_clock::now();

    predcorr(allfdistribu, nbiter);

    steady_clock::time_point const end = steady_clock::now();

    double const simulation_time = std::chrono::duration<double>(end - start).count();
    std::cout << "Simulation time: " << simulation_time << "s\n";

    PC_tree_destroy(&conf_pdi);

    PDI_finalize();

    PC_tree_destroy(&conf_voicexx);

    return EXIT_SUCCESS;
}
