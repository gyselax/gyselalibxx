#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <iostream>

#include <ddc/MCoord>
#include <ddc/ProductMDomain>
#include <ddc/pdi.hpp>

#include <paraconf.h>
#include <pdi.h>

#include "efieldfftsolver.h"
#include "fdistribu.h"
#include "fftw.h"
#include "geometry.h"
#include "ifftw.h"
#include "pdi_out.yml.h"
#include "predcorr.h"
#include "splineadvectionvx.h"
#include "splineadvectionx.h"
#include "splitvlasovsolver.h"

using std::cerr;
using std::cout;
using std::endl;
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

    // --> Equilibrium info
    long ion_charge;
    PC_int(PC_get(conf_voicexx, ".Equilibrium.ion_charge"), &ion_charge);
    double ion_mass;
    PC_double(PC_get(conf_voicexx, ".Equilibrium.ion_mass"), &ion_mass);
    double ion_density_eq;
    PC_double(PC_get(conf_voicexx, ".Equilibrium.ion_density_eq"), &ion_density_eq);
    double ion_temperature_eq;
    PC_double(PC_get(conf_voicexx, ".Equilibrium.ion_temperature_eq"), &ion_temperature_eq);
    double ion_mean_velocity_eq;
    PC_double(PC_get(conf_voicexx, ".Equilibrium.ion_mean_velocity_eq"), &ion_mean_velocity_eq);
    long electron_charge;
    PC_int(PC_get(conf_voicexx, ".Equilibrium.electron_charge"), &electron_charge);
    double electron_mass;
    PC_double(PC_get(conf_voicexx, ".Equilibrium.electron_mass"), &electron_mass);
    double electron_density_eq;
    PC_double(PC_get(conf_voicexx, ".Equilibrium.electron_density_eq"), &electron_density_eq);
    double electron_temperature_eq;
    PC_double(
            PC_get(conf_voicexx, ".Equilibrium.electron_temperature_eq"),
            &electron_temperature_eq);
    double electron_mean_velocity_eq;
    PC_double(
            PC_get(conf_voicexx, ".Equilibrium.electron_mean_velocity_eq"),
            &electron_mean_velocity_eq);

    // --> Perturbation info
    long init_perturb_mode;
    PC_int(PC_get(conf_voicexx, ".Perturbation.mode"), &init_perturb_mode);
    double init_perturb_amplitude;
    PC_double(PC_get(conf_voicexx, ".Perturbation.amplitude"), &init_perturb_amplitude);

    // --> Algorithm info
    double deltat;
    PC_double(PC_get(conf_voicexx, ".Algorithm.deltat"), &deltat);
    long nbiter;
    PC_int(PC_get(conf_voicexx, ".Algorithm.nbiter"), &nbiter);

    // --> Output info
    double time_diag;
    PC_double(PC_get(conf_voicexx, ".Output.time_diag"), &time_diag);
    int const nbstep_diag = int(time_diag / deltat);
    std::cout << "nbstep_diag= " << nbstep_diag << std::endl;

    PC_tree_t conf_pdi = PC_parse_string(PDI_CFG);

    PDI_init(conf_pdi);

    // Creating mesh & supports

    BSplinesX const bsplines_x(x_min, x_max, x_size);

    SplineXBuilder const builder_x(bsplines_x);

    BSplinesVx const bsplines_vx(vx_min, vx_max, vx_size);

    SplineVxBuilder const builder_vx(bsplines_vx);

    MDomainXVx const dom2d(builder_x.interpolation_domain(), builder_vx.interpolation_domain());

    // Creating operators

    SplineAdvectionX const advection_x(bsplines_x, builder_x);

    SplineAdvectionVx const advection_vx(bsplines_vx, builder_vx);

    SplitVlasovSolver const vlasov(advection_x, advection_vx);

    FftwFourierTransform<Dim::X> fft;

    FftwInverseFourierTransform<Dim::X> ifft;

    EfieldFftSolver efield(fft, ifft, bsplines_vx, builder_vx);

    PredCorr const predcorr(vlasov, efield, deltat, time_diag);

    // Creating of mesh for output saving
    MDomainX gridx = select<MeshX>(dom2d);
    BlockX<RCoordX> meshX_coord(gridx);
    for (MCoordX ix : gridx) {
        meshX_coord(ix) = gridx.to_real(ix);
    }

    MDomainVx gridvx = select<MeshVx>(dom2d);
    BlockVx<RCoordVx> meshVx_coord(gridvx);
    for (MCoordVx ivx : gridvx) {
        meshVx_coord(ivx) = gridvx.to_real(ivx);
    }

    // Initialization of the distribution function
    DistributionFunction felectron(
            electron_charge,
            electron_mass,
            electron_density_eq,
            electron_temperature_eq,
            electron_mean_velocity_eq,
            init_perturb_mode,
            init_perturb_amplitude,
            dom2d);
    felectron.init();

    // Starting the code
    expose_to_pdi("Nx", x_size);
    expose_to_pdi("Nvx", vx_size);
    expose_to_pdi("MeshX", meshX_coord);
    expose_to_pdi("MeshVx", meshVx_coord);
    expose_to_pdi("nbstep_diag", nbstep_diag);
    PdiEvent("initial_state");

    predcorr(felectron, electron_mass, nbiter);

    PC_tree_destroy(&conf_pdi);

    PDI_finalize();

    PC_tree_destroy(&conf_voicexx);

    return EXIT_SUCCESS;
}
