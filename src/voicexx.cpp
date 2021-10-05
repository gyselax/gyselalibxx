#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <iostream>

#include <ddc/MCoord>
#include <ddc/ProductMDomain>
#include <ddc/pdi.hpp>

#include <paraconf.h>
#include <pdi.h>

#include "fdistribu.h"
#include "geometry.h"
#include "nulladvectionvx.h"
#include "nullefieldsolver.h"
#include "pdi_out.yml.h"
#include "predcorr.h"
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
    // --> Mesh
    RCoordX x_min = [&]() {
        double x_min;
        PC_double(PC_get(conf_voicexx, ".Mesh.x_min"), &x_min);
        return x_min;
    }();
    RCoordX x_max = [&]() {
        double x_max;
        PC_double(PC_get(conf_voicexx, ".Mesh.x_max"), &x_max);
        return x_max;
    }();
    MLengthX x_size = [&]() {
        long x_size;
        PC_int(PC_get(conf_voicexx, ".Mesh.x_size"), &x_size);
        return x_size;
    }();
    RCoordVx vx_min = [&]() {
        double vx_min;
        PC_double(PC_get(conf_voicexx, ".Mesh.vx_min"), &vx_min);
        return vx_min;
    }();
    RCoordVx vx_max = [&]() {
        double vx_max;
        PC_double(PC_get(conf_voicexx, ".Mesh.vx_max"), &vx_max);
        return vx_max;
    }();
    MLengthVx vx_size = [&]() {
        double vx_size;
        PC_double(PC_get(conf_voicexx, ".Mesh.vx_size"), &vx_size);
        return vx_size;
    }();
    // --> Equilibrium
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

    double deltat;
    PC_double(PC_get(conf_voicexx, ".Algorithm.deltat"), &deltat);
    long nbiter;
    PC_int(PC_get(conf_voicexx, ".Algorithm.nbiter"), &nbiter);

    PC_tree_t conf_pdi = PC_parse_string(PDI_CFG);

    PDI_init(conf_pdi);

    // Creating mesh & supports

    BSplinesX const bsplines_x(x_min, x_max, x_size);

    SplineXBuilder const builder_x(bsplines_x);

    MeshX const mesh_x(x_min, x_max, x_size);

    MeshVx const mesh_vx(vx_min, vx_max, vx_size);


    MDomainXVx const
            dom2d(builder_x.interpolation_domain().mesh<MeshX>(),
                  mesh_vx,
                  MCoordXVx(builder_x.interpolation_domain().extents(), vx_size));

    // Creating operators

    SplineAdvectionX const advection_x(bsplines_x, builder_x);

    NullAdvectionVx const advection_vx;

    SplitVlasovSolver const vlasov(advection_x, advection_vx);

    NullEfieldSolver const efield;

    PredCorr const predcorr(vlasov, efield, deltat);

    // Creating of mesh for output saving
    MDomainX gridx = select<MeshX>(dom2d);
    BlockX<RCoordX> meshX_coord(gridx);
    for (MCoordX ix : gridx) {
        meshX_coord(ix) = gridx.to_real(ix);
        //cout << meshX_coord(ix);
    }
    //cout << endl;

    /*
    UniformMesh<Dim::Vx> const meshx_uniform(x_min, x_max, x_size);
    MDomainX const domainx_uniform(meshx_uniform);
    for (MCoordX ix1 : domainx_uniform()) {
        cout << domainx_uniform.to_real(ix1) ;
    }
    cout << endl;
    */

    MDomainVx gridvx = select<MeshVx>(dom2d);
    BlockVx<RCoordVx> meshVx_coord(gridvx);
    for (MCoordVx ivx : gridvx) {
        meshVx_coord(ivx) = gridvx.to_real(ivx);
    }

    // Initialization of the distribution function
    DistributionFunction
            fion(ion_charge,
                 ion_mass,
                 ion_density_eq,
                 ion_temperature_eq,
                 ion_mean_velocity_eq,
                 dom2d);
    cout << "ion charge=" << fion.charge << endl;
    fion.init();

    // Starting the code
    expose_to_pdi("Nx", x_size);
    expose_to_pdi("Nvx", vx_size);
    expose_to_pdi("MeshX", meshX_coord);
    expose_to_pdi("MeshVx", meshVx_coord);

    predcorr(fion, electron_mass, nbiter);

    PC_tree_destroy(&conf_pdi);

    PDI_finalize();

    PC_tree_destroy(&conf_voicexx);

    return EXIT_SUCCESS;
}
