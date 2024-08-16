// SPDX-License-Identifier: MIT
#include <cmath>

#include <ddc/ddc.hpp>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <pdi.h>

#include "constantfluidinitialization.hpp"
#include "ddc_alias_inline_functions.hpp"
#include "geometry.hpp"
#include "species_info.hpp"

/**
 * This test initializes a discrete space for moments (density,
 * particle_flux, stress) and initializes a neutral species
 * defined on the moment space and a spatial dimension with constant values
 * for the density, particle flux and stress. The test checks if the
 * initialization works properly.
 */
TEST(GeometryXM, MomentsInitialization)
{
    CoordX const x_min(0.0);
    CoordX const x_max(1.0);
    IdxStepX const x_size(50);

    // Creating mesh & supports
    ddc::init_discrete_space<BSplinesX>(x_min, x_max, x_size);
    ddc::init_discrete_space<GridX>(SplineInterpPointsX::get_sampling<GridX>());

    IdxRangeX meshX(SplineInterpPointsX::get_domain<GridX>());
    SplineXBuilder_1d const builder_x(meshX);

    // Kinetic species index range initialization
    IdxStepSp const nb_kinspecies(2);
    IdxRangeSp const idx_range_kinsp(IdxSp(0), nb_kinspecies);

    IdxSp const iion = idx_range_kinsp.front();
    IdxSp const ielec = idx_range_kinsp.back();

    host_t<DFieldMemSp> kinetic_charges(idx_range_kinsp);
    kinetic_charges(ielec) = -1.;
    kinetic_charges(iion) = 1.;

    host_t<DFieldMemSp> kinetic_masses(idx_range_kinsp);
    double const mass_ion(400.), mass_elec(1.);
    kinetic_masses(ielec) = mass_elec;
    kinetic_masses(iion) = mass_ion;

    // Neutral species index range initialization
    IdxStepSp const nb_fluidspecies(1);
    IdxRangeSp const idx_range_fluidsp(IdxSp(idx_range_kinsp.back() + 1), nb_fluidspecies);
    IdxSp const ifluid = idx_range_fluidsp.front();

    // neutrals charge is zero
    host_t<DFieldMemSp> fluid_charges(idx_range_fluidsp);
    ddc::parallel_fill(fluid_charges, 0.);

    host_t<DFieldMemSp> fluid_masses(idx_range_fluidsp);
    fluid_masses(ifluid) = kinetic_masses(iion);

    // Create the index range of kinetic species + fluid species
    IdxRangeSp const idx_range_allsp(IdxSp(0), nb_kinspecies + nb_fluidspecies);

    // Create a Field that contains charges of all species
    host_t<DFieldMemSp> charges(idx_range_allsp);

    // fill the Field with charges of kinetic species
    for (IdxSp isp : idx_range_kinsp) {
        charges(isp) = kinetic_charges(isp);
    }

    // fill the Field with charges of fluid species
    for (IdxSp isp : idx_range_fluidsp) {
        charges(isp) = fluid_charges(isp);
    }

    // Create a Field that contains masses of kinetic and fluid species
    host_t<DFieldMemSp> masses(idx_range_allsp);

    // fill the Field with masses of kinetic species
    for (IdxSp isp : idx_range_kinsp) {
        masses(isp) = kinetic_masses(isp);
    }

    // fill the Field with masses of fluid species
    for (IdxSp isp : idx_range_fluidsp) {
        masses(isp) = fluid_masses(isp);
    }

    ddc::init_discrete_space<Species>(std::move(charges), std::move(masses));

    // Moments index range initialization
    IdxStepMom const nb_fluid_moments(3);
    IdxRangeMom const meshM(IdxMom(0), nb_fluid_moments);
    ddc::init_discrete_space<GridMom>();

    IdxMom idensity(0);
    IdxMom iparticle_flux(1);
    IdxMom istress(2);

    // Neutral species initialization
    DFieldMemSpMomX neutrals_alloc(IdxRangeSpMomX(idx_range_fluidsp, meshM, meshX));
    DFieldSpMomX neutrals = get_field(neutrals_alloc);

    double const fluid_density_init(1.);
    double const fluid_particle_flux_init(0.5);
    double const fluid_stress_init(0.9);

    host_t<DFieldMemSpMom> moments_init_alloc(IdxRangeSpMom(idx_range_fluidsp, meshM));
    host_t<DFieldSpMom> moments_init = get_field(moments_init_alloc);

    moments_init(ifluid, idensity) = fluid_density_init;
    moments_init(ifluid, iparticle_flux) = fluid_particle_flux_init;
    moments_init(ifluid, istress) = fluid_stress_init;

    ConstantFluidInitialization fluid_init(moments_init);
    fluid_init(neutrals);

    auto neutrals_host = ddc::create_mirror_view_and_copy(neutrals);

    double const tolerance(1.e-12);
    ddc::for_each(get_idx_range<Species, GridX>(neutrals_host), [&](IdxSpX const ispx) {
        IdxSp const isp(ddc::select<Species>(ispx));
        IdxX const ix(ddc::select<GridX>(ispx));

        // test for equality of density
        IdxSpMomX const idensity_loc(isp, idensity, ix);
        EXPECT_LE(std::fabs(neutrals_host(idensity_loc) - fluid_density_init), tolerance);

        // test for equality of particle_flux
        IdxSpMomX const iparticle_flux_loc(isp, iparticle_flux, ix);
        EXPECT_LE(
                std::fabs(neutrals_host(iparticle_flux_loc) - fluid_particle_flux_init),
                tolerance);

        // test for equality of stresses
        IdxSpMomX const istress_loc(isp, istress, ix);
        EXPECT_LE(std::fabs(neutrals_host(istress_loc) - fluid_stress_init), tolerance);
    });
}
