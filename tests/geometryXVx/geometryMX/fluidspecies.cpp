// SPDX-License-Identifier: MIT

#include <cmath>

#include <ddc/ddc.hpp>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <pdi.h>

#include "constantfluidinitialization.hpp"
#include "geometry.hpp"
#include "species_info.hpp"

/**
 * This test initializes a discrete space for kinetic species 
 * and fluid species, with corresponding masses and charges. 
 * The test checks if the masses() and charges() attributes 
 * of the discrete space correspond to the masses and charges
 * attributes of the fluid and kinetic species.
 */
TEST(GeometryXM, KineticFluidSpecies)
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
    IdxRangeSp const dom_kinsp(IdxSp(0), nb_kinspecies);

    IdxSp const my_iion = dom_kinsp.front();
    IdxSp const my_ielec = dom_kinsp.back();

    host_t<DFieldMemSp> kinetic_charges(dom_kinsp);
    kinetic_charges(my_ielec) = -1.;
    kinetic_charges(my_iion) = 1.;

    host_t<DFieldMemSp> kinetic_masses(dom_kinsp);
    double const mass_ion(400.), mass_elec(1.);
    kinetic_masses(my_ielec) = mass_elec;
    kinetic_masses(my_iion) = mass_ion;

    // Fluid species index range initialization
    IdxStepSp const nb_fluidspecies(2);
    IdxRangeSp const dom_fluidsp(IdxSp(dom_kinsp.back() + 1), nb_fluidspecies);

    host_t<DFieldMemSp> fluid_charges(dom_fluidsp);
    fluid_charges(dom_fluidsp.front()) = 1.;
    fluid_charges(dom_fluidsp.back()) = -1.;

    host_t<DFieldMemSp> fluid_masses(dom_fluidsp);
    fluid_masses(dom_fluidsp.front()) = 5.;
    fluid_masses(dom_fluidsp.back()) = 8.;

    // Create the index range of kinetic species + fluid species
    IdxRangeSp const dom_allsp(IdxSp(0), nb_kinspecies + nb_fluidspecies);

    // Create a Field that contains charges of all species
    host_t<DFieldMemSp> charges(dom_allsp);

    // fill the Field with charges of kinetic species
    for (IdxSp isp : dom_kinsp) {
        charges(isp) = kinetic_charges(isp);
    }

    // fill the Field with charges of fluid species
    for (IdxSp isp : dom_fluidsp) {
        charges(isp) = fluid_charges(isp);
    }

    // Create a Field that contains masses of kinetic and fluid species
    host_t<DFieldMemSp> masses(dom_allsp);

    // fill the Field with masses of kinetic species
    for (IdxSp isp : dom_kinsp) {
        masses(isp) = kinetic_masses(isp);
    }

    // fill the Field with masses of fluid species
    for (IdxSp isp : dom_fluidsp) {
        masses(isp) = fluid_masses(isp);
    }

    ddc::init_discrete_space<Species>(std::move(charges), std::move(masses));

    ddc::for_each(dom_allsp, [&](IdxSp const isp) {
        if (isp.uid() < nb_kinspecies) {
            EXPECT_EQ(ddc::discrete_space<Species>().charges()(isp), kinetic_charges(isp));
            EXPECT_EQ(ddc::discrete_space<Species>().masses()(isp), kinetic_masses(isp));

        } else if (nb_kinspecies + 1 <= isp.uid() && isp.uid() < nb_kinspecies + nb_fluidspecies) {
            EXPECT_EQ(ddc::discrete_space<Species>().charges()(isp), fluid_charges(isp));
            EXPECT_EQ(ddc::discrete_space<Species>().masses()(isp), fluid_masses(isp));
        }
    });
}

/**
 * Same test as above with an adiabatic species
 * on top of the fluid and kinetic species.
 */
TEST(GeometryXM, KineticFluidAdiabaticSpecies)
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
    IdxRangeSp const dom_kinsp(IdxSp(0), nb_kinspecies);

    IdxSp const my_iion = dom_kinsp.front();
    IdxSp const my_ielec = dom_kinsp.back();

    host_t<DFieldMemSp> kinetic_charges(dom_kinsp);
    kinetic_charges(my_ielec) = -1.;
    kinetic_charges(my_iion) = 1.;

    host_t<DFieldMemSp> kinetic_masses(dom_kinsp);
    double const mass_ion(400.), mass_elec(1.);
    kinetic_masses(my_ielec) = mass_elec;
    kinetic_masses(my_iion) = mass_ion;

    // adiabatic species initialization
    int nb_ion_adiabspecies = 1;

    // Fluid species index range initialization
    IdxStepSp const nb_fluidspecies(2);
    IdxRangeSp const dom_fluidsp(IdxSp(dom_kinsp.back() + 1), nb_fluidspecies);

    host_t<DFieldMemSp> fluid_charges(dom_fluidsp);
    fluid_charges(dom_fluidsp.front()) = 1.;
    fluid_charges(dom_fluidsp.back()) = -1.;

    host_t<DFieldMemSp> fluid_masses(dom_fluidsp);
    fluid_masses(dom_fluidsp.front()) = 5.;
    fluid_masses(dom_fluidsp.back()) = 8.;

    // Create the index range of all species including kinetic species + fluid species + adiabatic species (if existing)
    // adiabatic species are placed at the back of the index range
    IdxRangeSp const dom_allsp(IdxSp(0), nb_kinspecies + nb_fluidspecies + nb_ion_adiabspecies);

    // Create a Field that contains charges of all species
    host_t<DFieldMemSp> charges(dom_allsp);

    // fill the Field with charges of kinetic species
    for (IdxSp isp : dom_kinsp) {
        charges(isp) = kinetic_charges(isp);
    }

    // fill the Field with charges of fluid species
    for (IdxSp isp : dom_fluidsp) {
        charges(isp) = fluid_charges(isp);
    }

    // fill the Field with charges of adiabatic species
    double const charge_adiabspecies(3.);
    charges(dom_allsp.back()) = nb_ion_adiabspecies * charge_adiabspecies;

    // Create the index range of kinetic and fluid species
    IdxRangeSp const dom_kinfluidsp(IdxSp(0), nb_kinspecies + nb_fluidspecies);

    // Create a Field that contains masses of kinetic and fluid species (adiabatic species do not have a mass)
    host_t<DFieldMemSp> masses(dom_kinfluidsp);

    // fill the Field with masses of kinetic species
    for (IdxSp isp : dom_kinsp) {
        masses(isp) = kinetic_masses(isp);
    }

    // fill the Field with masses of fluid species
    for (IdxSp isp : dom_fluidsp) {
        masses(isp) = fluid_masses(isp);
    }
    ddc::init_discrete_space<Species>(std::move(charges), std::move(masses));

    /**
     * checks that the masses and charges of dom_allsp are well-ordered:
     * kinetic species first, then fluid species, then adiabatic species.
     */
    ddc::for_each(dom_allsp, [&](IdxSp const isp) {
        if (isp.uid() < nb_kinspecies) {
            EXPECT_EQ(ddc::discrete_space<Species>().charges()(isp), kinetic_charges(isp));
            EXPECT_EQ(ddc::discrete_space<Species>().masses()(isp), kinetic_masses(isp));

        } else if (nb_kinspecies <= isp.uid() && isp.uid() < nb_kinspecies + nb_fluidspecies) {
            EXPECT_EQ(ddc::discrete_space<Species>().charges()(isp), fluid_charges(isp));
            EXPECT_EQ(ddc::discrete_space<Species>().masses()(isp), fluid_masses(isp));

        } else {
            EXPECT_EQ(isp, dom_allsp.back());
            EXPECT_EQ(ddc::discrete_space<Species>().charges()(isp), charge_adiabspecies);
        }
    });
}
