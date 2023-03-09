// SPDX-License-Identifier: MIT

#include <ddc/ddc.hpp>

#include <gtest/gtest.h>

#include <geometry.hpp>
#include <species_info.hpp>

TEST(SpeciesInfo, Ielec)
{
    IVectSp const nb_kinspecies(2);
    IDomainSp const dom_sp(IndexSp(0), nb_kinspecies);
    IndexSp my_iion = dom_sp.front();
    IndexSp my_ielec = dom_sp.back();

    FieldSp<int> charges(dom_sp);
    FieldSp<double> masses(dom_sp);
    FieldSp<int> init_perturb_mode(dom_sp);
    FieldSp<double> init_perturb_amplitude(dom_sp);
    charges(my_ielec) = -1;
    charges(my_iion) = 1;
    ddc::fill(masses, 1.);
    ddc::fill(init_perturb_amplitude, 0);
    ddc::fill(init_perturb_mode, 0);

    // Initialization of the distribution function
    ddc::init_discrete_space<IDimSp>(
            std::move(charges),
            std::move(masses),
            std::move(init_perturb_amplitude),
            std::move(init_perturb_mode));
    EXPECT_EQ(my_ielec, ielec());
}
