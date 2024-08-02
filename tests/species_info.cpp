// SPDX-License-Identifier: MIT

#include <ddc/ddc.hpp>

#include <gtest/gtest.h>

#include <species_info.hpp>

TEST(SpeciesInfo, Ielec)
{
    IdxStepSp const nb_kinspecies(2);
    IdxRangeSp const dom_sp(IdxSp(0), nb_kinspecies);
    IdxSp my_iion = dom_sp.front();
    IdxSp my_ielec = dom_sp.back();

    host_t<DFieldMemSp> charges(dom_sp);
    host_t<DFieldMemSp> masses(dom_sp);
    charges(my_ielec) = -1.;
    charges(my_iion) = 1.;
    ddc::parallel_fill(masses, 1.);

    // Initialization of the distribution function
    ddc::init_discrete_space<Species>(std::move(charges), std::move(masses));
    EXPECT_EQ(my_ielec, ielec());
}
