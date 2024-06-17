// SPDX-License-Identifier: MIT

#include <ddc/ddc.hpp>

#include <gtest/gtest.h>

#include <species_info.hpp>

TEST(SpeciesInfo, Ielec)
{
    IVectSp const nb_kinspecies(2);
    IDomainSp const dom_sp(IndexSp(0), nb_kinspecies);
    IndexSp my_iion = dom_sp.front();
    IndexSp my_ielec = dom_sp.back();

    FieldSp<int> charges(dom_sp);
    FieldSp<double> masses(dom_sp);
    charges(my_ielec) = -1;
    charges(my_iion) = 1;
    ddc::parallel_fill(masses, 1.);

    // Initialization of the distribution function
    ddc::init_discrete_space<IDimSp>(std::move(charges), std::move(masses));
    EXPECT_EQ(my_ielec, ielec());
}
