// SPDX-License-Identifier: MIT

#include <cmath>
#include <iostream>

#include "ddc/discretization"

#include "species_info.hpp"

using std::sqrt, std::exp;

SpeciesInformation::SpeciesInformation(
        FieldSp<int> charge,
        FieldSp<double> mass,
        FieldSp<double> perturb_amplitude,
        FieldSp<int> perturb_mode)
    : m_charge(std::move(charge))
    , m_mass(std::move(mass))
    , m_perturb_amplitude(std::move(perturb_amplitude))
    , m_perturb_mode(std::move(perturb_mode))
{
}
