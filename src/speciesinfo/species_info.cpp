// SPDX-License-Identifier: MIT

#include <ddc/ddc.hpp>

#include "species_info.hpp"

SpeciesInformation::SpeciesInformation(
        Chunk<int, idim_type> charge,
        Chunk<double, idim_type> mass,
        Chunk<double, idim_type> perturb_amplitude,
        Chunk<int, idim_type> perturb_mode)
    : m_charge(std::move(charge))
    , m_mass(std::move(mass))
    , m_perturb_amplitude(std::move(perturb_amplitude))
    , m_perturb_mode(std::move(perturb_mode))
{
}
