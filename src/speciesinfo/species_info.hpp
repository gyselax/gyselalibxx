// SPDX-License-Identifier: MIT

#pragma once

#include "geometry.hpp"

class SpeciesInformation
{
private:
    // charge of the particles (kinetic + adiabatic)
    FieldSp<int> const m_charge;

    // mass of the particles of all kinetic species
    FieldSp<double> const m_mass;

    // Initial perturbation amplitude of all kinetic species
    FieldSp<double> const m_perturb_amplitude;

    // Initial perturbation mode of all kinetic species
    FieldSp<int> const m_perturb_mode;

public:
    SpeciesInformation(
            FieldSp<int> charge,
            FieldSp<double> mass,
            FieldSp<double> perturb_amplitude,
            FieldSp<int> perturb_mode);

    // Consider that the electron species is always at the 0 position
    IndexSp ielec() const
    {
        return IndexSp(0);
    }

    ViewSp<int> charge() const
    {
        return m_charge;
    }

    ViewSp<double> mass() const
    {
        return m_mass;
    }

    ViewSp<double> perturb_amplitude() const
    {
        return m_perturb_amplitude;
    }

    ViewSp<int> perturb_mode() const
    {
        return m_perturb_mode;
    }
};
