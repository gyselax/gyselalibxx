// SPDX-License-Identifier: MIT

#pragma once

#include <ddc/ddc.hpp>

class SpeciesInformation
{
    using mcoord_type = DiscreteCoordinate<SpeciesInformation>;

    using dvect_type = DiscreteVector<SpeciesInformation>;

    using idim_type = DiscreteDomain<SpeciesInformation>;

public:
    template <class MemorySpace>
    class Impl
    {
    private:
        // charge of the particles (kinetic + adiabatic)
        Chunk<int, idim_type> const m_charge;

        // mass of the particles of all kinetic species
        Chunk<double, idim_type> const m_mass;

        // Initial perturbation amplitude of all kinetic species
        Chunk<double, idim_type> const m_perturb_amplitude;

        // Initial perturbation mode of all kinetic species
        Chunk<int, idim_type> const m_perturb_mode;

    public:
        using ddim_type = SpeciesInformation;

        Impl(Chunk<int, idim_type> charge,
             Chunk<double, idim_type> mass,
             Chunk<double, idim_type> perturb_amplitude,
             Chunk<int, idim_type> perturb_mode)
            : m_charge(std::move(charge))
            , m_mass(std::move(mass))
            , m_perturb_amplitude(std::move(perturb_amplitude))
            , m_perturb_mode(std::move(perturb_mode))
        {
        }

        // Consider that the electron species is always at the 0 position
        mcoord_type ielec() const
        {
            return mcoord_type(0);
        }

        ChunkSpan<const int, idim_type> charges() const
        {
            return m_charge;
        }

        ChunkSpan<const double, idim_type> masses() const
        {
            return m_mass;
        }

        ChunkSpan<const double, idim_type> perturb_amplitudes() const
        {
            return m_perturb_amplitude;
        }

        ChunkSpan<const int, idim_type> perturb_modes() const
        {
            return m_perturb_mode;
        }
    };
};

inline DiscreteCoordinate<SpeciesInformation> ielec()
{
    return discretization<SpeciesInformation>().ielec();
}

inline int charge(DiscreteCoordinate<SpeciesInformation> isp)
{
    return discretization<SpeciesInformation>().charges()(isp);
}

inline double mass(DiscreteCoordinate<SpeciesInformation> isp)
{
    return discretization<SpeciesInformation>().masses()(isp);
}
