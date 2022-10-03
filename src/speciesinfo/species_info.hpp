// SPDX-License-Identifier: MIT

#pragma once

#include <ddc/ddc.hpp>

class SpeciesInformation
{
public:
    using discrete_dimension_type = SpeciesInformation;

    using discrete_element_type = DiscreteElement<SpeciesInformation>;

    using discrete_domain_type = DiscreteDomain<SpeciesInformation>;

    using discrete_vector_type = DiscreteVector<SpeciesInformation>;

public:
    static constexpr std::size_t rank()
    {
        return 1;
    }

    template <class MemorySpace>
    class Impl
    {
        template <class OMemorySpace>
        friend class Impl;

    private:
        // charge of the particles (kinetic + adiabatic)
        Chunk<int, discrete_domain_type, KokkosAllocator<int, MemorySpace>> m_charge;

        // mass of the particles of all kinetic species
        Chunk<double, discrete_domain_type, KokkosAllocator<double, MemorySpace>> m_mass;

        // Initial perturbation amplitude of all kinetic species
        Chunk<double, discrete_domain_type, KokkosAllocator<double, MemorySpace>>
                m_perturb_amplitude;

        // Initial perturbation mode of all kinetic species
        Chunk<int, discrete_domain_type, KokkosAllocator<int, MemorySpace>> m_perturb_mode;

    public:
        using discrete_dimension_type = SpeciesInformation;

        template <class OMemorySpace>
        explicit Impl(Impl<OMemorySpace> const& impl)
            : m_charge(impl.m_charge.domain())
            , m_mass(impl.m_mass.domain())
            , m_perturb_amplitude(impl.m_perturb_amplitude.domain())
            , m_perturb_mode(impl.m_perturb_mode.domain())
        {
            deepcopy(m_charge, impl.m_charge);
            deepcopy(m_mass, impl.m_mass);
            deepcopy(m_perturb_amplitude, impl.m_perturb_amplitude);
            deepcopy(m_perturb_mode, impl.m_perturb_mode);
        }

        Impl(Chunk<int, discrete_domain_type> charge,
             Chunk<double, discrete_domain_type> mass,
             Chunk<double, discrete_domain_type> perturb_amplitude,
             Chunk<int, discrete_domain_type> perturb_mode)
            : m_charge(std::move(charge))
            , m_mass(std::move(mass))
            , m_perturb_amplitude(std::move(perturb_amplitude))
            , m_perturb_mode(std::move(perturb_mode))
        {
        }

        // Consider that the electron species is always at the 0 position
        discrete_element_type ielec() const
        {
            return discrete_element_type(0);
        }

        ChunkSpan<const int, discrete_domain_type> charges() const
        {
            return m_charge;
        }

        ChunkSpan<const double, discrete_domain_type> masses() const
        {
            return m_mass;
        }

        ChunkSpan<const double, discrete_domain_type> perturb_amplitudes() const
        {
            return m_perturb_amplitude;
        }

        ChunkSpan<const int, discrete_domain_type> perturb_modes() const
        {
            return m_perturb_mode;
        }
    };
};

inline DiscreteElement<SpeciesInformation> ielec()
{
    return discrete_space<SpeciesInformation>().ielec();
}

inline int charge(DiscreteElement<SpeciesInformation> const isp)
{
    return discrete_space<SpeciesInformation>().charges()(isp);
}

inline double mass(DiscreteElement<SpeciesInformation> const isp)
{
    return discrete_space<SpeciesInformation>().masses()(isp);
}
