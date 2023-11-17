// SPDX-License-Identifier: MIT

#pragma once

#include <ddc/ddc.hpp>

class SpeciesInformation
{
public:
    using discrete_dimension_type = SpeciesInformation;

    using discrete_element_type = ddc::DiscreteElement<SpeciesInformation>;

    using discrete_domain_type = ddc::DiscreteDomain<SpeciesInformation>;

    using discrete_vector_type = ddc::DiscreteVector<SpeciesInformation>;

public:
    static KOKKOS_FUNCTION constexpr std::size_t rank()
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
        ddc::Chunk<int, discrete_domain_type, ddc::KokkosAllocator<int, MemorySpace>> m_charge;

        // mass of the particles of all kinetic species
        ddc::Chunk<double, discrete_domain_type, ddc::KokkosAllocator<double, MemorySpace>> m_mass;

        // Initial perturbation amplitude of all kinetic species
        ddc::Chunk<double, discrete_domain_type, ddc::KokkosAllocator<double, MemorySpace>>
                m_perturb_amplitude;

        // Initial perturbation mode of all kinetic species
        ddc::Chunk<int, discrete_domain_type, ddc::KokkosAllocator<int, MemorySpace>>
                m_perturb_mode;

        discrete_element_type m_ielec;

    public:
        using discrete_dimension_type = SpeciesInformation;

        template <class OMemorySpace>
        explicit Impl(Impl<OMemorySpace> const& impl)
            : m_charge(impl.m_charge.domain())
            , m_mass(impl.m_mass.domain())
            , m_perturb_amplitude(impl.m_perturb_amplitude.domain())
            , m_perturb_mode(impl.m_perturb_mode.domain())
            , m_ielec(impl.m_ielec)
        {
            ddc::deepcopy(m_charge, impl.m_charge);
            ddc::deepcopy(m_mass, impl.m_mass);
            ddc::deepcopy(m_perturb_amplitude, impl.m_perturb_amplitude);
            ddc::deepcopy(m_perturb_mode, impl.m_perturb_mode);
        }

        Impl(ddc::Chunk<int, discrete_domain_type> charge,
             ddc::Chunk<double, discrete_domain_type> mass,
             ddc::Chunk<double, discrete_domain_type> perturb_amplitude,
             ddc::Chunk<int, discrete_domain_type> perturb_mode)
            : m_charge(std::move(charge))
            , m_mass(std::move(mass))
            , m_perturb_amplitude(std::move(perturb_amplitude))
            , m_perturb_mode(std::move(perturb_mode))
        {
            assert(charge.size() >= 2);
            bool electron_found = false;
            for (discrete_element_type const isp : m_charge.domain()) {
                if (m_charge(isp) == -1) {
                    electron_found = true;
                    m_ielec = isp;
                }
            }
            if (!electron_found) {
                throw std::runtime_error("electron not found");
            }
        }

        KOKKOS_FUNCTION discrete_element_type ielec() const
        {
            return m_ielec;
        }

        ddc::ChunkSpan<const int, discrete_domain_type> charges() const
        {
            return m_charge;
        }

        ddc::ChunkSpan<const double, discrete_domain_type> masses() const
        {
            return m_mass;
        }

        ddc::ChunkSpan<const double, discrete_domain_type> perturb_amplitudes() const
        {
            return m_perturb_amplitude;
        }

        ddc::ChunkSpan<const int, discrete_domain_type> perturb_modes() const
        {
            return m_perturb_mode;
        }
    };
};

KOKKOS_INLINE_FUNCTION ddc::DiscreteElement<SpeciesInformation> ielec()
{
    return ddc::discrete_space<SpeciesInformation>().ielec();
}

inline int charge(ddc::DiscreteElement<SpeciesInformation> const isp)
{
    return ddc::discrete_space<SpeciesInformation>().charges()(isp);
}

inline double mass(ddc::DiscreteElement<SpeciesInformation> const isp)
{
    return ddc::discrete_space<SpeciesInformation>().masses()(isp);
}
