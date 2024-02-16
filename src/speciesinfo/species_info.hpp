// SPDX-License-Identifier: MIT

#pragma once

#include <ddc/ddc.hpp>

/// @brief Species discrete dimension to access constant attributes related to species.
class SpeciesInformation
{
public:
    /// alias of the discrete dimension
    using discrete_dimension_type = SpeciesInformation;

    /// alias of the discrete element of this discrete dimension
    using discrete_element_type = ddc::DiscreteElement<SpeciesInformation>;

    /// alias of the discrete domain of this discrete dimension
    using discrete_domain_type = ddc::DiscreteDomain<SpeciesInformation>;

    /// alias of the discrete vector of this discrete dimension
    using discrete_vector_type = ddc::DiscreteVector<SpeciesInformation>;

public:
    /// @brief Impl object storing attributes in `MemorySpace`.
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

        // workaround to access charges on the device
        ddc::ChunkView<int, discrete_domain_type, std::experimental::layout_right, MemorySpace>
                m_charge_view;

        // workaround to access masses on the device
        ddc::ChunkView<double, discrete_domain_type, std::experimental::layout_right, MemorySpace>
                m_mass_view;

        discrete_element_type m_ielec;

    public:
        /// alias of the discrete dimension
        using discrete_dimension_type = SpeciesInformation;

        /**
         * @brief Conversion constructor between different memory spaces.
         * @param[in] impl object from `OMemorySpace` that will be used to initialize this object on `MemorySpace`
         */
        template <class OMemorySpace>
        explicit Impl(Impl<OMemorySpace> const& impl)
            : m_charge(impl.m_charge.domain())
            , m_mass(impl.m_mass.domain())
            , m_perturb_amplitude(impl.m_perturb_amplitude.domain())
            , m_perturb_mode(impl.m_perturb_mode.domain())
            , m_ielec(impl.m_ielec)
        {
            m_charge_view = m_charge.span_cview();
            m_mass_view = m_mass.span_cview();
            ddc::deepcopy(m_charge, impl.m_charge);
            ddc::deepcopy(m_mass, impl.m_mass);
            ddc::deepcopy(m_perturb_amplitude, impl.m_perturb_amplitude);
            ddc::deepcopy(m_perturb_mode, impl.m_perturb_mode);
        }

        /**
         * @brief Main constructor taking all attributes
         * @param[in] charge array storing both kinetic and adiabatic charges
         * @param[in] mass array storing both kinetic and adiabatic masses
         * @param[in] perturb_amplitude array storing kinetic pertubation amplitudes
         * @param[in] perturb_mode array storing kinetic pertubation modes
         */
        Impl(ddc::Chunk<int, discrete_domain_type, ddc::KokkosAllocator<int, MemorySpace>> charge,
             ddc::Chunk<double, discrete_domain_type, ddc::KokkosAllocator<double, MemorySpace>>
                     mass,
             ddc::Chunk<double, discrete_domain_type, ddc::KokkosAllocator<double, MemorySpace>>
                     perturb_amplitude,
             ddc::Chunk<int, discrete_domain_type, ddc::KokkosAllocator<int, MemorySpace>>
                     perturb_mode)
            : m_charge(std::move(charge))
            , m_mass(std::move(mass))
            , m_perturb_amplitude(std::move(perturb_amplitude))
            , m_perturb_mode(std::move(perturb_mode))
        {
            m_charge_view = m_charge.span_cview();
            m_mass_view = m_mass.span_cview();
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

        /// @return the discrete element representing the electron species
        KOKKOS_FUNCTION discrete_element_type ielec() const
        {
            return m_ielec;
        }

        /**
         * @param[in] isp a discrete element of either a kinetic or adiabatic species
         * @return the charge associated to the discrete element
         */
        KOKKOS_FUNCTION int charge(discrete_element_type const isp) const
        {
            return m_charge_view(isp);
        }

        /**
         * @param[in] isp a discrete element of either a kinetic or adiabatic species
         * @return the mass associated to the discrete element
         */
        KOKKOS_FUNCTION double mass(discrete_element_type const isp) const
        {
            return m_mass_view(isp);
        }

        /// @return kinetic and adiabatic charges array
        auto charges() const
        {
            return m_charge.span_view();
        }

        /// @return kinetic and adiabatic masses array
        auto masses() const
        {
            return m_mass.span_view();
        }

        /// @return kinetic perturbation amplitudes array
        auto perturb_amplitudes() const
        {
            return m_perturb_amplitude.span_view();
        }

        /// @return kinetic perturbation modes array
        auto perturb_modes() const
        {
            return m_perturb_mode.span_view();
        }
    };
};

/// @return the discrete element representing the electron species
KOKKOS_INLINE_FUNCTION ddc::DiscreteElement<SpeciesInformation> ielec()
{
    return ddc::discrete_space<SpeciesInformation>().ielec();
}

/**
 * @param[in] isp a discrete element of either a kinetic or adiabatic species
 * @return the charge associated to the discrete element
 */
KOKKOS_INLINE_FUNCTION int charge(ddc::DiscreteElement<SpeciesInformation> const isp)
{
    return ddc::discrete_space<SpeciesInformation>().charge(isp);
}

/**
 * @param[in] isp a discrete element of either a kinetic or adiabatic species
 * @return the mass associated to the discrete element
 */
KOKKOS_INLINE_FUNCTION double mass(ddc::DiscreteElement<SpeciesInformation> const isp)
{
    return ddc::discrete_space<SpeciesInformation>().mass(isp);
}
