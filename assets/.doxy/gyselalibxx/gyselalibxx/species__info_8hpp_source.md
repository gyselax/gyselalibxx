

# File species\_info.hpp

[**File List**](files.md) **>** [**speciesinfo**](dir_661be8452a62f1b4720eb6eb57123ae7.md) **>** [**species\_info.hpp**](species__info_8hpp.md)

[Go to the documentation of this file](species__info_8hpp.md)


```C++
// SPDX-License-Identifier: MIT

#pragma once

#include <ddc/ddc.hpp>

#include "ddc_alias_inline_functions.hpp"
#include "ddc_aliases.hpp"
#include "ddc_helper.hpp"

class SpeciesInformation
{
public:
    using discrete_dimension_type = SpeciesInformation;

public:
    template <class Grid1D, class MemorySpace>
    class Impl
    {
        template <class ODDim, class OMemorySpace>
        friend class Impl;

        using discrete_element_type = Idx<Grid1D>;

        using discrete_domain_type = IdxRange<Grid1D>;
        using index_range_type = discrete_domain_type;

        using discrete_vector_type = IdxStep<Grid1D>;

    private:
        // charge of the particles (kinetic + adiabatic)
        DFieldMem<index_range_type, MemorySpace> m_charge;

        // mass of the particles of all kinetic species
        DFieldMem<index_range_type, MemorySpace> m_mass;

        // workaround to access charges on the device
        DConstField<index_range_type, MemorySpace> m_charge_view;

        // workaround to access masses on the device
        DConstField<index_range_type, MemorySpace> m_mass_view;

        discrete_element_type m_ielec;

    public:
        using discrete_dimension_type = SpeciesInformation;

        template <class OMemorySpace>
        explicit Impl(Impl<Grid1D, OMemorySpace> const& impl)
            : m_charge(get_idx_range(impl.m_charge))
            , m_mass(get_idx_range(impl.m_mass))
            , m_ielec(impl.m_ielec)
        {
            m_charge_view = get_const_field(m_charge);
            m_mass_view = get_const_field(m_mass);
            ddc::parallel_deepcopy(m_charge, impl.m_charge);
            ddc::parallel_deepcopy(m_mass, impl.m_mass);
        }

        Impl(DFieldMem<index_range_type, MemorySpace> charge,
             DFieldMem<index_range_type, MemorySpace> mass)
            : m_charge(std::move(charge))
            , m_mass(std::move(mass))
        {
            m_charge_view = get_const_field(m_charge);
            m_mass_view = get_const_field(m_mass);
            assert(charge.size() >= 2);
            bool electron_found = false;
            for (discrete_element_type const isp : get_idx_range(m_charge)) {
                if (m_charge(isp) == -1.) {
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

        KOKKOS_FUNCTION double charge(discrete_element_type const isp) const
        {
            return m_charge_view(isp);
        }

        KOKKOS_FUNCTION double mass(discrete_element_type const isp) const
        {
            return m_mass_view(isp);
        }

        auto charges() const
        {
            return get_const_field(m_charge);
        }

        auto masses() const
        {
            return get_const_field(m_mass);
        }
    };
};

template <class Grid1D>
inline constexpr bool is_species_information_v
        = std::is_same_v<typename Grid1D::discrete_dimension_type, SpeciesInformation>;

// Species dimension
struct Species : SpeciesInformation
{
};

KOKKOS_INLINE_FUNCTION Idx<Species> ielec()
{
    return ddc::discrete_space<Species>().ielec();
}

KOKKOS_INLINE_FUNCTION double charge(Idx<Species> const isp)
{
    return ddc::discrete_space<Species>().charge(isp);
}

KOKKOS_INLINE_FUNCTION double mass(Idx<Species> const isp)
{
    return ddc::discrete_space<Species>().mass(isp);
}

using IdxSp = Idx<Species>;
using IdxRangeSp = IdxRange<Species>;
using IdxStepSp = IdxStep<Species>;

template <class ElementType>
using FieldMemSp = FieldMem<ElementType, IdxRangeSp>;
using DFieldMemSp = FieldMemSp<double>;
using IFieldMemSp = FieldMemSp<int>;

template <class ElementType>
using ConstFieldSp = ConstField<ElementType, IdxRangeSp>;
using DConstFieldSp = ConstFieldSp<double>;

template <class ElementType>
using FieldSp = Field<ElementType, IdxRangeSp>;
using DFieldSp = FieldSp<double>;
```


