

# File collision\_configuration.hpp

[**File List**](files.md) **>** [**collisions**](dir_d2f3522bfaadb2833cfc3c18da51fc33.md) **>** [**collision\_configuration.hpp**](collision__configuration_8hpp.md)

[Go to the documentation of this file](collision__configuration_8hpp.md)


```C++
// SPDX-License-Identifier: MIT

#pragma once

#include <paraconf.h>

#include "collision_common_configuration.hpp"
#include "geometry_vpar_mu.hpp"
#include "paraconfpp.hpp"

class CollisionConfigurationSpV2D
{
public:
    using IdxRangeDistributionFunctionType = IdxRange<Species, GridVpar, GridMu>;

    using CollisionConfigurationDataType = detail::CollisionConfigurationData<
            IdxRangeDistributionFunctionType,
            detail::InternalSpoofGridTor3,
            detail::InternalSpoofGridTor2,
            detail::InternalSpoofGridTor1,
            GridVpar,
            GridMu>;

public:
    CollisionConfigurationSpV2D(
            PC_tree_t const& yaml_input_file,
            IdxRangeDistributionFunctionType index_range_fdistribution,
            DConstField<IdxRange<GridMu>> coeff_intdmu,
            DConstField<IdxRange<GridVpar>> coeff_intdvpar,
            double B_norm,
            DConstField<IdxRange<Species, GridVpar>> Bstar_s)
        : m_operator_quantities {
                yaml_input_file,
                index_range_fdistribution,
                coeff_intdmu,
                coeff_intdvpar}
    {
        ddc::parallel_fill(
                Kokkos::DefaultExecutionSpace(),
                get_field(m_operator_quantities.m_B_norm_allocation),
                B_norm);
        CollisionConfigurationDataType::transpose_replicate_deep_copy(
                get_field(m_operator_quantities.m_Bstar_s_allocation),
                Bstar_s);
        /* When collisions only are tested, what matter is nu_ref * dt, we choose to only change dt.
         */
        if (m_operator_quantities.m_nu_ref != 1.0) {
            throw std::logic_error(
                    "nu_ref should be fixed at 1, use dt only to change collisionality.");
        }
        // nu_ref corresponds to the collision frequency normalised to the electron cyclotron frequency (as always in
        // Gysela), for a normalised density/temperature of 1. We use nu_ref to construct AD, nu_ref is physically more
        // intuitive than AD, in addition, nu_ref is a scalar where AD can be dependent on the geometry like the radial
        // profile we build below. Having AD on a radial profile can allow masking which has its uses for preventing
        // badly behaving collisions on the colder plasma edges. sqrt(2) is a subtlety between the intra-species and
        // inter-species collision frequencies.
        ddc::parallel_fill(
                Kokkos::DefaultExecutionSpace(),
                get_field(m_operator_quantities.m_coeff_AD_allocation),
                m_operator_quantities.m_nu_ref * Kokkos::sqrt(2.0));
        Kokkos::fence();
    }

    const CollisionConfigurationDataType& configuration() const
    {
        return m_operator_quantities;
    }

protected:
    CollisionConfigurationDataType m_operator_quantities;
};
```


