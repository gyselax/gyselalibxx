

# File collision\_operator.hpp

[**File List**](files.md) **>** [**collisions**](dir_64163437c27c8707f17f92558da22106.md) **>** [**collision\_operator.hpp**](collision__operator_8hpp.md)

[Go to the documentation of this file](collision__operator_8hpp.md)


```C++
// SPDX-License-Identifier: MIT
#pragma once

#include <ddc/ddc.hpp>

#include <KOLIOP/collisions.h>

#include "assert.hpp"
#include "ddc_alias_inline_functions.hpp"
#include "ddc_aliases.hpp"
#include "ddc_helper.hpp"
#include "species_info.hpp"

template <class CollisionConfiguration>
class CollisionOperator /* : public IRightHandSide */
{
public:
    using IdxRangeDistributionFunctionType =
            typename CollisionConfiguration::IdxRangeDistributionFunctionType;
    using CollisionConfigurationDataType =
            typename CollisionConfiguration::CollisionConfigurationDataType;

public:
    explicit CollisionOperator(CollisionConfiguration const& collision_configuration)
        : m_collision_configuration {collision_configuration}
        , m_operator {
                  m_collision_configuration.configuration().m_r_extent,
                  m_collision_configuration.configuration().m_theta_extent,
                  m_collision_configuration.configuration().m_phi_extent,
                  m_collision_configuration.configuration().m_collisions_interspecies,
                  CollisionConfigurationDataType::as_layout_left(
                          m_collision_configuration.configuration()
                                  .m_hat_As_allocation.allocation_kokkos_view()),
                  CollisionConfigurationDataType::as_layout_left(
                          m_collision_configuration.configuration()
                                  .m_hat_Zs_allocation.allocation_kokkos_view()),
                  CollisionConfigurationDataType::as_layout_left(
                          m_collision_configuration.configuration()
                                  .m_mug_allocation.allocation_kokkos_view()),
                  CollisionConfigurationDataType::as_layout_left(
                          m_collision_configuration.configuration()
                                  .m_vparg_allocation.allocation_kokkos_view()),
                  CollisionConfigurationDataType::as_layout_left(
                          m_collision_configuration.configuration()
                                  .m_coeff_intdmu.allocation_kokkos_view()),
                  CollisionConfigurationDataType::as_layout_left(
                          m_collision_configuration.configuration()
                                  .m_coeff_intdvpar.allocation_kokkos_view()),
                  CollisionConfigurationDataType::as_layout_left(
                          m_collision_configuration.configuration()
                                  .m_coeff_AD_allocation.allocation_kokkos_view()),
                  CollisionConfigurationDataType::as_layout_left(
                          m_collision_configuration.configuration()
                                  .m_B_norm_allocation.allocation_kokkos_view()),
                  CollisionConfigurationDataType::as_layout_left(
                          m_collision_configuration.configuration()
                                  .m_Bstar_s_allocation.allocation_kokkos_view())}
    {
    }

    void operator()(DField<IdxRangeDistributionFunctionType> all_f_distribution, double deltat_coll)
    {
        m_operator.Collisions(
                deltat_coll,
                m_collision_configuration.configuration().cast_fdistrib(
                        all_f_distribution.allocation_kokkos_view()));
    }

protected:
    CollisionConfiguration const& m_collision_configuration;

    koliop::collisions::Operator m_operator;
};
```


