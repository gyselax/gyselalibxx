

# File collision\_operator.hpp

[**File List**](files.md) **>** [**collisions**](dir_64163437c27c8707f17f92558da22106.md) **>** [**collision\_operator.hpp**](collision__operator_8hpp.md)

[Go to the documentation of this file](collision__operator_8hpp.md)


```C++
// SPDX-License-Identifier: MIT
#pragma once

#include <ddc/ddc.hpp>

#include "assert.hpp"
#include "ddc_alias_inline_functions.hpp"
#include "ddc_aliases.hpp"
#include "ddc_helper.hpp"
#include "species_info.hpp"

// SPDX-License-Identifier: MIT
#pragma once

#include <ddc/ddc.hpp>

#include <KOLIOP/koliop.h>

namespace detail {

::koliop_Operator do_operator_initialisation(
        std::size_t mu_extent,
        std::size_t vpar_extent,
        std::size_t tor1_extent,
        std::size_t tor2_extent,
        std::size_t tor3_extent,
        std::size_t species_extent,
        std::int8_t collision_interspecies,
        std::int64_t ir_sol_separatrix,
        double const* hat_As,
        double const* hat_Zs,
        double const* mug,
        double const* vparg,
        double const* coeff_intdmu,
        double const* coeff_intdvpar,
        double const* coeff_AD_r,
        double const* mask_buffer_r,
        double const* mask_limiter,
        double const* B_norm,
        double const* Bstar_s);

void do_operator_deinitialisation(::koliop_Operator operator_handle);

}; // namespace detail

template <class CollisionConfiguration>
class CollisionOperator /* : public IRightHandSide */
{
public:
    using CollisionConfigurationType = CollisionConfiguration;
    using IdxRangeDistributionFunctionType =
            typename CollisionConfigurationType::IdxRangeDistributionFunctionType;

public:
    explicit CollisionOperator(CollisionConfigurationType const& collision_configuration)
        : m_operator_handle {detail::do_operator_initialisation(
                collision_configuration.configuration().m_mu_extent,
                collision_configuration.configuration().m_vpar_extent,
                collision_configuration.configuration().m_r_extent,
                collision_configuration.configuration().m_theta_extent,
                collision_configuration.configuration().m_phi_extent,
                collision_configuration.configuration().m_sp_extent,
                collision_configuration.configuration().m_collisions_interspecies,
                /* One less than the number of point in the r direction. */
                0 /* Null offset */ + collision_configuration.configuration().m_r_extent - 1,
                collision_configuration.configuration().m_hat_As_allocation.data_handle(),
                collision_configuration.configuration().m_hat_Zs_allocation.data_handle(),
                collision_configuration.configuration().m_mug_allocation.data_handle(),
                collision_configuration.configuration().m_vparg_allocation.data_handle(),
                collision_configuration.configuration().m_coeff_intdmu.data_handle(),
                collision_configuration.configuration().m_coeff_intdvpar.data_handle(),
                collision_configuration.configuration().m_coeff_AD_allocation.data_handle(),
                collision_configuration.configuration().m_mask_buffer_r_allocation.data_handle(),
                collision_configuration.configuration().m_mask_LIM_allocation.data_handle(),
                collision_configuration.configuration().m_B_norm_allocation.data_handle(),
                collision_configuration.configuration().m_Bstar_s_allocation.data_handle())}
    {
        // EMPTY
    }

    ~CollisionOperator()
    {
        detail::do_operator_deinitialisation(m_operator_handle);
    }

    void operator()(DField<IdxRangeDistributionFunctionType> all_f_distribution, double deltat_coll)
            const
    {
        if (::koliop_Collision(
                    static_cast<::koliop_Operator>(m_operator_handle),
                    deltat_coll,
                    all_f_distribution.data_handle())
            != KOLIOP_STATUS_SUCCESS) {
            GSLX_ASSERT(false);
        }

        // NOTE: Koliop is asynchronous, fence to ensure the operator ended.
        if (::koliop_Fence(static_cast<::koliop_Operator>(m_operator_handle))
            != KOLIOP_STATUS_SUCCESS) {
            GSLX_ASSERT(false);
        }
    }

protected:
    ::koliop_Operator m_operator_handle;
};
```


