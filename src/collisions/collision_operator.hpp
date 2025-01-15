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

/**
 * @brief Initialise the collision operator defined in koliop.
 *
 * @param[in] mu_extent The number of points in the magnetic moment dimension.
 * @param[in] vpar_extent The number of points in the parallel velocity dimension.
 * @param[in] tor1_extent The number of points in the tor1 dimension.
 * @param[in] tor2_extent The number of points in the tor2 dimension.
 * @param[in] tor3_extent The number of points in the tor3 dimension.
 * @param[in] species_extent The number of species.
 * @param[in] collision_interspecies Boolean that is equal to true if inter-species collisions are taken into account.
 * @param[in] ir_sol_separatrix The index in the radial dimension where the SOL separatrix is found.
 * @param[in] hat_As The normalised masses for all species.
 * @param[in] hat_Zs The normalised charges for all species.
 * @param[in] mug The coordinates of the grid of magnetic moments.
 * @param[in] vparg The coordinates of the grid of parallel velocities.
 * @param[in] coeff_intdmu The quadrature coefficients to integrate on the grid of magnetic moments.
 * @param[in] coeff_intdvpar The quadrature coefficients to integrate on the grid of parallel velocities.
 * @param[in] coeff_AD_r
 * @param[in] mask_buffer_r The mask used to avoid applying collisions in certain radial region.
 * @param[in] mask_limiter The limiter mask defined over the polar slice.
 * @param[in] B_norm The norm of the magnetic field defined over the polar slice.
 * @param[in] Bstar_s Bstar defined on the coordinates (species,r,theta,vpar).
 *
 * @returns The handle which identifies the operator.
 */
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

/**
 * @brief Destructor for koliop.
 *
 * @param[in] operator_handle The handle which identifies the operator.
 */
void do_operator_deinitialisation(::koliop_Operator operator_handle);

}; // namespace detail

/**
 * @brief A class which computes the collision operator in (Sp,vpar,mu).
 */
template <class CollisionConfiguration>
class CollisionOperator /* : public IRightHandSide */
{
public:
    /**
     * @brief The type of collision configuration. Namely, what kind of geometry
     * the operator will ... operate on.
     */
    using CollisionConfigurationType = CollisionConfiguration;
    /**
     * @brief The distribution function that we should expect from the operator
     * user/caller.
     */
    using IdxRangeDistributionFunctionType =
            typename CollisionConfigurationType::IdxRangeDistributionFunctionType;

public:
    /**
     * @brief Construct a new Collision Operator object
     *
     * @param[in] collision_configuration This parametrise the operator. It contains the
     * data required for its initialization.
     */
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

    /**
     * @brief Destroy the Collision operator object.
     */
    ~CollisionOperator()
    {
        detail::do_operator_deinitialisation(m_operator_handle);
    }

    /**
     * @brief Apply the collision operator to the distribution functions of all
     * species on all species.
     *
     * @param[inout] all_f_distribution All the distribution function, depending
     * on the CollisionConfigurationType, the existence of the theta, r, phi
     * dimension may vary. At most, we have (species, phi, r, theta, vpar, mu)
     * in layout right.
     * @param[in] deltat_coll Collision time step.
     */
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
    /**
     * @brief Opaque C type representing the operator.
    */
    ::koliop_Operator m_operator_handle;
};
