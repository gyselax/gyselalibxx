// SPDX-License-Identifier: MIT
#pragma once

#include <ddc/ddc.hpp>

#include <KOLIOP/collisions.h>

#include "assert.hpp"
#include "ddc_alias_inline_functions.hpp"
#include "ddc_aliases.hpp"
#include "ddc_helper.hpp"
#include "species_info.hpp"

/**
 * @brief A class which computes the collision operator in (Sp,vpar,mu).
 */
template <class CollisionConfiguration>
class CollisionOperator /* : public IRightHandSide */
{
public:
    using IdxRangeDistributionFunctionType =
            typename CollisionConfiguration::IdxRangeDistributionFunctionType;

    using CollisionConfigurationDataType =
            typename CollisionConfiguration::CollisionConfigurationDataType;

public:
    /**
     * @brief Construct a new Collision Operator object
     *
     * @param[in] collision_configuration This parametrise the operator. It contains the data required for its
     * initialisation.
     */
    explicit CollisionOperator(CollisionConfiguration const& collision_configuration)
        : m_collision_configuration {collision_configuration}
        , m_operator {
                  m_collision_configuration.configuration().m_theta_extent,
                  m_collision_configuration.configuration().m_r_extent,
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

    /**
     * @brief Apply the collision operator to the distribution functions of all species on all species.
     *
     * @param[inout] all_f_distribution All the distribution function, depending on the CollisionConfiguration, the
     * existence of the theta, r, phi dimension may vary. At most, we have (species, phi, r, theta, vpar, mu) in layout
     * right.
     * @param[in] deltat_coll Collision time step.
     */
    void operator()(DField<IdxRangeDistributionFunctionType> all_f_distribution, double deltat_coll)
    {
        m_operator.Collisions(
                deltat_coll,
                m_collision_configuration.configuration().cast_fdistrib(
                        all_f_distribution.allocation_kokkos_view()));
    }

protected:
    /**
     * @brief The collision configuration.
    */
    CollisionConfiguration const& m_collision_configuration;

    /**
     * @brief Opaque C type representing the operator.
    */
    koliop::collisions::Operator m_operator;
};
