// SPDX-License-Identifier: MIT

#pragma once

#include <paraconf.h>

#include "collision_common_configuration.hpp"
#include "paraconfpp.hpp"

/**
 * @brief Class to collect information to initialise the collision operator for a SpVparMu geometry.
 */
class CollisionConfiguration
{
public:
    /**
     * @brief Distribution function index range.
     */
    using IdxRangeDistributionFunctionType = IdxRange<Species, GridVpar, GridMu>;

    /**
     * @brief Container for the operator input data.
     */
    using CollisionConfigurationDataType = detail::CollisionConfigurationData<
            IdxRangeDistributionFunctionType,
            Species,
            detail::InternalSpoofGridPhi,
            detail::InternalSpoofGridR,
            detail::InternalSpoofGridTheta,
            GridVpar,
            GridMu>;

public:
    /**
     * @brief The constructor for the CollisionConfiguration class.
     *
     * @param yaml_input_file yaml namelist file object.
     * @param index_range_fdistribution index range of the whole distribution function.
     * @param coeff_intdmu quadrature coefficients.
     * @param coeff_intdvpar quadrature coefficients.
     * @param[in] B_norm The norm of the magnetic field defined over the polar slice.
     * @param[in] Bstar_s Bstar defined on the coordinates (species,r,theta,vpar).
     */
    CollisionConfiguration(
            PC_tree_t const& yaml_input_file,
            IdxRange<Species, GridVpar, GridMu> index_range_fdistribution,
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
        /* B_norm is a scalar when there is no Theta nor R dimensions.
         */
        ddc::parallel_fill(
                Kokkos::DefaultExecutionSpace(),
                get_field(m_operator_quantities.m_B_norm_allocation),
                B_norm);
        /* Koliop's Bstar_s is needed in a format layout different from the rest of the code. It is needed in layout
         * right [sp, theta, r, vpar] where as the input is in [sp, vpar]. We simply iterate over [sp, theta, r, vpar]
         * and copy the input [sp, vpar] for every theta, r.
         */
        CollisionConfigurationDataType::transpose_replicate_deep_copy(
                get_field(m_operator_quantities.m_Bstar_s_allocation),
                Bstar_s);
        /* When collisions only are tested, what matter is nu_ref * dt, we choose to only change dt.
         */
        if (m_operator_quantities.m_nu_ref != 1.0) {
            throw std::logic_error(
                    "nu_ref should be fixed at 1, use dt only to change collisionality.");
        }
        Kokkos::fence();
    }

    /**
     * @brief Can be used to obtain the quantities the operator needs for its initialisation, modification is not
     * authorised.
     *
     * @return const CollisionConfigurationDataType&
     */
    const CollisionConfigurationDataType& configuration() const
    {
        return m_operator_quantities;
    }

protected:
    /**
     * @brief Container for the arrays the quantities the operator needs.
     */
    CollisionConfigurationDataType m_operator_quantities;
};
