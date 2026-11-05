// SPDX-License-Identifier: MIT

#pragma once

#include <paraconf.h>

#include "collision_common_configuration.hpp"
#include "geometry_vpar_mu.hpp"
#include "paraconfpp.hpp"

/**
 * @brief Class to collect information to initialise the collision operator for a SpVparMu geometry.
 */
class CollisionConfigurationSpV2D
{
public:
    /**
     * @brief Distribution function index range, define the distribution function layout (layout right).
     */
    using IdxRangeDistributionFunctionType = IdxRange<Species, GridVpar, GridMu>;

    /**
     * @brief Container for the operator input data.
     */
    using CollisionConfigurationDataType = detail::CollisionConfigurationData<
            IdxRangeDistributionFunctionType,
            detail::InternalSpoofGridTor3,
            detail::InternalSpoofGridTor2,
            detail::InternalSpoofGridTor1,
            GridVpar,
            GridMu>;

public:
    /**
     * @brief The constructor for the CollisionConfigurationSpV2D class.
     *
     * @param yaml_input_file yaml namelist file object.
     * @param index_range_fdistribution index range of the whole distribution function.
     * @param coeff_intdmu quadrature coefficients.
     * @param coeff_intdvpar quadrature coefficients.
     * @param[in] B_norm The norm of the magnetic field.
     * @param[in] Bstar_s Bstar defined on the coordinates.
     */
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
