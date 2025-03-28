

# File collision\_configuration.hpp

[**File List**](files.md) **>** [**collisions**](dir_d2f3522bfaadb2833cfc3c18da51fc33.md) **>** [**collision\_configuration.hpp**](collision__configuration_8hpp.md)

[Go to the documentation of this file](collision__configuration_8hpp.md)


```C++
// SPDX-License-Identifier: MIT

#pragma once

#include <paraconf.h>

#include "collision_common_configuration.hpp"
#include "paraconfpp.hpp"

template <class GridSp, class GridVpar, class GridMu>
class CollisionConfiguration
{
public:
    using GridMuType = GridMu;
    using GridVparType = GridVpar;
    using GridThetaType = detail::InternalSpoofGridTheta;
    using GridRType = detail::InternalSpoofGridR;
    using GridPhiType = detail::InternalSpoofGridPhi;
    using GridSpType = GridSp;

    using IdxRangeMuType = IdxRange<GridMuType>;
    using IdxRangeVparType = IdxRange<GridVparType>;
    using IdxRangeSpVparType = IdxRange<GridSpType, GridVparType>;
    using IdxRangeDistributionFunctionType = IdxRange<GridSpType, GridVparType, GridMuType>;

    using CollisionConfigurationDataType = detail::CollisionConfigurationData<
            IdxRangeDistributionFunctionType,
            GridSpType,
            GridPhiType,
            GridRType,
            GridThetaType,
            GridVparType,
            GridMuType>;

public:
    CollisionConfiguration(
            PC_tree_t const& yaml_input_file,
            IdxRangeDistributionFunctionType index_range_fdistribution,
            DConstField<IdxRangeMuType> coeff_intdmu,
            DConstField<IdxRangeVparType> coeff_intdvpar,
            double B_norm,
            DConstField<IdxRangeSpVparType> Bstar_s)
        : m_operator_quantities {yaml_input_file, index_range_fdistribution, coeff_intdmu, coeff_intdvpar}
        , m_nustar0 {PCpp_double(yaml_input_file, ".Collisions.nustar0_rpeak")}
    {
        /* nustar fixed to 1. => normalisation of time is equal to 1./nustar
        */
        if (m_nustar0 != 1.0) {
            throw std::invalid_argument("nustar0 must be equal to 1");
        }

        do_configuration_data_initialisation(B_norm, Bstar_s);
    }

    const CollisionConfigurationDataType& configuration() const
    {
        return m_operator_quantities;
    }

protected:
    void do_configuration_data_initialisation(
            double input_B_norm,
            DConstField<IdxRangeSpVparType> input_Bstar_s)
    {
        {
            /* input_B_norm is a scalar when there is no Theta nor R dimensions.
             */
            ddc::parallel_fill(
                    Kokkos::DefaultExecutionSpace(),
                    get_field(m_operator_quantities.m_B_norm_allocation),
                    input_B_norm);
        }

        {
            /* Koliop's Bstar_s is needed in a format layout different from the
             * rest of the code. It is needed in layout right
             * [sp, theta, r, vpar] where as the input is in [sp, vpar].
             * We simply iterate over [sp, theta, r, vpar] and copy the input
             * [sp, vpar] for every theta, r.
             */
            CollisionConfigurationDataType::transpose_replicate_deep_copy(
                    get_field(m_operator_quantities.m_Bstar_s_allocation),
                    input_Bstar_s);
        }

        {
            /* Computation of the radial profile of ADi that is required for
             * koliop interface
             * @f[ AD(rpeak) = \sqrt(2)*eps_rpeak^(3/2)/(q_rpeak R0)*nustar0_rpeak @f]
             * eps_rpeak = eps(rpeak) with eps(r) = r/R0 is the aspect ratio.
             *
             * In this specific case, where AD is a scalar all values are forced
             * to 1, so coeff_AD_atrpeak = sqrt(2)
             */
            double const coeff_AD = Kokkos::sqrt(2.0);

            /* [TODO] According to Peter D., the places where coeff_AD is needed
             * can resynthesise the quantities from a scalar and some other
             * fields that koliop is already aware of. This change has not
             * landed in Gysela Fortran yet. For koliop to make do with a scalar
             * and not a radial profile, it would require koliop to modify its
             * interface and the Cpar koliop operator computation to be modified
             * slightly. For now, the uniform coeff_AD_rpeak is applied to
             * produce a koliop compatible radial profile.
             */
            ddc::parallel_fill(
                    Kokkos::DefaultExecutionSpace(),
                    get_field(m_operator_quantities.m_coeff_AD_allocation),
                    coeff_AD);
        }

        Kokkos::fence();
    }

    CollisionConfigurationDataType m_operator_quantities;

    double m_nustar0;
};
```


