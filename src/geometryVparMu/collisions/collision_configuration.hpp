// SPDX-License-Identifier: MIT

#pragma once

#include <paraconf.h>

#include "collision_common_configuration.hpp"
#include "paraconfpp.hpp"

/**
 * @brief Class to collect information to initialise the collision operator for
 * a SpVparMu geometry.
 *
 * NOTE: Thanks to C++17 template constructor argument type deduction, we do not
 * need to give the 3 template arguments. Giving a correct
 * index_range_fdistribution is sufficient.
 */
template <class GridSp, class GridVpar, class GridMu>
class CollisionConfiguration
{
public:
    /**
     * @brief Mu grid.
     */
    using GridMuType = GridMu;
    /**
     * @brief Vpar grid.
     */
    using GridVparType = GridVpar;
    /**
     * @brief Theta grid.
     */
    using GridThetaType = detail::InternalSpoofGridTheta;
    /**
     * @brief R grid.
     */
    using GridRType = detail::InternalSpoofGridR;
    /**
     * @brief Phi grid.
     */
    using GridPhiType = detail::InternalSpoofGridPhi;
    /**
     * @brief Sp grid.
     */
    using GridSpType = GridSp;

    /**
     * @brief Mu index range.
     */
    using IdxRangeMuType = IdxRange<GridMuType>;
    /**
     * @brief Vpar index range.
     */
    using IdxRangeVparType = IdxRange<GridVparType>;
    /**
     * @brief Sp,Vpar index range.
     */
    using IdxRangeSpVparType = IdxRange<GridSpType, GridVparType>;
    /**
     * @brief Distribution function index range.
     */
    using IdxRangeDistributionFunctionType = IdxRange<GridSpType, GridVparType, GridMuType>;

    /**
     * @brief Container for the operator input data.
     */
    using CollisionConfigurationDataType = detail::CollisionConfigurationData<
            IdxRangeDistributionFunctionType,
            GridSpType,
            GridPhiType,
            GridRType,
            GridThetaType,
            GridVparType,
            GridMuType>;

public:
    /**
     * @brief The constructor for the CollisionConfiguration class.
     *
     * @param yaml_input_file yaml namelist file object.
     * @param index_range_fdistribution index range of the whole distribution
     * function.
     * @param coeff_intdmu quadrature coefficients.
     * @param coeff_intdvpar quadrature coefficients.
     * @param[in] B_norm The norm of the magnetic field defined over the polar
     * slice.
     * @param[in] Bstar_s Bstar defined on the coordinates
     * (species,r,theta,vpar).
     */
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

    /**
     * @brief Can be used to obtain the quantities the operator needs for its
     * initialisation, modification is not authorised.
     *
     * @return const CollisionConfigurationDataType&
     */
    const CollisionConfigurationDataType& configuration() const
    {
        return m_operator_quantities;
    }

protected:
    /**
     * @brief Does the initialisation specific to this geometry.
     *
     * @param[in] input_B_norm
     * @param[in] input_Bstar_s
     */
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

    /**
     * @brief Container for the arrays the quantities the operator needs.
     */
    CollisionConfigurationDataType m_operator_quantities;

    /**
     * @brief Value of nustar0 = to nustar0_rpeak read in the YAML input file.
     */
    double m_nustar0;
};
