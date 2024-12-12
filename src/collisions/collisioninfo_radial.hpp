// SPDX-License-Identifier: MIT

#pragma once
#include <paraconf.h>

#include "ddc_alias_inline_functions.hpp"
#include "ddc_aliases.hpp"
#include "paraconfpp.hpp"
#include "species_info.hpp"

/// Class to collect information to initialize the collision operator
template <class GridR>
class CollisionInfoRadial
{
private:
    using IdxR = Idx<GridR>;
    using IdxRangeR = IdxRange<GridR>;

public:
    /// Type alias for a field memory block on a grid of radial values
    using DFieldMemR = DFieldMem<IdxRangeR>;
    /// Type alias for a field defined on a grid of radial values
    using DFieldR = DField<IdxRangeR>;
    /// Type alias for a field defined on a grid of radial values
    using DConstFieldR = DConstField<IdxRangeR>;
    /// radial_chunk_type used to treat the 0D case for radial profile
    using radial_chunk_type = DConstFieldR;

private:
    /// value of nustar0_rpeak read in the YAML input file
    double m_nustar0_rpeak;

    /// boolean that is equal to true if inter-species collisions are taken into account
    // read in the YAML input file
    std::int8_t m_collisions_interspecies;

    // [TODO] Values that are temporarily need for the interface with koliop
    // but that will be deleted as soon as no more required
    double const m_rpeak;
    double const m_q_rpeak;
    double const m_R0;

    DConstFieldR m_rg; //VG ?
    DConstFieldR m_safety_factor; //VG?

    /// radial coefficients of AD
    DFieldMemR m_coeff_AD;

public:
    /**
     * Computation of the radial profile of AD
     *  @f[ AD(rpeak) = \sqrt(2)*eps_rpeak^(3/2)/(q_rpeak R0)*nustar0_rpeak @f]
     * where R0, q_rpeak, nustar0_rpeak are given as input data and
     * eps_rpeak = eps(rpeak) with eps(r) = r/R0 is the aspect ratio.
     * 
     * @param[in] idxrange_r radial index range
     * 
     * This function should be private but is public due to Kokkos restrictions.
    */
    void compute_coeff_AD(IdxRangeR idxrange_r)
    {
        double const rpeak = m_rpeak;
        double const q_rpeak = m_q_rpeak;
        double const R0 = m_R0;
        double const nustar0_rpeak = m_nustar0_rpeak;

        // Compute coeff_AD_rpeak
        DFieldR coeff_AD = get_field(m_coeff_AD);
        double const eps_rpeak = rpeak / R0;
        double const coeff_AD_rpeak
                = std::sqrt(2.) * eps_rpeak * std::sqrt(eps_rpeak) * nustar0_rpeak / (q_rpeak * R0);

        // [TODO] coeff_AD should be a scalar and not a radial profile
        //  but to change that koliop interface must be changed
        //  so coeff_AD(r) is fixed equal to coeff_AD_rpeak
        ddc::parallel_for_each(
                Kokkos::DefaultExecutionSpace(),
                idxrange_r,
                KOKKOS_LAMBDA(Idx<GridR> idx_r) { coeff_AD(idx_r) = coeff_AD_rpeak; });
    }

public:
    /**
     * @brief The constructor for the CollisionInfoRadial class.
     * @param[in] rpeak the radial reference position r = rpeak
     * @param[in] q_rpeak the safety factor value at r = rpeak
     * @param[in] R0 the major radius at the magnetic axis
     * @param[in] nustar0_rpeak the value of nustar at r = rpeak
     * @param[in] collisions_interspecies logical value equal to true if interspecies collisions are taken into account
     * @param[in] idxrange_r radial index range
     */
    CollisionInfoRadial(
            double const rpeak,
            double const q_rpeak,
            double const R0,
            double const nustar0_rpeak,
            std::int8_t const collisions_interspecies,
            IdxRangeR idxrange_r)
        : m_rpeak(rpeak)
        , m_q_rpeak(q_rpeak)
        , m_R0(R0)
        , m_nustar0_rpeak(nustar0_rpeak)
        , m_collisions_interspecies(collisions_interspecies)
        , m_coeff_AD(idxrange_r)
    {
        compute_coeff_AD(idxrange_r);
    };

    ~CollisionInfoRadial() = default;


    /**
     * @brief A method for accessing the collisions_interspecies member variable of the class.
     * @return A boolean containing the collisions_interspecies value. 
     */
    std::int8_t collisions_interspecies() const
    {
        return m_collisions_interspecies;
    }


    /**
     * @brief A method for accessing coeff_AD (the radial profile of AD coeff) variable of the class.
     * @return A field containing the radial profile of the AD coefficients.
     */
    DConstFieldR coeff_AD() const
    {
        return get_const_field(m_coeff_AD);
    }
};
