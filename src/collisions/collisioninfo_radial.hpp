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
    DConstFieldR m_rg;
    DConstFieldR m_safety_factor;

    /// radial profile of nustar0
    DFieldMemR m_nustar0_r;

public:
    /**
     * Computation of the radial profile of nustar0 as
     *  @f[ nustar0(r) = nustar0(rp) * q(r)/q(rp)*(rp/r)^{3/2} @f]
     * where rp is equal to rpeak and where nustar0(rp) is given as input data.
     * 
     * This function should be private but is public due to Kokkos restrictions.
    */
    void compute_nustar0_r()
    {
        double const rpeak = m_rpeak;
        double const q_rpeak = m_q_rpeak;
        DConstFieldR radial_profile = get_const_field(m_rg);
        DConstFieldR safety_factor = get_const_field(m_safety_factor);
        double const nustar0_rpeak = m_nustar0_rpeak;

        DFieldR nustar0_r = get_field(m_nustar0_r);
        ddc::parallel_for_each(
                Kokkos::DefaultExecutionSpace(),
                get_idx_range(radial_profile),
                KOKKOS_LAMBDA(Idx<GridR> idx) {
                    double rpeak_on_r = rpeak / radial_profile(idx);
                    double q_on_qrpeak = safety_factor(idx) / q_rpeak;
                    nustar0_r(idx)
                            = nustar0_rpeak * q_on_qrpeak * rpeak_on_r * std::sqrt(rpeak_on_r);
                });
    }

public:
    /**
     * @brief The constructor for the CollisionInfoRadial class.
     * @param[in] yaml_input_file YAML input file containing CollisionsInfo
     * @param[in] rpeak radial position of the gradient
     * @param[in] q_rpeak safety factor value at r=rpeak
     * @param[in] radial_profile value of the grid points in the radial direction
     * @param[in] safety_factor radial profile of the safety factor
     */
    CollisionInfoRadial(
            PC_tree_t const& yaml_input_file,
            double const rpeak,
            double const q_rpeak,
            DConstFieldR radial_profile,
            DConstFieldR safety_factor)
        : m_nustar0_rpeak {PCpp_double(yaml_input_file, ".CollisionsInfo.nustar0_rpeak")}
        , m_collisions_interspecies {PCpp_bool(
                  yaml_input_file,
                  ".CollisionsInfo.collisions_interspecies")}
        , m_rpeak(rpeak)
        , m_q_rpeak(q_rpeak)
        , m_rg(radial_profile)
        , m_safety_factor(safety_factor)
        , m_nustar0_r(get_idx_range(radial_profile))
    {
        compute_nustar0_r();
    };

    /**
     * @brief The constructor for the CollisionInfoRadial class.
     * @param[in] nustar0_rpeak value of nustar0 at r=rpeak
     * @param[in] collisions_interspecies
     *     boolean that is equal to true if inter-species collisions are taken into account
     * @param[in] rpeak radial position of the gradient
     * @param[in] q_rpeak safety factor value at r=rpeak
     * @param[in] radial_profile value of the grid points in the radial direction
     * @param[in] safety_factor radial profile of the safety factor
     */
    CollisionInfoRadial(
            double const nustar0_rpeak,
            std::int8_t const collisions_interspecies,
            double const rpeak,
            double const q_rpeak,
            DConstFieldR radial_profile,
            DConstFieldR safety_factor)
        : m_nustar0_rpeak(nustar0_rpeak)
        , m_collisions_interspecies(collisions_interspecies)
        , m_rpeak(rpeak)
        , m_q_rpeak(q_rpeak)
        , m_rg(radial_profile)
        , m_safety_factor(safety_factor)
        , m_nustar0_r(get_idx_range(radial_profile))
    {
        compute_nustar0_r();
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
     * @brief A method for accessing the radial profile variable of the class.
     * @return A field containing the radial profile of the grid. 
     */
    DConstFieldR rg() const
    {
        return m_rg;
    }

    /**
     * @brief A method for accessing the safety factor variable of the class.
     * @return A field containing the radial profile of the safety_factor. 
     */
    DConstFieldR safety_factor() const
    {
        return m_safety_factor;
    }

    /**
     * @brief A method for accessing nustar0_r (the radial profile of nustar0) variable of the class.
     * @return A field containing the radial profile of nustar0.
     */
    DConstFieldR nustar0() const
    {
        return m_nustar0_r;
    }
};
