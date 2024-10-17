// SPDX-License-Identifier: MIT

#pragma once

#include <paraconf.h>

#include "paraconfpp.hpp"
#include "species_info.hpp"

/// Class to collect information to initialize the collision operator
class CollisionInfo
{
    // value of nustar0 = to nustar0_rpeak read in the YAML input file
    double const m_nustar0;

    // boolean that is equal to true if inter-species collisions are taken into account
    // read in the YAML input file
    std::int8_t m_collisions_interspecies;

    // [TODO] Values that are temporarily needed for the interface with koliop
    // but that will be deleted as soon as they are no longer required
    // rg fixed to 1. => because already included in nustar
    double const m_rpeak = 1.;
    double const m_rg = 1.;
    // safety factor fixed to 1. => because already included in nustar
    double const m_q_rpeak = 1.;
    double const m_safety_factor = 1.;
    // AD coefficient
    double m_coeff_AD;

public:
    /// radial_chunk_type used to treat the 0D case for radial profile
    using radial_chunk_type = double;

    /**
     * @brief The constructor for the CollisionFrequency class.
     * @param[in] yaml_input_file YAML input file containing CollisionsInfo
     */
    explicit CollisionInfo(PC_tree_t const& yaml_input_file)
        : m_nustar0 {PCpp_double(yaml_input_file, ".CollisionsInfo.nustar0_rpeak")}
        , m_collisions_interspecies {
                  PCpp_bool(yaml_input_file, ".CollisionsInfo.collisions_interspecies")}
    {
        // nustar fixed to 1. => normalisation of time is equal to 1./nustar
        if (m_nustar0 != 1.0) {
            throw std::invalid_argument("nustar0 must be equal to 1");
        }

        m_coeff_AD = Kokkos::sqrt(2.);
    };

    ~CollisionInfo() = default;

    /**
     * @brief A method for accessing the collisions_interspecies member variable of the class.
     * @return A boolean containing the collisions_interspecies value. 
     */
    std::int8_t collisions_interspecies() const
    {
        return m_collisions_interspecies;
    }

    /**
     * @brief A method for accessing the rg variable of the class.
     * @return A double containing the rg value. 
     */
    double rg() const
    {
        return m_rg;
    }

    /**
     * @brief A method for accessing the safety factor variable of the class.
     * @return A double containing the safety factor value. 
     */
    double safety_factor() const
    {
        return m_safety_factor;
    }

    /**
     * @brief A method for accessing the nustar0 variable of the class.
     * @return A double containing the nustar0 value. 
     */
    double nustar0() const
    {
        return m_nustar0;
    }

    /**
     * @brief A method for accessing the coeff_AD variable of the class.
     * @return A double containing the coeff_AD value. 
     */
    double coeff_AD() const
    {
        return m_coeff_AD;
    }
};
