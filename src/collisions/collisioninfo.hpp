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

    // AD coefficient almost equivalent to the collision frequency
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

        /**
         * Computation of the radial profile of ADi that is required for koliop interface
         *  @f[ AD(rpeak) = \sqrt(2)*eps_rpeak^(3/2)/(q_rpeak R0)*nustar0_rpeak @f]
         * eps_rpeak = eps(rpeak) with eps(r) = r/R0 is the aspect ratio.
         * 
         * In this specific case, where AD is a scalar all values are forced to 1, so 
         *   coeff_AD_atrpeak = sqrt(2)
         **/
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
     * @brief A method for accessing the coeff_AD variable of the class.
     * @return A double containing the coeff_AD value. 
     */
    double coeff_AD() const
    {
        return m_coeff_AD;
    }
};
