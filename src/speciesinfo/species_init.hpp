// SPDX-License-Identifier: MIT

#pragma once

#include <ddc/ddc.hpp>

#include <paraconf.h>

#include "ddc_helper.hpp"
#include "paraconfpp.hpp"
#include "species_info.hpp"

/**
 * @brief Initialise the species domain
 * @param[out] kinetic species domain 
 * @param[out] fluid species domain 
 * @param[in] conf_voicexx is the YAML input file
 * @param[in] nb_kinspecies number of kinetic species
 * @param[in] nb_fluidspecies number of fluid species
 */
void init_all_species(
        IdxRangeSp& dom_kinsp,
        IdxRangeSp& dom_fluidsp,
        PC_tree_t conf_voicexx,
        int nb_kinspecies,
        int nb_fluidspecies);

/**
 * @brief Initialise the species domain in the case of adiabatic and kinetic species
 * @param[in] conf_voicexx is the YAML input file
 * @return the kinetic species domain 
 */
IdxRangeSp init_species(PC_tree_t conf_voicexx);

/**
 * @brief Initialise the species domain in the specific case of fluid species 
 *   added to adiabatic and kinetic species
 * @param[out] kinetic species domain 
 * @param[out] fluid species domain 
 * @param[in] conf_voicexx is the YAML input file
 */
void init_species_withfluid(IdxRangeSp& dom_kinsp, IdxRangeSp& dom_fluidsp, PC_tree_t conf_voicexx);
