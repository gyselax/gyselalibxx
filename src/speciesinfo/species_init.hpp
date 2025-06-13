// SPDX-License-Identifier: MIT

#pragma once

#include <ddc/ddc.hpp>

#include <paraconf.h>

#include "ddc_helper.hpp"
#include "paraconfpp.hpp"
#include "pdi_helper.hpp"
#include "species_info.hpp"

/**
 * @brief Initialise the species domain
 * @param[out] idx_range_kinsp kinetic species domain 
 * @param[out] idx_range_fluidsp fluid species domain 
 * @param[in] conf_gyselalibxx is the YAML input file
 * @param[in] nb_kinspecies number of kinetic species
 * @param[in] nb_fluidspecies number of fluid species
 */
void init_all_species(
        IdxRangeSp& idx_range_kinsp,
        IdxRangeSp& idx_range_fluidsp,
        PC_tree_t conf_gyselalibxx,
        int nb_kinspecies,
        int nb_fluidspecies);

/**
 * @brief Initialise the species domain in the case of adiabatic and kinetic species
 * @param[in] conf_gyselalibxx is the YAML input file
 * @return the kinetic species domain 
 */
IdxRangeSp init_species(PC_tree_t conf_gyselalibxx);

/**
 * @brief Initialise the species domain in the specific case of fluid species 
 *   added to adiabatic and kinetic species
 * @param[out] idx_range_kinsp kinetic species domain 
 * @param[out] idx_range_fluidsp fluid species domain 
 * @param[in] conf_gyselalibxx is the YAML input file
 */
void init_species_withfluid(
        IdxRangeSp& idx_range_kinsp,
        IdxRangeSp& idx_range_fluidsp,
        PC_tree_t conf_gyselalibxx);

/**
 * @brief Initialise the kinetic species domain 
 * @return the kinetic species domain 
 */
IdxRangeSp init_kinetic_species();
