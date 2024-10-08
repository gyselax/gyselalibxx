// SPDX-License-Identifier: MIT

#include "species_init.hpp"

/**
 * @brief Initialise the species domain
 * @param[out] kinetic species domain 
 * @param[out] fluid species domain 
 * @param[in] conf_voicexx is the YAML input file
 * @param[in] nb_kinspecies number of kinetic species
 * @param[in] nb_fluidspecies number of fluid species
 */
void init_all_species(
        IdxRangeSp& idx_range_kinsp,
        IdxRangeSp& idx_range_fluidsp,
        PC_tree_t conf_voicexx,
        int nb_kinspecies,
        int nb_fluidspecies)
{
    // Define the domain for kinetic species
    idx_range_kinsp = IdxRangeSp(IdxSp(0), IdxStepSp(nb_kinspecies));
    host_t<DFieldMemSp> kinetic_charges(idx_range_kinsp);
    host_t<DFieldMemSp> kinetic_masses(idx_range_kinsp);

    // Define the domain of kinetic species
    for (IdxSp const isp : idx_range_kinsp) {
        PC_tree_t const conf_isp = PCpp_get(conf_voicexx, ".SpeciesInfo[%d]", isp.uid());

        kinetic_charges(isp) = PCpp_double(conf_isp, ".charge");
        kinetic_masses(isp) = PCpp_double(conf_isp, ".mass");
    }

    // Deduce adiabatic species
    int nb_elec_adiabspecies = 1;
    int nb_ion_adiabspecies = 1;
    for (IdxSp const isp : idx_range_kinsp) {
        if (kinetic_charges(isp) == -1.) {
            nb_elec_adiabspecies = 0;
        } else {
            nb_ion_adiabspecies = 0;
        }
    }
    if (nb_elec_adiabspecies + nb_ion_adiabspecies > 1) {
        throw "Error: Not possible to have both adiabatic electrons and ions";
    }

    // Define the domain of fluid species
    idx_range_fluidsp = IdxRangeSp(IdxSp(nb_kinspecies), IdxStepSp(nb_fluidspecies));
    host_t<DFieldMemSp> fluid_charges(idx_range_fluidsp);
    host_t<DFieldMemSp> fluid_masses(idx_range_fluidsp);
    for (IdxSp const isp : idx_range_fluidsp) {
        PC_tree_t const conf_nisp = PCpp_get(
                conf_voicexx,
                ".NeutralSpeciesInfo[%d]",
                isp.uid() - idx_range_fluidsp.front().uid());
        fluid_charges(isp) = 0.; // neutral charge is zero
        fluid_masses(isp) = PCpp_double(conf_nisp, ".mass");
    }

    // Create the domain of all species including kinetic species + fluid species + adiabatic species (if existing)
    // adiabatic species are placed at the back of the domain
    IdxRangeSp const idx_range_allsp(
            IdxSp(0),
            IdxStepSp(
                    nb_kinspecies + nb_fluidspecies + nb_elec_adiabspecies + nb_ion_adiabspecies));
    host_t<DFieldMemSp> charges(idx_range_allsp);
    for (IdxSp isp : idx_range_kinsp) {
        charges(isp) = kinetic_charges(isp);
    }
    for (IdxSp isp : idx_range_fluidsp) {
        charges(isp) = fluid_charges(isp);
    }
    if (nb_elec_adiabspecies + nb_ion_adiabspecies > 0) {
        charges(idx_range_allsp.back()) = nb_ion_adiabspecies - nb_elec_adiabspecies;
    }

    // Create the domain of kinetic and fluid species
    IdxRangeSp const idx_range_kinfluidsp(IdxSp(0), IdxStepSp(nb_kinspecies + nb_fluidspecies));
    host_t<DFieldMemSp> masses(idx_range_kinfluidsp);
    for (IdxSp isp : idx_range_kinsp) {
        masses(isp) = kinetic_masses(isp);
    }
    for (IdxSp isp : idx_range_fluidsp) {
        masses(isp) = fluid_masses(isp);
    }

    ddc::init_discrete_space<Species>(std::move(charges), std::move(masses));
}


/**
 * @brief Initialise the species domain in the case of adiabatic and kinetic species
 * @param[in] conf_voicexx is the YAML input file
 * @return the kinetic species domain 
 */
IdxRangeSp init_species(PC_tree_t conf_voicexx)
{
    IdxStepSp const nb_kinspecies(PCpp_len(conf_voicexx, ".SpeciesInfo"));
    IdxStepSp const nb_fluidspecies(0);
    IdxRangeSp idx_range_kinsp;
    IdxRangeSp idx_range_fluidsp;
    init_all_species(
            idx_range_kinsp,
            idx_range_fluidsp,
            conf_voicexx,
            nb_kinspecies,
            nb_fluidspecies);
    return idx_range_kinsp;
}


/**
 * @brief Initialise the species domain in the specific case of fluid species 
 *   added to adiabatic and kinetic species
 * @param[out] kinetic species domain 
 * @param[out] fluid species domain 
 * @param[in] conf_voicexx is the YAML input file
 */
void init_species_withfluid(
        IdxRangeSp& idx_range_kinsp,
        IdxRangeSp& idx_range_fluidsp,
        PC_tree_t conf_voicexx)
{
    IdxStepSp const nb_kinspecies(PCpp_len(conf_voicexx, ".SpeciesInfo"));
    IdxStepSp const nb_fluidspecies(PCpp_len(conf_voicexx, ".NeutralSpeciesInfo"));
    init_all_species(
            idx_range_kinsp,
            idx_range_fluidsp,
            conf_voicexx,
            nb_kinspecies,
            nb_fluidspecies);
}
