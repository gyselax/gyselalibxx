// SPDX-License-Identifier: MIT

#include "species_init.hpp"

void init_all_species(
        IdxRangeSp& idxrange_kinsp,
        IdxRangeSp& idxrange_fluidsp,
        PC_tree_t conf_gyselalibxx,
        int nb_kinspecies,
        int nb_fluidspecies)
{
    // Define the domain for kinetic species
    idxrange_kinsp = IdxRangeSp(IdxSp(0), IdxStepSp(nb_kinspecies));
    host_t<DFieldMemSp> kinetic_charges(idxrange_kinsp);
    host_t<DFieldMemSp> kinetic_masses(idxrange_kinsp);

    // Define the domain of kinetic species
    for (IdxSp const isp : idxrange_kinsp) {
        PC_tree_t const conf_isp = PCpp_get(
                conf_gyselalibxx,
                ".SpeciesInfo[%d]",
                (isp - idxrange_kinsp.front()).value());

        kinetic_charges(isp) = PCpp_double(conf_isp, ".charge");
        kinetic_masses(isp) = PCpp_double(conf_isp, ".mass");
    }

    // Deduce adiabatic species
    int nb_elec_adiabspecies = 1;
    int nb_ion_adiabspecies = 1;
    for (IdxSp const isp : idxrange_kinsp) {
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
    idxrange_fluidsp = IdxRangeSp(IdxSp(nb_kinspecies), IdxStepSp(nb_fluidspecies));
    host_t<DFieldMemSp> fluid_charges(idxrange_fluidsp);
    host_t<DFieldMemSp> fluid_masses(idxrange_fluidsp);
    for (IdxSp const isp : idxrange_fluidsp) {
        PC_tree_t const conf_nisp = PCpp_get(
                conf_gyselalibxx,
                ".NeutralSpeciesInfo[%d]",
                (isp - idxrange_fluidsp.front()).value());
        fluid_charges(isp) = 0.; // neutral charge is zero
        fluid_masses(isp) = PCpp_double(conf_nisp, ".mass");
    }

    // Create the domain of all species including kinetic species + fluid species + adiabatic species (if existing)
    // adiabatic species are placed at the back of the domain
    IdxRangeSp const idxrange_allsp(
            IdxSp(0),
            IdxStepSp(
                    nb_kinspecies + nb_fluidspecies + nb_elec_adiabspecies + nb_ion_adiabspecies));
    host_t<DFieldMemSp> charges(idxrange_allsp);
    for (IdxSp isp : idxrange_kinsp) {
        charges(isp) = kinetic_charges(isp);
    }
    for (IdxSp isp : idxrange_fluidsp) {
        charges(isp) = fluid_charges(isp);
    }
    if (nb_elec_adiabspecies + nb_ion_adiabspecies > 0) {
        charges(idxrange_allsp.back()) = nb_ion_adiabspecies - nb_elec_adiabspecies;
    }

    // Create the domain of kinetic and fluid species
    IdxRangeSp const idxrange_kinfluidsp(IdxSp(0), IdxStepSp(nb_kinspecies + nb_fluidspecies));
    host_t<DFieldMemSp> masses(idxrange_kinfluidsp);
    for (IdxSp isp : idxrange_kinsp) {
        masses(isp) = kinetic_masses(isp);
    }
    for (IdxSp isp : idxrange_fluidsp) {
        masses(isp) = fluid_masses(isp);
    }

    ddc::init_discrete_space<Species>(std::move(charges), std::move(masses));
}


/**
 * @brief Initialise the species domain in the case of adiabatic and kinetic species
 * @param[in] conf_gyselalibxx is the YAML input file
 * @return the kinetic species domain 
 */
IdxRangeSp init_species(PC_tree_t conf_gyselalibxx)
{
    IdxStepSp const nb_kinspecies(PCpp_len(conf_gyselalibxx, ".SpeciesInfo"));
    IdxStepSp const nb_fluidspecies(0);
    IdxRangeSp idxrange_kinsp;
    IdxRangeSp idxrange_fluidsp;
    init_all_species(
            idxrange_kinsp,
            idxrange_fluidsp,
            conf_gyselalibxx,
            nb_kinspecies,
            nb_fluidspecies);
    return idxrange_kinsp;
}


/**
 * @brief Initialise the species domain in the specific case of fluid species 
 *   added to adiabatic and kinetic species
 * @param[out] kinetic species domain 
 * @param[out] fluid species domain 
 * @param[in] conf_gyselalibxx is the YAML input file
 */
void init_species_withfluid(
        IdxRangeSp& idxrange_kinsp,
        IdxRangeSp& idxrange_fluidsp,
        PC_tree_t conf_gyselalibxx)
{
    IdxStepSp const nb_kinspecies(PCpp_len(conf_gyselalibxx, ".SpeciesInfo"));
    IdxStepSp const nb_fluidspecies(PCpp_len(conf_gyselalibxx, ".NeutralSpeciesInfo"));
    init_all_species(
            idxrange_kinsp,
            idxrange_fluidsp,
            conf_gyselalibxx,
            nb_kinspecies,
            nb_fluidspecies);
}


/**
 * @brief Initialise the kinetic species domain 
 * @return the kinetic species domain 
 */
IdxRangeSp init_kinetic_species()
{
    std::vector<int> species;
    std::vector<double> charges;
    std::vector<double> masses;
    PDI_get_arrays("read_species", "species", species, "charges", charges, "masses", masses);
    IdxStepSp const idxstep_kinsp(charges.size());
    IdxRangeSp const idxrange_kinsp(IdxSp(0), idxstep_kinsp);
    host_t<DFieldMemSp> kinetic_charges(idxrange_kinsp);
    host_t<DFieldMemSp> kinetic_masses(idxrange_kinsp);
    // Define the domain of kinetic species
    for (IdxSp const isp : idxrange_kinsp) {
        kinetic_charges(isp) = charges[isp - idxrange_kinsp.front()];
        kinetic_masses(isp) = masses[isp - idxrange_kinsp.front()];
    }
    ddc::init_discrete_space<Species>(std::move(kinetic_charges), std::move(kinetic_masses));

    return idxrange_kinsp;
}
