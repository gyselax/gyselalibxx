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
        ddc::DiscreteDomain<IDimSp>& dom_kinsp,
        ddc::DiscreteDomain<IDimSp>& dom_fluidsp,
        PC_tree_t conf_voicexx,
        int nb_kinspecies,
        int nb_fluidspecies)
{
    using IndexSp = ddc::DiscreteElement<IDimSp>;
    using IDomainSp = ddc::DiscreteDomain<IDimSp>;
    using IVectSp = ddc::DiscreteVector<IDimSp>;
    using DFieldSp = device_t<ddc::Chunk<double, IDomainSp>>;
    using IFieldSp = device_t<ddc::Chunk<int, IDomainSp>>;

    // Define the domain for kinetic species
    dom_kinsp = IDomainSp(IndexSp(0), IVectSp(nb_kinspecies));
    host_t<IFieldSp> kinetic_charges(dom_kinsp);
    host_t<DFieldSp> kinetic_masses(dom_kinsp);

    // Define the domain of kinetic species
    for (IndexSp const isp : dom_kinsp) {
        PC_tree_t const conf_isp = PCpp_get(conf_voicexx, ".SpeciesInfo[%d]", isp.uid());

        kinetic_charges(isp) = static_cast<int>(PCpp_int(conf_isp, ".charge"));
        kinetic_masses(isp) = PCpp_double(conf_isp, ".mass");
    }

    // Deduce adiabatic species
    int nb_elec_adiabspecies = 1;
    int nb_ion_adiabspecies = 1;
    for (IndexSp const isp : dom_kinsp) {
        if (kinetic_charges(isp) == -1) {
            nb_elec_adiabspecies = 0;
        } else {
            nb_ion_adiabspecies = 0;
        }
    }
    if (nb_elec_adiabspecies + nb_ion_adiabspecies > 1) {
        throw "Error: Not possible to have both adiabatic electrons and ions";
    }

    // Define the domain of fluid species
    dom_fluidsp = IDomainSp(IndexSp(nb_kinspecies), IVectSp(nb_fluidspecies));
    host_t<IFieldSp> fluid_charges(dom_fluidsp);
    host_t<DFieldSp> fluid_masses(dom_kinsp);
    for (IndexSp const isp : dom_fluidsp) {
        PC_tree_t const conf_nisp = PCpp_get(
                conf_voicexx,
                ".NeutralSpeciesInfo[%d]",
                isp.uid() - dom_fluidsp.front().uid());
        fluid_charges(isp) = 0; // neutral charge is zero
        fluid_masses(isp) = PCpp_double(conf_nisp, ".mass");
    }

    // Create the domain of all species including kinetic species + fluid species + adiabatic species (if existing)
    // adiabatic species are placed at the back of the domain
    IDomainSp const dom_allsp(
            IndexSp(0),
            IVectSp(nb_kinspecies + nb_fluidspecies + nb_elec_adiabspecies + nb_ion_adiabspecies));
    host_t<IFieldSp> charges(dom_allsp);
    for (IndexSp isp : dom_kinsp) {
        charges(isp) = kinetic_charges(isp);
    }
    for (IndexSp isp : dom_fluidsp) {
        charges(isp) = fluid_charges(isp);
    }
    if (nb_elec_adiabspecies + nb_ion_adiabspecies > 0) {
        charges(dom_allsp.back()) = nb_ion_adiabspecies - nb_elec_adiabspecies;
    }

    // Create the domain of kinetic and fluid species
    IDomainSp const dom_kinfluidsp(IndexSp(0), IVectSp(nb_kinspecies + nb_fluidspecies));
    host_t<DFieldSp> masses(dom_kinfluidsp);
    for (IndexSp isp : dom_kinsp) {
        masses(isp) = kinetic_masses(isp);
    }
    for (IndexSp isp : dom_fluidsp) {
        masses(isp) = fluid_masses(isp);
    }

    ddc::init_discrete_space<IDimSp>(std::move(charges), std::move(masses));
}


/**
 * @brief Initialise the species domain in the case of adiabatic and kinetic species
 * @param[in] conf_voicexx is the YAML input file
 * @return the kinetic species domain 
 */
ddc::DiscreteDomain<IDimSp> init_species(PC_tree_t conf_voicexx)
{
    using IVectSp = ddc::DiscreteVector<IDimSp>;

    IVectSp const nb_kinspecies(PCpp_len(conf_voicexx, ".SpeciesInfo"));
    IVectSp const nb_fluidspecies(0);
    ddc::DiscreteDomain<IDimSp> dom_kinsp;
    ddc::DiscreteDomain<IDimSp> dom_fluidsp;
    init_all_species(dom_kinsp, dom_fluidsp, conf_voicexx, nb_kinspecies, nb_fluidspecies);
    return dom_kinsp;
}


/**
 * @brief Initialise the species domain in the specific case of fluid species 
 *   added to adiabatic and kinetic species
 * @param[out] kinetic species domain 
 * @param[out] fluid species domain 
 * @param[in] conf_voicexx is the YAML input file
 */
void init_species_withfluid(
        ddc::DiscreteDomain<IDimSp>& dom_kinsp,
        ddc::DiscreteDomain<IDimSp>& dom_fluidsp,
        PC_tree_t conf_voicexx)
{
    using IVectSp = ddc::DiscreteVector<IDimSp>;

    IVectSp const nb_kinspecies(PCpp_len(conf_voicexx, ".SpeciesInfo"));
    IVectSp const nb_fluidspecies(PCpp_len(conf_voicexx, ".NeutralSpeciesInfo"));
    init_all_species(dom_kinsp, dom_fluidsp, conf_voicexx, nb_kinspecies, nb_fluidspecies);
}
