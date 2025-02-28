// SPDX-License-Identifier: MIT

#pragma once

#include <paraconf.h>

#include "geometry.hpp"
#include "iinitialisation.hpp"
#include "paraconfpp.hpp"
#include "species_info.hpp"

/// Initialisation operator with a sinusoidal perturbation of a Maxwellian. This initialises all species.
class SingleModePerturbInitialisation : public IInitialisation
{
    DConstFieldSpVxVy m_fequilibrium;

    host_t<IFieldMemSp> m_init_perturb_mode;

    host_t<DFieldMemSp> m_init_perturb_amplitude;

public:
    /**
     * @brief Initialisation of the perturbation.
     * @param[in, out] perturbation On input: an uninitialized array
     *                              On output: an array containing a values that has a 
     *                              sinusoidal variation with given amplitude and mode. 
     * @param[in] perturb_mode The mode of the perturbation. 
     * @param[in] perturb_amplitude The amplitude of the perturbation. 
     */
    void perturbation_initialisation(
            DFieldXY perturbation,
            int const perturb_mode,
            double const perturb_amplitude) const;

    /**
     * @brief Creates an instance of the SingleModePerturbInitialisation class.
     * @param[in] fequilibrium A Maxwellian. 
     * @param[in] init_perturb_mode The perturbation mode. 
     * @param[in] init_perturb_amplitude The perturbation amplitude. 
     */
    SingleModePerturbInitialisation(
            DConstFieldSpVxVy fequilibrium,
            host_t<IFieldMemSp> init_perturb_mode,
            host_t<DFieldMemSp> init_perturb_amplitude);

    ~SingleModePerturbInitialisation() override = default;

    /**
     * @brief Initialises the distribution function as as a perturbed Maxwellian. 
     * @param[in, out] allfdistribu The initialised distribution function.
     * @return The initialised distribution function.
     */
    DFieldSpXYVxVy operator()(DFieldSpXYVxVy allfdistribu) const override;

    /**
     * @brief Read init_perturb_mode and init_perturb amplitude in a YAML input file 
     *      to initialise the perturbation. 
     * @param[in] allfequilibrium equilibrium distribution function.
     * @param[in] idx_range_kinsp Index range for the kinetic species.
     * @param[in] yaml_input_file YAML input file.
     * @return an instance of SingleModePerturbInitialisation class.
     */
    static SingleModePerturbInitialisation init_from_input(
            DConstFieldSpVxVy allfequilibrium,
            IdxRangeSp idx_range_kinsp,
            PC_tree_t const& yaml_input_file);
};
