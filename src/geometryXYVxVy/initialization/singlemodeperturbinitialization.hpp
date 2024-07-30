// SPDX-License-Identifier: MIT

#pragma once

#include <geometry.hpp>
#include <paraconf.h>
#include <species_info.hpp>

#include "iinitialization.hpp"
#include "paraconfpp.hpp"

/// Initialization operator with a sinusoidal perturbation of a Maxwellian. This initializes all species.
class SingleModePerturbInitialization : public IInitialization
{
    DConstFieldSpVxVy m_fequilibrium;

    host_t<IFieldSp> m_init_perturb_mode;

    host_t<DFieldMemSp> m_init_perturb_amplitude;

public:
    /**
     * @brief Initialization of the perturbation.
     * @param[in, out] perturbation On input: an uninitialized array
     *                              On output: an array containing a values that has a 
     *                              sinusoidal variation with given amplitude and mode. 
     * @param[in] perturb_mode The mode of the perturbation. 
     * @param[in] perturb_amplitude The amplitude of the perturbation. 
     */
    void perturbation_initialization(
            DFieldXY perturbation,
            int const perturb_mode,
            double const perturb_amplitude) const;

    /**
     * @brief Creates an instance of the SingleModePerturbInitialization class.
     * @param[in] fequilibrium A Maxwellian. 
     * @param[in] init_perturb_mode The perturbation mode. 
     * @param[in] init_perturb_amplitude The perturbation amplitude. 
     */
    SingleModePerturbInitialization(
            DConstFieldSpVxVy fequilibrium,
            host_t<IFieldSp> init_perturb_mode,
            host_t<DFieldMemSp> init_perturb_amplitude);

    ~SingleModePerturbInitialization() override = default;

    /**
     * @brief Initializes the distribution function as as a perturbed Maxwellian. 
     * @param[in, out] allfdistribu The initialized distribution function.
     * @return The initialized distribution function.
     */
    DFieldSpXYVxVy operator()(DFieldSpXYVxVy allfdistribu) const override;

    /**
     * @brief Read init_perturb_mode and init_perturb amplitude in a YAML input file 
     *      to initialize the perturbation. 
     * @param[in] allfequilibrium equilibrium distribution function.
     * @param[in] dom_kinsp Discrete Domain for the kinetic species.
     * @param[in] yaml_input_file YAML input file.
     * @return an instance of SingleModePerturbInitialization class.
     */
    static SingleModePerturbInitialization init_from_input(
            DConstFieldSpVxVy allfequilibrium,
            IdxRangeSp dom_kinsp,
            PC_tree_t const& yaml_input_file);
};