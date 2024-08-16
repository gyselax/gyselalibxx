// SPDX-License-Identifier: MIT

#pragma once

#include <paraconf.h>

#include "geometry.hpp"
#include "iinitialization.hpp"
#include "paraconfpp.hpp"
#include "species_info.hpp"

/**
 * @brief A class that initializes the distribution function as a perturbed Maxwellian.
 *
 * A class that initializes the distribution function as a 
 * perturbed Maxwellian defined as $f = f_{maxw}(v) * (1 + perturb(x))$,
 * where $f_{maxw}(v)$ is a Maxwellian, and $perturb(x)$ is a sinusoidal perturbation.
 */
class SingleModePerturbInitialization : public IInitialization
{
    DConstFieldSpVx m_fequilibrium;

    host_t<IFieldMemSp> m_init_perturb_mode;

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
            DFieldX perturbation,
            int const perturb_mode,
            double const perturb_amplitude) const;

    /**
     * @brief Creates an instance of the SingleModePerturbInitialization class.
     * @param[in] fequilibrium A Maxwellian. 
     * @param[in] init_perturb_mode The perturbation mode. 
     * @param[in] init_perturb_amplitude The perturbation amplitude. 
     */
    SingleModePerturbInitialization(
            DConstFieldSpVx fequilibrium,
            host_t<IFieldMemSp> init_perturb_mode,
            host_t<DFieldMemSp> init_perturb_amplitude);

    ~SingleModePerturbInitialization() override = default;

    /**
     * @brief Read init_perturb_mode and init_perturb amplitude in a YAML input file 
     *      to initialize the perturbation. 
     * @param[in] allfequilibrium equilibrium distribution function.
     * @param[in] idx_range_kinsp Index range for the kinetic species.
     * @param[in] yaml_input_file YAML input file.
     * @return an instance of SingleModePerturbInitialization class.
     */
    static SingleModePerturbInitialization init_from_input(
            DConstFieldSpVx allfequilibrium,
            IdxRangeSp idx_range_kinsp,
            PC_tree_t const& yaml_input_file);

    /**
     * @brief Initializes the distribution function as as a perturbed Maxwellian. 
     * @param[in, out] allfdistribu The initialized distribution function.
     * @return The initialized distribution function.
     */
    DFieldSpXVx operator()(DFieldSpXVx allfdistribu) const override;
};
