// SPDX-License-Identifier: MIT

#pragma once

#include <geometry.hpp>
#include <species_info.hpp>

#include "iinitialization.hpp"

/// Initialization operator with a sinusoidal perturbation of a Maxwellian. This initializes all species.
class SingleModePerturbInitialization : public IInitialization
{
    DViewSpVxVy m_fequilibrium;

    host_t<ViewSp<int>> m_init_perturb_mode;

    host_t<DViewSp> m_init_perturb_amplitude;

public:
    /**
     * @brief Initialization of the perturbation.
     * @param[in, out] perturbation On input: an uninitialized array
     *                              On output: an array containing a values that has a 
     *                              sinusoidal variation with given amplitude and mode. 
     * @param[in] mode The mode of the perturbation. 
     * @param[in] perturb_amplitude The amplitude of the perturbation. 
     */
    void perturbation_initialization(DSpanXY perturbation, int mode, double perturb_amplitude)
            const;

    /**
     * @brief Creates an instance of the SingleModePerturbInitialization class.
     * @param[in] fequilibrium A Maxwellian. 
     * @param[in] init_perturb_mode The perturbation mode. 
     * @param[in] init_perturb_amplitude The perturbation amplitude. 
     */
    SingleModePerturbInitialization(
            DViewSpVxVy fequilibrium,
            host_t<ViewSp<int>> init_perturb_mode,
            host_t<DViewSp> init_perturb_amplitude);

    ~SingleModePerturbInitialization() override = default;

    /**
     * @brief Initializes the distribution function as as a perturbed Maxwellian. 
     * @param[in, out] allfdistribu The initialized distribution function.
     * @return The initialized distribution function.
     */
    DSpanSpXYVxVy operator()(DSpanSpXYVxVy allfdistribu) const override;
};
