// SPDX-License-Identifier: MIT

#pragma once

#include <geometry.hpp>
#include <species_info.hpp>

#include "iinitialization.hpp"

/**
 * @brief A class that initializes the distribution function as a perturbed Maxwellian.
 *
 * A class that initializes the distribution function as a 
 * perturbed Maxwellian defined as $f = f_{maxw}(v) * (1 + perturb(x))$,
 * where $f_{maxw}(v)$ is a Maxwellian, and $perturb(x)$ is a sinusoidal perturbation.
 */
class SingleModePerturbInitialization : public IInitialization
{
    device_t<DViewSpVx> m_fequilibrium;

    ViewSp<int> m_init_perturb_mode;

    DViewSp m_init_perturb_amplitude;

public:
    /**
     * @brief Initialization of the perturbation.
     * @param[in, out] perturbation On input: an uninitialized array
     *                              On output: an array containing a values that has a 
     *                              sinusoidal variation with given amplitude and mode. 
     * @param[in] mode The mode of the perturbation. 
     * @param[in] perturb_amplitude The amplitude of the perturbation. 
     */
    void perturbation_initialization(
            device_t<DSpanX> perturbation,
            int mode,
            double perturb_amplitude) const;

    /**
     * @brief Creates an instance of the SingleModePerturbInitialization class.
     * @param[in] fequilibrium A Maxwellian. 
     * @param[in] init_perturb_mode The perturbation mode. 
     * @param[in] init_perturb_amplitude The perturbation amplitude. 
     */
    SingleModePerturbInitialization(
            device_t<DViewSpVx> fequilibrium,
            ViewSp<int> init_perturb_mode,
            DViewSp init_perturb_amplitude);

    ~SingleModePerturbInitialization() override = default;

    /**
     * @brief Initializes the distribution function as as a perturbed Maxwellian. 
     * @param[in, out] allfdistribu The initialized distribution function.
     * @return The initialized distribution function.
     */
    device_t<DSpanSpXVx> operator()(device_t<DSpanSpXVx> allfdistribu) const override;
};
