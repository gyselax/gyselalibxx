// SPDX-License-Identifier: MIT

#pragma once

#include <geometry.hpp>
#include <species_info.hpp>

#include "iinitialization.hpp"

/// Initialization operator with a sinusoidal perturbation of a Maxwellian. This initializes all species.
class SingleModePerturbInitialization : public IInitialization
{
    device_t<DViewSpVx> m_fequilibrium;

    ViewSp<int> m_init_perturb_mode;

    DViewSp m_init_perturb_amplitude;

public:
    /*
      Initialization of the perturbation
    */
    void perturbation_initialization(
            device_t<DSpanX> perturbation,
            int mode,
            double perturb_amplitude) const;

    // Warning: all variables shall remain valid until the last call to `operator()`
    SingleModePerturbInitialization(
            device_t<DViewSpVx> fequilibrium,
            ViewSp<int> init_perturb_mode,
            DViewSp init_perturb_amplitude);

    ~SingleModePerturbInitialization() override = default;

    device_t<DSpanSpXVx> operator()(device_t<DSpanSpXVx> allfdistribu) const override;
};
