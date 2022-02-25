// SPDX-License-Identifier: MIT

#pragma once

#include <geometry.hpp>
#include <species_info.hpp>

#include "iinitialization.hpp"

/// Initialization operator with a sinusoidal perturbation of a Maxwellian. This initializes all species.
class SingleModePerturbInitialization : public IInitialization
{
    SpeciesInformation const& m_species_info;

    ViewSp<int> m_init_perturb_mode;

    DViewSp m_init_perturb_amplitude;

private:
    /*
      Initialization of the perturbation
    */
    void perturbation_initialization(DSpanX perturbation, int mode, double perturb_amplitude) const;

public:
    // Warning: all variables shall remain valid until the last call to `operator()`
    SingleModePerturbInitialization(
            SpeciesInformation const& species_info,
            ViewSp<int> init_perturb_mode,
            DViewSp init_perturb_amplitude);

    ~SingleModePerturbInitialization() override = default;

    DSpanSpXVx operator()(DSpanSpXVx allfdistribu) const override;
};
