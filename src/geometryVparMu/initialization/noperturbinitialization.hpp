// SPDX-License-Identifier: MIT

#pragma once

#include <paraconf.h>

#include "geometry.hpp"
#include "iinitialization.hpp"
#include "paraconfpp.hpp"
#include "species_info.hpp"

/// Initialisation operator with no perturbation, i.e the distribution function equal to the Maxwellian
class NoPerturbInitialisation : public IInitialisation
{
    DConstFieldSpVparMu m_fequilibrium;

public:
    /**
     * @brief Creates an instance of the NoPerturbInitialisation class.
     * @param[in] fequilibrium A Maxwellian. 
     */
    NoPerturbInitialisation(DConstFieldSpVparMu fequilibrium);

    ~NoPerturbInitialisation() override = default;

    /**
     * @brief Initialises the distribution function as as a  Maxwellian. 
     * @param[in, out] allfdistribu The initialised distribution function.
     * @return The initialised distribution function.
     */
    DFieldSpVparMu operator()(DFieldSpVparMu allfdistribu) const override;
};
