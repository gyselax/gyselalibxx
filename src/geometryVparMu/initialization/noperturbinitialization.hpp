// SPDX-License-Identifier: MIT

#pragma once

#include <paraconf.h>
#include <species_info.hpp>

#include "geometry.hpp"
#include "iinitialization.hpp"
#include "paraconfpp.hpp"

/// Initialization operator with no perturbation, i.e the distribution function equal to the Maxwellian
class NoPerturbInitialization : public IInitialization
{
    DConstFieldSpVparMu m_fequilibrium;

public:
    /**
     * @brief Creates an instance of the NoPerturbInitialization class.
     * @param[in] fequilibrium A Maxwellian. 
     */
    NoPerturbInitialization(DConstFieldSpVparMu fequilibrium);

    ~NoPerturbInitialization() override = default;

    /**
     * @brief Initializes the distribution function as as a  Maxwellian. 
     * @param[in, out] allfdistribu The initialized distribution function.
     * @return The initialized distribution function.
     */
    DFieldSpVparMu operator()(DFieldSpVparMu allfdistribu) const override;
};
