// SPDX-License-Identifier: MIT

#pragma once

#include <geometry.hpp>

#include "ireactionrate.hpp"

/**
 * @brief A class that describes a recombination reaction rate.
 */
class RecombinationRate : public IReactionRate
{
private:
    /**
     * @brief Stores the coefficients depending on electron density for calculating the recombination reaction rate polynomial as a function of electron temperature.
     */
    static constexpr std::size_t s_coefficients_size = 2;
    DKokkosView<s_coefficients_size> m_slope_coefficient;
    DKokkosView<s_coefficients_size> m_intercept_coefficient;

    double const m_norm_coeff_rate;

public:
    /**
     * @brief Creates an instance of the ConstantRecombinationRate class. 
     * A polynomial of reaction rate vs. electron temperature is fitted.
     * Coefficients are fitted in semi-log space, while the rate vs. electron temperature polynomial is fitted in log-log space.
     * The polynomial fits data from the OPEN-ADAS database for Hydrogen.
     * 
     * @param[in] norm_coeff_rate All rate coefficients are normalized so that CX is of order unity at normalized density (n) = 1 and normalized temperature (T) = 1 with n_0 = 1e20 and T_0 = 10 eV.
     * These are typical values for a SOL plasma.
     * The norm_coeff_rate parameter shifts all the reaction rates to modify the source dynamics.
     */
    explicit RecombinationRate(double norm_coeff_rate);

    ~RecombinationRate() override = default;

    /**
     * @brief Compute the recombination reaction rate depending on density and temperature.
     * 
     * @param[out] rate The recombination reaction rates.
     * @param[in] density The plasma density at which the reaction rate is computed.
     * @param[in] temperature The plasma temperature at which the reaction rate is computed.
     * @return The recombination reaction rate.
     */
    DSpanSpX operator()(DSpanSpX rate, DViewSpX density, DViewSpX temperature) const override;
};
