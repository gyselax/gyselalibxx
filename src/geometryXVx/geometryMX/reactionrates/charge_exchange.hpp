// SPDX-License-Identifier: MIT

#pragma once

#include <geometry.hpp>

#include "ireactionrate.hpp"

/**
 * @brief A class that describes a charge-exchange reaction rate.
 */
class ChargeExchangeRate : public IReactionRate
{
private:
    /*
     * @brief Stores the linear coefficients to compute the polynomial for the charge-exchange reaction rate vs. ion temperature.
    */
    static constexpr std::size_t s_coefficients_size = 5;
    DKokkosView<s_coefficients_size> m_cx_coefficients;
    double const m_norm_coeff_rate;

public:
    /**
     * @brief Creates an instance of the ConstantChargeExchangeRate class. 
     * A polynomial of reaction rate vs. ion temperature is fitted in log-log space.
     * The polynomial fits data from the OPEN-ADAS database for Hydrogen.
     * 
     * @param[in] norm_coeff_rate All rate coefficients coefficients are normalized so that CX is
     * of order unity at normalized density (n) = 1 and normalized temperature (T) = 1
     * with n_0 = 1e20 and T_0 = 10 eV.
     * These are typical values for a SOL plasma.
     * The norm_coeff_rate parameter shifts all the reaction rates to modify the source dynamics.
     */
    explicit ChargeExchangeRate(double norm_coeff_rate);

    ~ChargeExchangeRate() override = default;

    /**
     * @brief Compute the charge-exchange reaction rate depending on density and temperature.
     * 
     * @param[out] rate The charge-exchange reaction rates.
     * @param[in] density The plasma density at which the reaction rate is computed. Although inputted, charge-exchange reaction rate is not affected by density.
     * @param[in] temperature The plasma temperature at which the reaction rate is computed.
     * @return The charge-exchange reaction rate.
     */
    DFieldSpX operator()(DFieldSpX rate, DConstFieldSpX density, DConstFieldSpX temperature)
            const override;
};
