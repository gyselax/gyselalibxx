// SPDX-License-Identifier: MIT

#pragma once

#include <geometry.hpp>
#include <quadrature.hpp>

/**
 * @brief A class that computes fluid moments of the distribution function.
 * 
 * These fluid moments are the density, mean velocity and temperature of 
 * the distribution function.
 */
class FluidMoments
{
private:
    Quadrature<IDimVx> m_integrate_v;

public:
    /**
     * A tag type to indicate that the density should be calculated.
     */
    struct MomentDensity
    {
    };

    /**
     * A tag type to indicate that the velocity should be calculated.
     */
    struct MomentVelocity
    {
    };

    /**
     * A tag type to indicate that the temperature should be calculated.
     */
    struct MomentTemperature
    {
    };

    /**
     * A static instance of MomentDensity that can be used to indicated to the operator()
     * that the density should be calculated.
     */
    static constexpr MomentDensity s_density = MomentDensity();
    /**
     * A static instance of MomentVelocity that can be used to indicated to the operator()
     * that the velocity should be calculated.
     */
    static constexpr MomentVelocity s_velocity = MomentVelocity();
    /**
     * A static instance of MomentTemperature that can be used to indicated to the operator()
     * that the temperature should be calculated.
     */
    static constexpr MomentTemperature s_temperature = MomentTemperature();

    /**
     * The constructor for the operator.
     *
     * @param[in] integrate_v A quadrature method which integrates over the velocity space.
     */
    FluidMoments(Quadrature<IDimVx> integrate_v);

    ~FluidMoments() = default;

    /**
     * Calculate the density at a specific point of the distribution function.
     *
     * @param[out] density The density at the point.
     * @param[in] fdistribu A slice in velocity space of the distribution function
     *                          at the given point.
     * @param[in] moment_density A tag to ensure that the correct operator is called.
     */
    void operator()(double& density, DViewVx fdistribu, MomentDensity moment_density);

    /**
     * Calculate the density of the distribution function.
     *
     * @param[out] density The density at various points for different species.
     * @param[in] allfdistribu The distribution function.
     * @param[in] moment_density A tag to ensure that the correct operator is called.
     */
    void operator()(DSpanSpX density, DViewSpXVx allfdistribu, MomentDensity moment_density);

    /**
     * Calculate the mean velocity at a specific point of the distribution function.
     *
     * @param[out] mean_velocity The mean velocity at the point.
     * @param[in] fdistribu A slice in velocity space of the distribution function
     *                          at the given point.
     * @param[in] density The density at the point.
     * @param[in] moment_velocity A tag to ensure that the correct operator is called.
     */
    void operator()(
            double& mean_velocity,
            DViewVx fdistribu,
            double density,
            MomentVelocity moment_velocity);

    /**
     * Calculate the mean velocity of the distribution function.
     *
     * @param[out] mean_velocity The mean velocity at various points for different species.
     * @param[in] allfdistribu The distribution function.
     * @param[in] density The density at various points for different species.
     * @param[in] moment_velocity A tag to ensure that the correct operator is called.
     */
    void operator()(
            DSpanSpX mean_velocity,
            DViewSpXVx allfdistribu,
            DViewSpX density,
            MomentVelocity moment_velocity);

    /**
     * Calculate the temperature at a specific point of the distribution function.
     *
     * @param[out] temperature The mean temperature at the point.
     * @param[in] fdistribu A slice in velocity space of the distribution function
     *                          at the given point.
     * @param[in] density The density at the point.
     * @param[in] mean_velocity The mean velocity at the point.
     * @param[in] moment_temperature A tag to ensure that the correct operator is called.
     */
    void operator()(
            double& temperature,
            DViewVx fdistribu,
            double density,
            double mean_velocity,
            MomentTemperature moment_temperature);

    /**
     * Calculate the mean temperature of the distribution function.
     *
     * @param[out] temperature The mean temperature at various points for different species.
     * @param[in] allfdistribu The distribution function.
     * @param[in] density The density at various points for different species.
     * @param[in] mean_velocity The mean velocity at various points for different species.
     * @param[in] moment_temperature A tag to ensure that the correct operator is called.
     */
    void operator()(
            DSpanSpX temperature,
            DViewSpXVx allfdistribu,
            DViewSpX density,
            DViewSpX mean_velocity,
            MomentTemperature moment_temperature);
};