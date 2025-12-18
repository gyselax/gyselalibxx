// SPDX-License-Identifier: MIT

#pragma once

#include "geometry.hpp"
#include "quadrature.hpp"

/**
 * @brief A class that computes fluid moments of the distribution function.
 * 
 * These fluid moments are the density, mean velocity and temperature of 
 * the distribution function integrated over the 2D velocity space (vpar, mu).
 */
class FluidMoments
{
private:
    Quadrature<IdxRangeVparMu, IdxRangeSpTor3DV2D, Kokkos::DefaultExecutionSpace::memory_space>
            m_integrate_vparmu;

public:
    struct MomentDensity
    {
    };
    struct MomentVelocity
    {
    };
    struct MomentTemperature
    {
    };
    static constexpr MomentDensity s_density = MomentDensity();
    static constexpr MomentVelocity s_velocity = MomentVelocity();
    static constexpr MomentTemperature s_temperature = MomentTemperature();

    /**
     * The constructor for the operator.
     *
     * @param[in] integrate_vparmu A quadrature method which integrates over the velocity space (vpar, mu).
     */
    explicit FluidMoments(Quadrature<
                 IdxRangeVparMu,
                 IdxRangeSpTor3DV2D,
                 Kokkos::DefaultExecutionSpace::memory_space> integrate_vparmu);

    ~FluidMoments() = default;

    /**
     * Calculate the density of the distribution function.
     *
     * @param[out] density The density for each species at each spatial point (tor1, tor2, tor3).
     * @param[in] allfdistribu The distribution function (species × tor1 × tor2 × tor3 × vpar × mu).
     * @param[in] moment_density A tag to ensure that the correct operator is called.
     */
    void operator()(
            DField<IdxRangeSpTor3D> const density,
            DConstField<IdxRangeSpTor3DV2D> const allfdistribu,
            MomentDensity moment_density);

    /**
     * Calculate the mean velocity (vpar) of the distribution function.
     *
     * @param[out] mean_velocity The mean velocity for each species at each spatial point.
     * @param[in] allfdistribu The distribution function.
     * @param[in] density The density for each species at each spatial point.
     * @param[in] moment_velocity A tag to ensure that the correct operator is called.
     */
    void operator()(
            DField<IdxRangeSpTor3D> const mean_velocity,
            DConstField<IdxRangeSpTor3DV2D> const allfdistribu,
            DConstField<IdxRangeSpTor3D> const density,
            MomentVelocity moment_velocity);

    /**
     * Calculate the temperature of the distribution function.
     *
     * @param[out] temperature The mean temperature for each species at each spatial point.
     * @param[in] allfdistribu The distribution function.
     * @param[in] density The density for each species at each spatial point.
     * @param[in] mean_velocity The mean velocity for each species at each spatial point.
     * @param[in] moment_temperature A tag to ensure that the correct operator is called.
     */
    void operator()(
            DField<IdxRangeSpTor3D> const temperature,
            DConstField<IdxRangeSpTor3DV2D> const allfdistribu,
            DConstField<IdxRangeSpTor3D> const density,
            DConstField<IdxRangeSpTor3D> const mean_velocity,
            MomentTemperature moment_temperature);
};
