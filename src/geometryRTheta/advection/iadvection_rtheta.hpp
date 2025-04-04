// SPDX-License-Identifier: MIT

#pragma once
#include <ddc/ddc.hpp>

#include "ddc_aliases.hpp"
#include "geometry.hpp"

/**
 * @brief Define the base class of 2D advection operators in polar index range.
 */
class IAdvectionRTheta
{
public:
    virtual ~IAdvectionRTheta() = default;

    /**
     * @brief Advect a function along the advection field given on dt
     * with a given advection field along XY.
     *
     * @param[in, out] allfdistribu
     *      The function to be advected.
     * @param[in] advection_field
     *      The advection field along the physical index range axes, XY.
     * @param[in] dt
     *      The time step.
     *
     * @return A Field of the advected function (allfdistribu).
     */
    virtual DFieldRTheta operator()(
            DFieldRTheta allfdistribu,
            DConstVectorFieldRTheta<X, Y> advection_field,
            double const dt) const = 0;

    /**
     * @brief Advect a function along the advection field given on dt
     * with a given advection field along RTheta.
     *
     * @param[in, out] allfdistribu
     *      The function to be advected.
     * @param[in] advection_field
     *      The advection field on the contravariant basis of the logical domain.
     * @param [in] advection_field_xy_centre
     *      A vector in the Cartesian basis, containing the value of the advection
     *      field at the O-point.
     * @param[in] dt
     *      The time step.
     *
     * @return A Field of the advected function (allfdistribu).
     */
    virtual DFieldRTheta operator()(
            DFieldRTheta allfdistribu,
            DConstVectorFieldRTheta<R, Theta> advection_field,
            DVector<X, Y> const& advection_field_xy_centre,
            double const dt) const = 0;
};
