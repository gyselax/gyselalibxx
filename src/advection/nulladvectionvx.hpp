// SPDX-License-Identifier: MIT

#pragma once
#include "ddc_aliases.hpp"
#include "iadvectionvx.hpp"

/**
 * @brief This is a class which imitates a velocity advection.
 * It inherits from IAdvectionV and can be used as an advection operator but does not
 * actually modify the distribution function. This can be useful for debugging purposes.
 */
template <class FdistribuIdxRange, class SpatialIdxRange>
class NullAdvectionVelocity : public IAdvectionV<FdistribuIdxRange, SpatialIdxRange>
{
public:
    NullAdvectionVelocity() = default;

    ~NullAdvectionVelocity() override = default;

    /**
     * @brief Do nothing instead of advecting fdistribu along GridV for a duration dt.
     *
     * @param[in, out] allfdistribu Reference to an array containing the value of the distribution function.
     * @param[in] electric_field The electric field which derives from electrostatic potential and is the advection speed.
     * @param[in] dt Time step
     *
     * @return A reference to the allfdistribu array containing the value of the function at the coordinates.
     */
    Field<double, FdistribuIdxRange> operator()(
            Field<double, FdistribuIdxRange> allfdistribu,
            [[maybe_unused]] Field<const double, SpatialIdxRange> electric_field,
            [[maybe_unused]] double dt) const override
    {
        return allfdistribu;
    }
};
