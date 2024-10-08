#pragma once
#include <functional>

#include "ddc_aliases.hpp"
#include "geometry.hpp"

/**
 * @brief Define a base class for all the time integration methods used for the advection.
 *
 * @see BslAdvectionRTheta
 */
class IFootFinder
{
public:
    virtual ~IFootFinder() = default;

    /**
     * @brief Advect the feet over @f$ dt @f$.
     *
     * @param[in, out] feet
     *      On input: the mesh points.
     *      On output: the characteristic feet.
     * @param[in] advection_field
     *      The advection field in the physical index range.
     * @param[in] dt
     *      The time step.
     */
    virtual void operator()(
            host_t<FieldRTheta<CoordRTheta>> feet,
            host_t<DConstVectorFieldRTheta<X, Y>> advection_field,
            double dt) const = 0;
};
