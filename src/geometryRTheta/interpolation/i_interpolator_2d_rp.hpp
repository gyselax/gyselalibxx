// SPDX-License-Identifier: MIT
#pragma once

#include <memory>

#include "geometry.hpp"



/**
 * @brief A class which provides an interpolating function.
 *
 * An abstract class which implements a function allowing
 * the value of a function to be approximated at a set of
 * coordinates from a set of known values of the function.
 */
class IInterpolatorRTheta
{
public:
    virtual ~IInterpolatorRTheta() = default;

    /**
     * @brief Approximate the value of a function at a set of coordinates using the
     * current values at a known set of interpolation points.
     *
     * @param[in, out] inout_data
     * 					On input: an array containing the value of the function at the interpolation points.
     * 					On output: an array containing the value of the function at the coordinates.
     * @param[in] coordinates
     * 					The coordinates where the function should be evaluated.
     *
     * @return A reference to the inout_data array containing the value of the function at the coordinates.
     */
    virtual host_t<DFieldRTheta> operator()(
            host_t<DFieldRTheta> inout_data,
            host_t<ConstFieldRTheta<CoordRTheta>> coordinates) const = 0;
};



/**
 * @brief A class which provides access to an interpolating function which can be preallocated where useful.
 *
 * An abstract class which implements a preallocate function returning a pointer to an InterpolatorRTheta.
 * A pointer to an InterpolatorRTheta is used so that the returned object can be any sub-class of IInterpolatorRTheta.
 * The type (and thus the implementation of the operator) will be determined when the pointer is
 * dereferenced.
 *
 * The preallocate function should be used to allocate an instance of the IInterpolatorRTheta before
 * using it repeatedly. Once the preallocated object goes out of scope it will be deallocated.
 * This means that objects of this class take up little or no space in memory.
 *
 */
class IPreallocatableInterpolatorRTheta : public IInterpolatorRTheta
{
public:
    ~IPreallocatableInterpolatorRTheta() override = default;

    /**
     * @brief Allocate an instance of a pointer to an InterpolatorRTheta.
     *
     * Allocate and return an instance of a pointer to a sub-class
     * of IInterpolatorRTheta.
     *
     * @return A pointer to an InterpolatorRTheta.
     *
     * @see IInterpolatorRTheta
     */
    virtual std::unique_ptr<IInterpolatorRTheta> preallocate() const = 0;

    host_t<DFieldRTheta> operator()(
            host_t<DFieldRTheta> const inout_data,
            host_t<ConstFieldRTheta<CoordRTheta>> const coordinates) const override
    {
        return (*preallocate())(inout_data, coordinates);
    }
};
