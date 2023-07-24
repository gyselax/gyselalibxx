// SPDX-License-Identifier: MIT

#pragma once

#include <memory>

#include <ddc/ddc.hpp>

/**
 * @brief A class which provides an interpolating function.
 *
 * An abstract class which implements a function allowing
 * the value of a function to be approximated at a set of
 * coordinates from a set of known values of the function.
 */
template <class DDim>
class IInterpolator
{
    using CDim = typename DDim::continuous_dimension_type;

public:
    virtual ~IInterpolator() = default;

    /**
     * @brief Approximate the value of a function at a set of coordinates using the
     * current values at a known set of interpolation points.
     *
     * @param[in, out] inout_data On input: an array containing the value of the function at the interpolation points.
     * 			 On output: an array containing the value of the function at the coordinates.
     * @param[in] coordinates The coordinates where the function should be evaluated.
     *
     * @return A reference to the inout_data array containing the value of the function at the coordinates.
     */
    virtual ddc::ChunkSpan<double, ddc::DiscreteDomain<DDim>> operator()(
            ddc::ChunkSpan<double, ddc::DiscreteDomain<DDim>> inout_data,
            ddc::ChunkSpan<const ddc::Coordinate<CDim>, ddc::DiscreteDomain<DDim>> coordinates)
            const = 0;
};

/**
 * @brief A class which provides access to an interpolating function which can be preallocated where useful.
 *
 * An abstract class which implements a preallocate function returning an unique pointer to an IInterpolator.
 * A pointer is used so that the returned object can be any sub-class of IInterpolator.
 * The type (and thus the implementation of the operator) will be determined when the pointer is
 * dereferenced.
 *
 * The preallocate function should be used to allocate an instance of the IInterpolator before
 * using it repeatedly. Once the preallocated object goes out of scope it will be deallocated.
 * This means that objects of this class take up little or no space in memory.
 *
 * An example of this is seen in BslAdvectionVelocity. The IPreallocatableInterpolator stored in the
 * BslAdvectionVelocity takes up no memory between advections, however during the execution of the
 * BslAdvectionVelocity::operator() function the IPreallocatableInterpolator::preallocate() function
 * is called. This leads to the creation of an IInterpolator instance, ensuring that all buffers necessary
 * for the interpolation during the advection are allocated before the IInterpolator is used for
 * interpolation in the advection loop. This ensures that these buffers are only allocated once per
 * advection at the start of the BslAdvectionVelocity::operator() function. At the end of this function
 * the unique pointer goes out of scope and the buffers are deallocated.
 */
template <class DDim>
class IPreallocatableInterpolator : public IInterpolator<DDim>
{
    using CDim = typename DDim::continuous_dimension_type;

public:
    ~IPreallocatableInterpolator() override = default;

    /**
     * @brief Allocate an instance of an InterpolatorProxy to use as an IInterpolator.
     *
     * Allocate and return an unique pointer to an instance of an IInterpolator.
     *
     * @return An allocated instance of an InterpolatorProxy.
     *
     * @see InterpolatorProxy
     * @see IInterpolator
     */
    virtual std::unique_ptr<IInterpolator<DDim>> preallocate() const = 0;

    ddc::ChunkSpan<double, ddc::DiscreteDomain<DDim>> operator()(
            ddc::ChunkSpan<double, ddc::DiscreteDomain<DDim>> const inout_data,
            ddc::ChunkSpan<const ddc::Coordinate<CDim>, ddc::DiscreteDomain<DDim>> const
                    coordinates) const override
    {
        return (*preallocate())(inout_data, coordinates);
    }
};
