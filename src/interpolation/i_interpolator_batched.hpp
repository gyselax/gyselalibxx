// SPDX-License-Identifier: MIT

#pragma once

#include <memory>

#include <ddc/ddc.hpp>

#include <ddc_helper.hpp>

// TODO: Generalize (IDimI -> Tags...) and make it usable for all Gysela operators ?
template <template <class...> class Interp, class DDimI, class Domain>
struct interpolator_on_domain
{
};

/**
 * A structure which builds an interpolation type.
 *
 * @tparam Interp The batched interpolator class being built.
 * @tparam DDimI The dimension along which the operator will interpolate.
 * @tparam DDim... The dimensions on which the data being interpolated are defined.
 */
template <template <class...> class Interp, class DDimI, class... DDim>
struct interpolator_on_domain<Interp, DDimI, ddc::DiscreteDomain<DDim...>>
{
    /// The type of the interpolator
    using type = Interp<DDimI, DDim...>;
};

/**
 * A template function which returns an interpolation type.
 *
 * @tparam Interp The batched interpolator class being built.
 * @tparam DDimI The dimension along which the operator will interpolate.
 * @tparam Domain The domain on which the data being interpolated is defined.
 */
template <template <class...> class Interp, class DDimI, class Domain>
using interpolator_on_domain_t = typename interpolator_on_domain<Interp, DDimI, Domain>::type;

/**
 * @brief A class which provides an interpolating function.
 *
 * An abstract class which implements a function allowing
 * the value of a function to be approximated at a set of
 * coordinates from a set of known values of the function.
 */
template <class DDimI, class... DDim>
class IInterpolatorBatched
{
public:
    virtual ~IInterpolatorBatched() = default;

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
    virtual device_t<ddc::ChunkSpan<double, ddc::DiscreteDomain<DDim...>>> operator()(
            device_t<ddc::ChunkSpan<double, ddc::DiscreteDomain<DDim...>>> inout_data,
            device_t<ddc::ChunkSpan<
                    const ddc::Coordinate<typename DDimI::continuous_dimension_type>,
                    ddc::DiscreteDomain<DDim...>>> coordinates) const = 0;
};

/**
 * @brief A class which provides access to an interpolating function which can be preallocated where useful.
 *
 * An abstract class which implements a preallocate function returning an unique pointer to an IInterpolatorBatched.
 * A pointer is used so that the returned object can be any sub-class of IInterpolatorBatched.
 * The type (and thus the implementation of the operator) will be determined when the pointer is
 * dereferenced.
 *
 * The preallocate function should be used to allocate an instance of the IInterpolatorBatched before
 * using it repeatedly. Once the preallocated object goes out of scope it will be deallocated.
 * This means that objects of this class take up little or no space in memory.
 *
 * An example of this is seen in BslAdvectionVelocity. The IPreallocatableInterpolator stored in the
 * BslAdvectionVelocity takes up no memory between advections, however during the execution of the
 * BslAdvectionVelocity::operator() function the IPreallocatableInterpolator::preallocate() function
 * is called. This leads to the creation of an IInterpolatorBatched instance, ensuring that all buffers necessary
 * for the interpolation during the advection are allocated before the IInterpolatorBatched is used for
 * interpolation in the advection loop. This ensures that these buffers are only allocated once per
 * advection at the start of the BslAdvectionVelocity::operator() function. At the end of this function
 * the unique pointer goes out of scope and the buffers are deallocated.
 */
template <class DDimI, class... DDim>
class IPreallocatableInterpolatorBatched : public IInterpolatorBatched<DDimI, DDim...>
{
public:
    ~IPreallocatableInterpolatorBatched() override = default;

    /**
     * @brief Allocate an instance of an InterpolatorProxy to use as an IInterpolatorBatched.
     *
     * Allocate and return an unique pointer to an instance of an IInterpolatorBatched.
     *
     * @return An allocated instance of an InterpolatorProxy.
     *
     * @see InterpolatorProxy
     * @see IInterpolatorBatched
     */
    virtual std::unique_ptr<IInterpolatorBatched<DDimI, DDim...>> preallocate() const = 0;

    /**
     * @brief Approximate the value of a function at a set of coordinates using the
     * current values at a known set of interpolation points by temporarily preallocating an IInterpolatorBatched.
     *
     * @param[in, out] inout_data On input: an array containing the value of the function at the interpolation points.
     * 			 On output: an array containing the value of the function at the coordinates.
     * @param[in] coordinates The coordinates where the function should be evaluated.
     *
     * @return A reference to the inout_data array containing the value of the function at the coordinates.
     */
    device_t<ddc::ChunkSpan<double, ddc::DiscreteDomain<DDim...>>> operator()(
            device_t<ddc::ChunkSpan<double, ddc::DiscreteDomain<DDim...>>> const inout_data,
            device_t<ddc::ChunkSpan<
                    const ddc::Coordinate<typename DDimI::continuous_dimension_type>,
                    ddc::DiscreteDomain<DDim...>>> const coordinates) const override
    {
        return (*preallocate())(inout_data, coordinates);
    }
};
