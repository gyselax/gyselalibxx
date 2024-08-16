// SPDX-License-Identifier: MIT

#pragma once

#include <memory>

#include <ddc/ddc.hpp>
#include <ddc/kernels/splines/deriv.hpp>

#include "ddc_aliases.hpp"
#include "ddc_helper.hpp"

// TODO: Generalize (IDimI -> Tags...) and make it usable for all Gysela operators ?
template <template <class...> class Interp, class GridInterp, class IdxRange>
struct interpolator_on_idx_range
{
};

/**
 * A structure which builds an interpolation type.
 *
 * @tparam Interp The interpolator class being built.
 * @tparam GridInterp The dimension along which the operator will interpolate.
 * @tparam Grid1D... The dimensions on which the data being interpolated are defined.
 */
template <template <class...> class Interp, class GridInterp, class... Grid1D>
struct interpolator_on_idx_range<Interp, GridInterp, IdxRange<Grid1D...>>
{
    /// The type of the interpolator
    using type = Interp<GridInterp, Grid1D...>;
};

/**
 * A template function which returns an interpolation type.
 *
 * @tparam Interp The interpolator class being built.
 * @tparam GridInterp The dimension along which the operator will interpolate.
 * @tparam IdxRange The index range on which the data being interpolated is defined.
 */
template <template <class...> class Interp, class GridInterp, class IdxRange>
using interpolator_on_idx_range_t =
        typename interpolator_on_idx_range<Interp, GridInterp, IdxRange>::type;

/**
 * @brief A class which provides an interpolating function.
 *
 * An abstract class which implements a function allowing
 * the value of a function to be approximated at a set of
 * coordinates from a set of known values of the function.
 */
template <class GridInterp, class... Grid1D>
class IInterpolator
{
public:
    virtual ~IInterpolator() = default;

    /// @brief The type of the dimension representing derivatives.
    using deriv_type = ddc::Deriv<typename GridInterp::continuous_dimension_type>;
    /// @brief The type of the whole index range on which derivatives are defined.
    using batched_derivs_idx_range_type
            = ddc::replace_dim_of_t<IdxRange<Grid1D...>, GridInterp, deriv_type>;

    /**
     * @brief Get the batched derivs index range on lower boundaries.
     *
     * Dimension of interest IDimI is replaced with ddc::Deriv<IDimI::continuous_dimensions_type>.
     * This is the index range on which derivatives on lower boundaries are defined.
     *
     * @param[in] idx_range The index range of a single-species distribution function.
     * @return idx_range The lower boundaries of this index range.
     */
    virtual batched_derivs_idx_range_type batched_derivs_idx_range_xmin(
            IdxRange<Grid1D...> idx_range) const = 0;

    /**
     * @brief Get the batched derivs index range on upper boundaries.
     *
     * Dimension of interest IDimI is replaced with ddc::Deriv<IDimI::continuous_dimensions_type>.
     * This is the index range on which derivatives on upper boundaries are defined.
     *
     * @param[in] idx_range The index range of a single-species distribution function.
     * @return idx_range The upper boundaries of this index range.
     */
    virtual batched_derivs_idx_range_type batched_derivs_idx_range_xmax(
            IdxRange<Grid1D...> idx_range) const = 0;

    /**
     * @brief Approximate the value of a function at a set of coordinates using the
     * current values at a known set of interpolation points.
     *
     * @param[in, out] inout_data On input: an array containing the value of the function at the interpolation points.
     * 			 On output: an array containing the value of the function at the coordinates.
     * @param[in] coordinates The coordinates where the function should be evaluated.
     * @param[in] derivs_xmin The values of the derivatives at the lower boundary
     * (used only with splines and ddc::BoundCond::HERMITE lower boundary condition).
     * @param[in] derivs_xmax The values of the derivatives at the upper boundary
     * (used only with splines and ddc::BoundCond::HERMITE upper boundary condition).
     *
     * @return A reference to the inout_data array containing the value of the function at the coordinates.
     */
    virtual Field<double, IdxRange<Grid1D...>> operator()(
            Field<double, IdxRange<Grid1D...>> inout_data,
            Field<const Coord<typename GridInterp::continuous_dimension_type>, IdxRange<Grid1D...>>
                    coordinates,
            std::optional<Field<double const, batched_derivs_idx_range_type>> derivs_xmin
            = std::nullopt,
            std::optional<Field<double const, batched_derivs_idx_range_type>> derivs_xmax
            = std::nullopt) const = 0;
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
template <class GridInterp, class... Grid1D>
class IPreallocatableInterpolator : public IInterpolator<GridInterp, Grid1D...>
{
public:
    ~IPreallocatableInterpolator() override = default;

    /// @brief The type of the dimension representing derivatives.
    using deriv_type = typename IInterpolator<GridInterp, Grid1D...>::deriv_type;
    /// @brief The type of the whole index range on which derivatives are defined.
    using batched_derivs_idx_range_type =
            typename IInterpolator<GridInterp, Grid1D...>::batched_derivs_idx_range_type;

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
    virtual std::unique_ptr<IInterpolator<GridInterp, Grid1D...>> preallocate() const = 0;

    /**
     * @brief Get the batched derivs index range on lower boundaries.
     *
     * Dimension of interest IDimI is replaced with ddc::Deriv<IDimI::continuous_dimensions_type>.
     * This is the index range on which derivatives on lower boundaries are defined.
     *
     * @param[in] idx_range The index range of a single-species distribution function.
     * @return idx_range The lower boundaries of this index range.
     */
    batched_derivs_idx_range_type batched_derivs_idx_range_xmin(
            IdxRange<Grid1D...> idx_range) const override
    {
        return (*preallocate()).batched_derivs_idx_range_xmin(idx_range);
    }

    /**
     * @brief Get the batched derivs index range on upper boundaries.
     *
     * Dimension of interest IDimI is replaced with ddc::Deriv<IDimI::continuous_dimensions_type>.
     * This is the index range on which derivatives on upper boundaries are defined.
     *
     * @param[in] idx_range The index range of a single-species distribution function.
     * @return idx_range The upper boundaries of this index range.
     */
    batched_derivs_idx_range_type batched_derivs_idx_range_xmax(
            IdxRange<Grid1D...> idx_range) const override
    {
        return (*preallocate()).batched_derivs_idx_range_xmax(idx_range);
    }

    /**
     * @brief Approximate the value of a function at a set of coordinates using the
     * current values at a known set of interpolation points by temporarily preallocating an IInterpolator.
     *
     * @param[in, out] inout_data On input: an array containing the value of the function at the interpolation points.
     * 			 On output: an array containing the value of the function at the coordinates.
     * @param[in] coordinates The coordinates where the function should be evaluated.
     * @param[in] derivs_xmin The values of the derivatives at the lower boundary
     * (used only with splines and ddc::BoundCond::HERMITE lower boundary condition).
     * @param[in] derivs_xmax The values of the derivatives at the upper boundary
     * (used only with splines and ddc::BoundCond::HERMITE upper boundary condition).
     *
     * @return A reference to the inout_data array containing the value of the function at the coordinates.
     */
    Field<double, IdxRange<Grid1D...>> operator()(
            Field<double, IdxRange<Grid1D...>> const inout_data,
            Field<const Coord<typename GridInterp::continuous_dimension_type>,
                  IdxRange<Grid1D...>> const coordinates,
            std::optional<Field<double const, batched_derivs_idx_range_type>> derivs_xmin
            = std::nullopt,
            std::optional<Field<double const, batched_derivs_idx_range_type>> derivs_xmax
            = std::nullopt) const override
    {
        return (*preallocate())(inout_data, coordinates);
    }
};
