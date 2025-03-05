// SPDX-License-Identifier: MIT
#pragma once

#include <memory>

#include "ddc_aliases.hpp"



/**
 * @brief A class which provides an interpolating function.
 *
 * An abstract class which implements a function allowing
 * the value of a function to be approximated at a set of
 * coordinates from a set of known values of the function.
 */
template <class IdxRange2D, class IdxRangeBatched>
class IInterpolator2D
{
    static_assert(IdxRange2D::rank() == 2);
    using Dim1 = typename ddc::type_seq_element_t<0, ddc::to_type_seq_t<IdxRange2D>>::
            continuous_dimension_type;
    using Dim2 = typename ddc::type_seq_element_t<1, ddc::to_type_seq_t<IdxRange2D>>::
            continuous_dimension_type;

public:
    using CoordType = Coord<Dim1, Dim2>;
    using DFieldType = DField<IdxRangeBatched>;
    using CConstFieldType = ConstField<CoordType, IdxRangeBatched>;

public:
    virtual ~IInterpolator2D() = default;

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
    virtual DField<IdxRangeBatched> operator()(
            DField<IdxRangeBatched> inout_data,
            ConstField<CoordType, IdxRangeBatched> coordinates) const = 0;
};



/**
 * @brief A class which provides access to an interpolating function which can be preallocated where useful.
 *
 * An abstract class which implements a preallocate function returning a pointer to an InterpolatorRTheta.
 * A pointer to an InterpolatorRTheta is used so that the returned object can be any sub-class of IInterpolator2D.
 * The type (and thus the implementation of the operator) will be determined when the pointer is
 * dereferenced.
 *
 * The preallocate function should be used to allocate an instance of the IInterpolator2D before
 * using it repeatedly. Once the preallocated object goes out of scope it will be deallocated.
 * This means that objects of this class take up little or no space in memory.
 *
 */
template <class IdxRange2D, class IdxRangeBatched>
class IPreallocatableInterpolator2D : public IInterpolator2D<IdxRange2D, IdxRangeBatched>
{
public:
    using typename IInterpolator2D<IdxRange2D, IdxRangeBatched>::CoordType;

public:
    ~IPreallocatableInterpolator2D() override = default;

    /**
     * @brief Allocate an instance of a pointer to an InterpolatorRTheta.
     *
     * Allocate and return an instance of a pointer to a sub-class
     * of IInterpolator2D.
     *
     * @return A pointer to an InterpolatorRTheta.
     *
     * @see IInterpolator2D
     */
    virtual std::unique_ptr<IInterpolator2D<IdxRange2D, IdxRangeBatched>> preallocate() const = 0;

    DField<IdxRangeBatched> operator()(
            DField<IdxRangeBatched> inout_data,
            ConstField<CoordType, IdxRangeBatched> coordinates) const override
    {
        return (*preallocate())(inout_data, coordinates);
    }
};
