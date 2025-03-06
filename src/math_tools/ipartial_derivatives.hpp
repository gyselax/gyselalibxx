// SPDX-License-Identifier: MIT
#pragma once

#include <ddc/ddc.hpp>

#include "ddc_alias_inline_functions.hpp"
#include "ddc_aliases.hpp"
#include "type_seq_tools.hpp"


/**
 * @brief A super class for a partial derivative operator.
 * @tparam DFieldType The type of the field on which the operator acts.
 * @tparam DerivDirection The dimension Xi on which the partial derivative is calculated.
 */
template <class IdxRange, class DerivativeDirection>
class IPartialDerivative
{
    static_assert(ddc::is_borrowed_chunk_v<DFieldType>);

public:
    using DerivativeDirection = typename DerivativeDirection;
    /// The index range on which this operator acts.
    using IdxRange = typename IdxRange;

    /// The type of the object that will be differentiated.
    using DFieldType = DField<IdxRange>;

    /// The type of the calculated derivative.
    using DConstFieldType = typename DFieldType::view_type;

    /// The type of the grid on the dimension Xi on which the partial derivative is calculated.
    using GridDerivativeDirection
            = find_grid_t<DerivativeDirection, ddc::to_type_seq_t<IdxRange>>;

    /// The index range of the dimension Xi on which the partial derivative is calculated.
    using IdxRangeDeriv = IdxRange<GridDerivativeDirection>;

    /// The index range of all dimensions except Xi.
    using IdxRangeBatch = ddc::remove_dims_of_t<IdxRange, GridDerivativeDirection>;

public:
    /**
    * @brief Compute the partial derivative of @f$ F(X1,..,Xn)@f$ in Xi direction.
    *
    * @param[out] dfieldval_dxi Partial derivatives in Xi direction.
    * @param[in] fieldval Values of the field @f$ F(X1,..,Xn)@f$.
    */
    virtual DFieldType operator()(DField<IdxRangeBatched> fieldval,
            ConstField<CoordType, IdxRange> coordinates) const = 0;
};

/**
 * @brief A class which provides access to a partial derivative function which can be preallocated where useful.
 *
 * An abstract class which implements a preallocate function returning a pointer to an IPartialDerivative.
 * A pointer to an IPartialDerivative is used so that the returned object can be any sub-class of IPartialDerivative.
 * The type (and thus the implementation of the operator) will be determined when the pointer is
 * dereferenced.
 *
 * The preallocate function should be used to allocate an instance of the IPartialDerivative before
 * using it repeatedly. Once the preallocated object goes out of scope it will be deallocated.
 * This means that objects of this class take up little or no space in memory.
 *
 */
template <class IdxRange, class DerivDirection>
class IPreallocatablePartialDerivative : public IPartialDerivative<IdxRange, DerivativeDirection>
{

public:
    /**
     * @brief Allocate an instance of a pointer to an IPartialDerivative.
     *
     * Allocate and return an instance of a pointer to a sub-class
     * of IPartialDerivative.
     *
     * @return A pointer to an IPartialDerivative.
     *
     * @see IPartialDerivative
     */
    virtual std::unique_ptr<IPartialDerivative<IdxRange, DerivativeDirection>> preallocate() const = 0;

    /**
     * @brief Computes the partial derivative of a function at each coordinate 
     * where the value of the function is known.
     *
     * @param[out] dfieldval_dxi
     * 					On output: an array containing the value of the partial derivative.
     * @param[in] fieldval
     * 					An array containing the value of the function to be differentiated.
     *
     * @return A reference to the inout_data array containing the value of the function at the coordinates.
     */
    DField<IdxRange> operator()(
            DField<IdxRange> fieldval,
            ConstField<CoordType, IdxRange> coordinates) const override
    {
        return (*preallocate())(fieldval, coordinates);
    }
};
