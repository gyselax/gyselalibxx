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
template <class IdxRangeType, class DerivativeDim>
class IPartialDerivative
{
    static_assert(ddc::is_borrowed_chunk_v<DField<IdxRangeType>>);

public:
    /// The type of the object that will be differentiated.
    using DFieldType = DField<IdxRangeType>;
    using DFieldMemType = DFieldMem<IdxRangeType>;

    /// The type of the calculated derivative.
    using DConstFieldType = typename DFieldType::view_type;

    /// The type of the grid on the dimension Xi on which the partial derivative is calculated.
    using GridDerivativeDim = find_grid_t<DerivativeDim, ddc::to_type_seq_t<IdxRangeType>>;

    /// The index range of the dimension Xi on which the partial derivative is calculated.
    using IdxRangeDeriv = IdxRange<GridDerivativeDim>;

    /// The index range of all dimensions except Xi.
    using IdxRangeBatched = ddc::remove_dims_of_t<IdxRangeType, GridDerivativeDim>;

public:
    /**
    * @brief Compute the partial derivative of @f$ F(X1,..,Xn)@f$ in Xi direction.
    *
    * @param[in] fieldval Values of the field @f$ F(X1,..,Xn)@f$.
    */
    virtual DFieldType operator()(DConstFieldType fieldval) const = 0;
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
template <class IdxRangeType, class DerivativeDim>
class IPartialDerivativeCreator : public IPartialDerivative<IdxRangeType, DerivativeDim>
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
    virtual std::unique_ptr<IPartialDerivative<IdxRangeType, DerivativeDim>> preallocate()
            const = 0;

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
    DField<IdxRangeType> operator()(DConstField<IdxRangeType> fieldval) const override
    {
        return (*preallocate())(fieldval);
    }
};
