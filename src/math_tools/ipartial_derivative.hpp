// SPDX-License-Identifier: MIT
#pragma once

#include <ddc/ddc.hpp>

#include "ddc_alias_inline_functions.hpp"
#include "ddc_aliases.hpp"
#include "type_seq_tools.hpp"


/**
 * @brief An abstract class for a partial derivative operator.
 *
 * @tparam IdxRangeFull The index range of the field on which the operator acts 
 * (with all dimensions, batched and dimension of interest).
 * @tparam DerivativeDimension The dimension on which the partial derivative is calculated.
 */
template <class IdxRangeFull, class DerivativeDimension>
class IPartialDerivative
{
public:
    /// The type of a reference to the field to be differentiated.
    using DFieldType = DField<IdxRangeFull>;

    /// The type of a constant reference to the field to be differentiated.
    using DConstFieldType = DConstField<IdxRangeFull>;

    /// The type of the grid on the dimension on which the partial derivative is calculated.
    using GridDerivativeDimension
            = find_grid_t<DerivativeDimension, ddc::to_type_seq_t<IdxRangeFull>>;

    /// The index range of the dimension on which the partial derivative is calculated.
    using IdxRangeDeriv = IdxRange<GridDerivativeDimension>;

    /// The index range of all dimensions except DerivativeDimension.
    using IdxRangeBatch = ddc::remove_dims_of_t<IdxRangeFull, GridDerivativeDimension>;

    /**
     * @brief Compute the partial derivative of a field in a given direction.
     *
     * @param[out] differentiated_field On output, contains values of the differentiated field.
     */
    virtual void operator()(DFieldType differentiated_field) const = 0;
};

/**
 * @brief An abstract class which provides a create_instance function to instantiate 
 * an object of the IPartialDerivative class where required. 
 * 
 */
template <class IdxRangeFull, class DerivativeDimension>
class IPartialDerivativeCreator
{
public:
    /**
     * @brief Create an instance of a pointer to an IPartialDerivative object.
     *
     * @param[in] field A field to be passed to the constructor of IPartialDerivative. 
     *
     * @return A pointer to an IPartialDerivative object.
     *
     * @see IPartialDerivative
     */
    virtual std::unique_ptr<IPartialDerivative<IdxRangeFull, DerivativeDimension>> create_instance(
            DConstField<IdxRangeFull> field) const = 0;
};
