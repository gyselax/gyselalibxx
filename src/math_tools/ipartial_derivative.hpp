// SPDX-License-Identifier: MIT
#pragma once

#include <ddc/ddc.hpp>

#include "ddc_alias_inline_functions.hpp"
#include "ddc_aliases.hpp"


/**
 * @brief An abstract class for a partial derivative operator.
 * @tparam IdxRangeBatched The index range of the field on which the operator acts 
 * (with all dimensions, including batched).
 * @tparam DerivativeDimension The dimension on which the partial derivative is calculated.
 */
template <class IdxRangeBatched, class DerivativeDimension>
class IPartialDerivative
{
public:
    /// The type of a reference to the field to be differentiated.
    using DFieldType = DField<IdxRangeBatched>;

    /// The type of a constant reference to the field to be differentiated.
    using DConstFieldType = DConstField<IdxRangeBatched>;

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
template <class IdxRangeBatched, class DerivativeDimension>
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
    virtual std::unique_ptr<IPartialDerivative<IdxRangeBatched, DerivativeDimension>>
    create_instance(DConstField<IdxRangeBatched> field) const = 0;
};
