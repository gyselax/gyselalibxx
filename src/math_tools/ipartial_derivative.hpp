// SPDX-License-Identifier: MIT
#pragma once

#include <ddc/ddc.hpp>

#include "ddc_alias_inline_functions.hpp"
#include "ddc_aliases.hpp"


/**
 * @brief An abstract class for a partial derivative operator.
 * @tparam IdxRangeType The index range of the field on which the operator acts.
 * @tparam DerivativeDimension The dimension on which the partial derivative is calculated.
 */
template <class IdxRangeType, class DerivativeDimension>
class IPartialDerivative
{
public:
    using DFieldMemType = DFieldMem<IdxRangeType>;
    using DFieldType = DField<IdxRangeType>;
    using DConstFieldType = DConstField<IdxRangeType>;
    ;

public:
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
template <class IdxRangeType, class DerivativeDimension>
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
    virtual std::unique_ptr<IPartialDerivative<IdxRangeType, DerivativeDimension>> create_instance(
            DConstField<IdxRangeType> field) const = 0;
};
