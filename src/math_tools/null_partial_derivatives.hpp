// SPDX-License-Identifier: MIT
#pragma once

/**
 * @brief A class to get the derivative of a constant function.
 * When the derivative of a function is known to be 0 but the dimension
 * is still needed this class can be used to avoid unnecessary calculations.
 * @tparam IdxRangeFull The index range of the field on which the operator acts 
 * (with all dimensions, batched and dimension of interest, used for inheritance).
 * @tparam DerivativeDimension The dimension on which the partial derivative is calculated.
 * (used for inheritance).
 */
template <class IdxRangeFull, class DerivativeDimension>
class NullPartialDerivative : public IPartialDerivative<IdxRangeFull, DerivativeDimension>
{
private:
    using base_type = IPartialDerivative<IdxRangeFull, DerivativeDimension>;

    /// The type of a reference to the field to be differentiated.
    using typename base_type::DFieldType;

public:
    /**
     * @brief Set the partial derivative of a field to 0.
     *
     * @param[out] differentiated_field On output, contains values of the differentiated field.
     */
    void operator()(DFieldType differentiated_field) const final
    {
        ddc::parallel_fill(differentiated_field, 0.0);
    }
};

/**
 * @brief A class to create a NullPartialDerivative via a create_instance function.
 *
 * @tparam IdxRangeFull The index range of the field on which the operator acts 
 * (with all dimensions, batched and dimension of interest, used for inheritance).
 * @tparam DerivativeDimension The dimension on which the partial derivative is calculated.
 * (used for inheritance).
 */
template <class IdxRangeFull, class DerivativeDimension>
class NullPartialDerivativeCreator
    : public IPartialDerivativeCreator<IdxRangeFull, DerivativeDimension>
{
public:
    /**
     * @brief Create an instance of a pointer to an IPartialDerivative object.
     *
     * @param[in] field The field whose derivative should be calculated.
     *
     * @return A pointer to an IPartialDerivative object.
     */
    std::unique_ptr<IPartialDerivative<IdxRangeFull, DerivativeDimension>> create_instance(
            DConstField<IdxRangeFull> field_ref) const final
    {
        return std::make_unique<NullPartialDerivative<IdxRangeFull, DerivativeDimension>>();
    }
};
