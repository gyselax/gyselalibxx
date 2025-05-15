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
class ConstantPartialDerivative : public IPartialDerivative<IdxRangeFull, DerivativeDimension>
{
private:
    using base_type = IPartialDerivative<IdxRangeFull, DerivativeDimension>;

    /// The type of a reference to the field to be differentiated.
    using typename base_type::DFieldType;

private:
    double m_deriv_value;

public:
    /**
     * @brief Create an instance of ConstantPartialDerivative.
     *
     * @param[in] deriv_value The value that should be returned as the constant value of
     *          the derivative.
     */
    ConstantPartialDerivative(double deriv_value) : m_deriv_value(deriv_value) {}

    /**
     * @brief Set the partial derivative of a field to 0.
     *
     * @param[out] differentiated_field On output, contains values of the differentiated field.
     */
    void operator()(DFieldType differentiated_field) const final
    {
        ddc::parallel_fill(differentiated_field, m_deriv_value);
    }
};

/**
 * @brief A class to create a ConstantPartialDerivative via a create_instance function.
 *
 * @tparam IdxRangeFull The index range of the field on which the operator acts 
 * (with all dimensions, batched and dimension of interest, used for inheritance).
 * @tparam DerivativeDimension The dimension on which the partial derivative is calculated.
 * (used for inheritance).
 */
template <class IdxRangeFull, class DerivativeDimension>
class ConstantPartialDerivativeCreator
    : public IPartialDerivativeCreator<IdxRangeFull, DerivativeDimension>
{
private:
    double m_deriv_value;

public:
    /**
     * @brief Create an instance of ConstantPartialDerivativeCreator.
     *
     * @param[in] deriv_value The value that should be returned as the constant value of
     *          the derivative.
     */
    ConstantPartialDerivativeCreator(double deriv_value) : m_deriv_value(deriv_value) {}

    /**
     * @brief Create an instance of a pointer to an IPartialDerivative object.
     *
     * @param[in] field The field whose derivative should be calculated.
     *
     * @return A pointer to an IPartialDerivative object.
     */
    std::unique_ptr<IPartialDerivative<IdxRangeFull, DerivativeDimension>> create_instance(
            DConstField<IdxRangeFull> field) const final
    {
        return std::make_unique<ConstantPartialDerivative<IdxRangeFull, DerivativeDimension>>(
                m_deriv_value);
    }
};
