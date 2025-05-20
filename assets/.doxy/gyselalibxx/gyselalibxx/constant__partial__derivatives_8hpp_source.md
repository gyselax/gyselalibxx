

# File constant\_partial\_derivatives.hpp

[**File List**](files.md) **>** [**math\_tools**](dir_3ced5d1c6eac490d7704c2e023d148d8.md) **>** [**constant\_partial\_derivatives.hpp**](constant__partial__derivatives_8hpp.md)

[Go to the documentation of this file](constant__partial__derivatives_8hpp.md)


```C++
// SPDX-License-Identifier: MIT
#pragma once

template <class IdxRangeFull, class DerivativeDimension>
class ConstantPartialDerivative : public IPartialDerivative<IdxRangeFull, DerivativeDimension>
{
private:
    using base_type = IPartialDerivative<IdxRangeFull, DerivativeDimension>;

    using typename base_type::DFieldType;

private:
    double m_deriv_value;

public:
    explicit ConstantPartialDerivative(double deriv_value) : m_deriv_value(deriv_value) {}

    void operator()(DFieldType differentiated_field) const final
    {
        ddc::parallel_fill(differentiated_field, m_deriv_value);
    }
};

template <class IdxRangeFull, class DerivativeDimension>
class ConstantPartialDerivativeCreator
    : public IPartialDerivativeCreator<IdxRangeFull, DerivativeDimension>
{
private:
    double m_deriv_value;

public:
    explicit ConstantPartialDerivativeCreator(double deriv_value) : m_deriv_value(deriv_value) {}

    std::unique_ptr<IPartialDerivative<IdxRangeFull, DerivativeDimension>> create_instance(
            DConstField<IdxRangeFull> field) const final
    {
        return std::make_unique<ConstantPartialDerivative<IdxRangeFull, DerivativeDimension>>(
                m_deriv_value);
    }
};
```


