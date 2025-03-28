

# File spline\_1d\_partial\_derivative.hpp

[**File List**](files.md) **>** [**math\_tools**](dir_3ced5d1c6eac490d7704c2e023d148d8.md) **>** [**spline\_1d\_partial\_derivative.hpp**](spline__1d__partial__derivative_8hpp.md)

[Go to the documentation of this file](spline__1d__partial__derivative_8hpp.md)


```C++
// SPDX-License-Identifier: MIT

#pragma once

#include "ddc_alias_inline_functions.hpp"
#include "ddc_aliases.hpp"
#include "ipartial_derivative.hpp"


template <class Spline1DBuilder, class Spline1DEvaluator>
class Spline1DPartialDerivative
    : public IPartialDerivative<
              typename Spline1DBuilder::batched_interpolation_domain_type,
              typename Spline1DBuilder::continuous_dimension_type>
{
    static_assert(std::is_same_v<
                  typename Spline1DBuilder::batched_spline_domain_type,
                  typename Spline1DEvaluator::batched_spline_domain_type>);
    static_assert(std::is_same_v<
                  typename Spline1DBuilder::batched_interpolation_domain_type,
                  typename Spline1DEvaluator::batched_evaluation_domain_type>);

private:
    using base_type = IPartialDerivative<
            typename Spline1DBuilder::batched_interpolation_domain_type,
            typename Spline1DBuilder::continuous_dimension_type>;

    using typename base_type::DConstFieldType;
    using typename base_type::DFieldType;

    using IdxRangeBS = typename Spline1DBuilder::batched_spline_domain_type;
    using DFieldBSMem = DFieldMem<IdxRangeBS>;
    using DFieldBS = DField<IdxRangeBS>;

    Spline1DBuilder const& m_builder;
    Spline1DEvaluator const& m_evaluator;
    DFieldBSMem m_spline_coefs;

public:
    explicit Spline1DPartialDerivative(
            Spline1DBuilder const& builder,
            Spline1DEvaluator const& evaluator,
            DConstFieldType const field)
        : m_builder(builder)
        , m_evaluator(evaluator)
        , m_spline_coefs(builder.batched_spline_domain())
    {
        m_builder(get_field(m_spline_coefs), field);
    }

    void operator()(DFieldType differentiated_field) const final
    {
        m_evaluator.deriv(differentiated_field, get_const_field(m_spline_coefs));
    }
};


template <class Spline1DBuilder, class Spline1DEvaluator>
class Spline1DPartialDerivativeCreator
    : public IPartialDerivativeCreator<
              typename Spline1DBuilder::batched_interpolation_domain_type,
              typename Spline1DBuilder::continuous_dimension_type>
{
private:
    using DConstFieldType
            = DConstField<typename Spline1DBuilder::batched_interpolation_domain_type>;

    Spline1DBuilder const& m_builder;
    Spline1DEvaluator const& m_evaluator;

public:
    Spline1DPartialDerivativeCreator(
            Spline1DBuilder const& builder,
            Spline1DEvaluator const& evaluator)
        : m_builder(builder)
        , m_evaluator(evaluator)
    {
    }

    std::unique_ptr<IPartialDerivative<
            typename Spline1DBuilder::batched_interpolation_domain_type,
            typename Spline1DBuilder::continuous_dimension_type>>
    create_instance(DConstFieldType field) const final
    {
        return std::make_unique<Spline1DPartialDerivative<
                Spline1DBuilder,
                Spline1DEvaluator>>(m_builder, m_evaluator, field);
    }
};
```


