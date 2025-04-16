

# File central\_fdm\_partial\_derivatives\_with\_boundary\_values.hpp

[**File List**](files.md) **>** [**math\_tools**](dir_3ced5d1c6eac490d7704c2e023d148d8.md) **>** [**central\_fdm\_partial\_derivatives\_with\_boundary\_values.hpp**](central__fdm__partial__derivatives__with__boundary__values_8hpp.md)

[Go to the documentation of this file](central__fdm__partial__derivatives__with__boundary__values_8hpp.md)


```C++
// SPDX-License-Identifier: MIT
#pragma once

#include "ddc/uniform_point_sampling.hpp"

#include "ddc_alias_inline_functions.hpp"
#include "ddc_aliases.hpp"
#include "ipartial_derivative.hpp"


template <class IdxRangeFull, class DerivativeDimension>
class CentralFDMPartialDerivativeWithBValue
    : public IPartialDerivative<IdxRangeFull, DerivativeDimension>
{
private:
    using base_type = IPartialDerivative<IdxRangeFull, DerivativeDimension>;

    using DFieldMemType = DFieldMem<IdxRangeFull>;

    using typename base_type::DFieldType;

    using typename base_type::DConstFieldType;

    using typename base_type::IdxRangeDeriv;

    using typename base_type::IdxRangeBatch;

private:
    DFieldMemType m_field;
    double m_b_value_left;
    double m_b_value_right;

public:
    explicit CentralFDMPartialDerivativeWithBValue(
            DConstFieldType const field_ref,
            double bvalue_left,
            double bvalue_right)
        : m_field(get_idx_range(field_ref))
        , m_b_value_left(bvalue_left)
        , m_b_value_right(bvalue_right)
    {
        ddc::parallel_deepcopy(get_field(m_field), field_ref);
    }

    void operator()(DFieldType differentiated_field) const final
    {
        using IdxFull = typename IdxRangeFull::discrete_element_type;
        using IdxDeriv = typename IdxRangeDeriv::discrete_element_type;
        using IdxBatch = typename IdxRangeBatch::discrete_element_type;
        using IdxStepDeriv = typename IdxRangeDeriv::discrete_vector_type;

        IdxRangeFull idxrange_full = get_idx_range(m_field);
        IdxRangeDeriv idxrange_deriv(idxrange_full);
        IdxRangeBatch idxrange_batch(idxrange_full);

        DConstFieldType const field_proxy = get_const_field(m_field);

        IdxStepDeriv const step(1);

        // front batched derivative
        // we compute second order left and right decentred FDM (using boundary values)
        // and we take the average
        {
            IdxDeriv const ix(idxrange_deriv.front());
            double const h1 = ddc::coordinate(ix + step) - ddc::coordinate(ix);
            double const h2 = ddc::coordinate(ix + 2 * step) - ddc::coordinate(ix + step);
            double const c3 = -h1 / (h2 * (h1 + h2));
            double const c2 = 1. / h1 + 1. / h2;
            double bvalue_left(m_b_value_left);
            ddc::parallel_for_each(
                    Kokkos::DefaultExecutionSpace(),
                    idxrange_batch,
                    KOKKOS_LAMBDA(IdxBatch ib) {
                        IdxFull ibx(ib, ix);
                        double value_left = c2 * field_proxy(ib, ix + step)
                                            + c3 * field_proxy(ib, ix + 2 * step);
                        double value_right = -c2 * bvalue_left - c3 * bvalue_left;
                        differentiated_field(ibx) = (value_right + value_left) / 2;
                    });
        }

        // back batched derivative
        {
            IdxDeriv const ix = idxrange_deriv.back();
            double const h1 = ddc::coordinate(ix) - ddc::coordinate(ix - step);
            double const h2 = ddc::coordinate(ix - step) - ddc::coordinate(ix - 2 * step);
            double const c3 = h1 / (h2 * (h1 + h2));
            double const c2 = -(h1 + h2) / (h1 * h2);
            double bvalue_right(m_b_value_right);
            ddc::parallel_for_each(
                    Kokkos::DefaultExecutionSpace(),
                    idxrange_batch,
                    KOKKOS_LAMBDA(IdxBatch ib) {
                        IdxFull ibx(ib, ix);
                        double value_left = c2 * field_proxy(ib, ix - step)
                                            + c3 * field_proxy(ib, ix - 2 * step);
                        double value_right = -c2 * bvalue_right - c3 * bvalue_right;
                        differentiated_field(ibx) = (value_left + value_right) / 2;
                    });
        }

        // central domain batched derivative
        {
            IdxRangeDeriv idxrange_deriv_central = idxrange_deriv.remove(step, step);
            IdxRangeFull idxrange_central(idxrange_deriv_central, idxrange_batch);
            ddc::parallel_for_each(
                    Kokkos::DefaultExecutionSpace(),
                    idxrange_central,
                    KOKKOS_LAMBDA(IdxFull ibx) {
                        IdxBatch ib(ibx);
                        IdxDeriv ix(ibx);
                        double const h1 = ddc::coordinate(ix) - ddc::coordinate(ix - step);
                        double const h2 = ddc::coordinate(ix + step) - ddc::coordinate(ix);
                        double const c3 = h1 / (h2 * (h1 + h2));
                        double const c2 = 1. / h1 - 1. / h2;
                        double const c1 = -c3 - c2;
                        differentiated_field(ibx) = c1 * field_proxy(ib, ix - step)
                                                    + c2 * field_proxy(ibx)
                                                    + c3 * field_proxy(ib, ix + step);
                    });
        }
    }
};

template <class IdxRangeFull, class DerivativeDimension>
class CentralFDMPartialDerivativeWithBValueCreator
    : public IPartialDerivativeCreator<IdxRangeFull, DerivativeDimension>
{
private:
    using DConstFieldType = DConstField<IdxRangeFull>;

private:
    // The boundary values for the derivative
    double m_b_value_left;
    double m_b_value_right;

public:
    CentralFDMPartialDerivativeWithBValueCreator(double bvalue_left, double bvalue_right)
        : m_b_value_left(bvalue_left)
        , m_b_value_right(bvalue_right)
    {
    }

    std::unique_ptr<IPartialDerivative<IdxRangeFull, DerivativeDimension>> create_instance(
            DConstFieldType field_ref) const final
    {
        return std::make_unique<CentralFDMPartialDerivativeWithBValue<
                IdxRangeFull,
                DerivativeDimension>>(field_ref, m_b_value_left, m_b_value_right);
    }
};
```


