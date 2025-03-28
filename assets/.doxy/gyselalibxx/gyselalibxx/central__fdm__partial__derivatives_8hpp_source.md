

# File central\_fdm\_partial\_derivatives.hpp

[**File List**](files.md) **>** [**math\_tools**](dir_3ced5d1c6eac490d7704c2e023d148d8.md) **>** [**central\_fdm\_partial\_derivatives.hpp**](central__fdm__partial__derivatives_8hpp.md)

[Go to the documentation of this file](central__fdm__partial__derivatives_8hpp.md)


```C++
// SPDX-License-Identifier: MIT
#pragma once

#include <iostream>

#include "ddc/discrete_domain.hpp"
#include "ddc/uniform_point_sampling.hpp"

#include "ddc_alias_inline_functions.hpp"
#include "ddc_aliases.hpp"
#include "ipartial_derivative.hpp"


template <class IdxRangeFull, class DerivativeDimension>
class CentralFDMPartialDerivative : public IPartialDerivative<IdxRangeFull, DerivativeDimension>
{
private:
    using base_type = IPartialDerivative<IdxRangeFull, DerivativeDimension>;

    using DFieldMemType = DFieldMem<IdxRangeFull>;

    using typename base_type::DFieldType;

    using typename base_type::DConstFieldType;

    using typename base_type::IdxRangeDeriv;

    using typename base_type::IdxRangeBatch;

    DFieldMemType m_field;

public:
    explicit CentralFDMPartialDerivative(DConstFieldType const field_ref)
        : m_field(get_idx_range(field_ref))
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
        {
            IdxDeriv const ix(idxrange_deriv.front());
            double const h1 = ddc::coordinate(ix + step) - ddc::coordinate(ix);
            double const h2 = ddc::coordinate(ix + 2 * step) - ddc::coordinate(ix + step);
            double const c3 = -h1 / (h2 * (h1 + h2));
            double const c2 = 1. / h1 + 1. / h2;
            double const c1 = -c3 - c2;
            ddc::parallel_for_each(
                    Kokkos::DefaultExecutionSpace(),
                    idxrange_batch,
                    KOKKOS_LAMBDA(IdxBatch ib) {
                        IdxFull ibx(ib, ix);
                        differentiated_field(ibx) = c1 * field_proxy(ibx)
                                                    + c2 * field_proxy(ib, ix + step)
                                                    + c3 * field_proxy(ib, ix + 2 * step);
                    });
        }

        // back batched derivative
        {
            IdxDeriv const ix = idxrange_deriv.back();
            double const h1 = ddc::coordinate(ix) - ddc::coordinate(ix - step);
            double const h2 = ddc::coordinate(ix - step) - ddc::coordinate(ix - 2 * step);
            double const c3 = h1 / (h2 * (h1 + h2));
            double const c2 = -(h1 + h2) / (h1 * h2);
            double const c1 = -c3 - c2;
            ddc::parallel_for_each(
                    Kokkos::DefaultExecutionSpace(),
                    idxrange_batch,
                    KOKKOS_LAMBDA(IdxBatch ib) {
                        IdxFull ibx(ib, ix);
                        differentiated_field(ibx) = c1 * field_proxy(ibx)
                                                    + c2 * field_proxy(ib, ix - step)
                                                    + c3 * field_proxy(ib, ix - 2 * step);
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
class CentralFDMPartialDerivativeCreator
    : public IPartialDerivativeCreator<IdxRangeFull, DerivativeDimension>
{
private:
    using DConstFieldType = DConstField<IdxRangeFull>;

public:
    CentralFDMPartialDerivativeCreator() = default;

    std::unique_ptr<IPartialDerivative<IdxRangeFull, DerivativeDimension>> create_instance(
            DConstFieldType field_ref) const final
    {
        return std::make_unique<CentralFDMPartialDerivative<IdxRangeFull, DerivativeDimension>>(
                field_ref);
    }
};
```


