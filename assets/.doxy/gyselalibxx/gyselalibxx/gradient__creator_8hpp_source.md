

# File gradient\_creator.hpp

[**File List**](files.md) **>** [**math\_tools**](dir_3ced5d1c6eac490d7704c2e023d148d8.md) **>** [**gradient\_creator.hpp**](gradient__creator_8hpp.md)

[Go to the documentation of this file](gradient__creator_8hpp.md)


```C++
// SPDX-License-Identifier: MIT
#pragma once
#include <tuple>

#include "coord_transformation_tools.hpp"
#include "ipartial_derivative.hpp"
#include "vector_mapper.hpp"

template <typename IdxRangeFull, typename... DerivativeDims>
class GradientCreator
{
    static_assert((DerivativeDims::IS_CONTRAVARIANT && ...));

private:
    std::tuple<IPartialDerivativeCreator<IdxRangeFull, DerivativeDims> const&...>
            m_derivative_creators;

public:
    explicit GradientCreator(IPartialDerivativeCreator<
                             IdxRangeFull,
                             DerivativeDims> const&... partial_derivative_operator)
        : m_derivative_creators(std::tie(partial_derivative_operator...))
    {
    }

    void operator()(
            DVectorField<IdxRangeFull, get_covariant_dims_t<VectorIndexSet<DerivativeDims...>>>
                    grad_func_cov,
            DConstField<IdxRangeFull> func) const
    {
        (((*std::get<IPartialDerivativeCreator<IdxRangeFull, DerivativeDims> const&>(
                    m_derivative_creators)
                    .create_instance(get_const_field(func)))(
                 ddcHelper::get<typename DerivativeDims::Dual>(grad_func_cov))),
         ...);
    }
};
```


