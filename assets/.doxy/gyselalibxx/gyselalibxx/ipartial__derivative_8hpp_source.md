

# File ipartial\_derivative.hpp

[**File List**](files.md) **>** [**math\_tools**](dir_3ced5d1c6eac490d7704c2e023d148d8.md) **>** [**ipartial\_derivative.hpp**](ipartial__derivative_8hpp.md)

[Go to the documentation of this file](ipartial__derivative_8hpp.md)


```C++
// SPDX-License-Identifier: MIT
#pragma once

#include <ddc/ddc.hpp>

#include "ddc_alias_inline_functions.hpp"
#include "ddc_aliases.hpp"
#include "type_seq_tools.hpp"


template <class IdxRangeFull, class DerivativeDimension>
class IPartialDerivative
{
public:
    using DFieldType = DField<IdxRangeFull>;

    using DConstFieldType = DConstField<IdxRangeFull>;

    using GridDerivativeDimension
            = find_grid_t<DerivativeDimension, ddc::to_type_seq_t<IdxRangeFull>>;

    using IdxRangeDeriv = IdxRange<GridDerivativeDimension>;

    using IdxRangeBatch = ddc::remove_dims_of_t<IdxRangeFull, GridDerivativeDimension>;

    virtual void operator()(DFieldType differentiated_field) const = 0;
};

template <class IdxRangeFull, class DerivativeDimension>
class IPartialDerivativeCreator
{
public:
    virtual std::unique_ptr<IPartialDerivative<IdxRangeFull, DerivativeDimension>> create_instance(
            DConstField<IdxRangeFull> field) const = 0;
};
```


