

# File constant\_identity\_interpolation\_extrapolation\_rule.hpp

[**File List**](files.md) **>** [**interpolation**](dir_264890e5c091f8c8d7fe1f842870c25e.md) **>** [**constant\_identity\_interpolation\_extrapolation\_rule.hpp**](constant__identity__interpolation__extrapolation__rule_8hpp.md)

[Go to the documentation of this file](constant__identity__interpolation__extrapolation__rule_8hpp.md)


```C++
// SPDX-License-Identifier: MIT
#pragma once
#include "ddc_aliases.hpp"

template <class CoeffGrid, class DataType = double>
class ConstantIdentityInterpolationExtrapolationRule
{
private:
    Idx<CoeffGrid> m_coeff_idx;

public:
    explicit ConstantIdentityInterpolationExtrapolationRule(Idx<CoeffGrid> coeff_idx)
        : m_coeff_idx(coeff_idx)
    {
    }

    template <class CoordType, class Layout, class MemorySpace>
    KOKKOS_FUNCTION double operator()(
            [[maybe_unused]] CoordType pos,
            ConstField<DataType, IdxRange<CoeffGrid>, Layout, MemorySpace> const interp_coef) const
    {
        return interp_coef(m_coeff_idx);
    }
};
```


