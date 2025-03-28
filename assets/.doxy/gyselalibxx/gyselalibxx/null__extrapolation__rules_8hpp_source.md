

# File null\_extrapolation\_rules.hpp

[**File List**](files.md) **>** [**multipatch**](dir_7740c6927b2da0a836b00bedb040a06d.md) **>** [**spline**](dir_729d943c83b6b5573a69e28a4db4673a.md) **>** [**null\_extrapolation\_rules.hpp**](null__extrapolation__rules_8hpp.md)

[Go to the documentation of this file](null__extrapolation__rules_8hpp.md)


```C++
// SPDX-License-Identifier: MIT

#pragma once
#include <stdexcept>

#include <ddc/ddc.hpp>
#include <ddc/kernels/splines.hpp>

#include "ddc_aliases.hpp"
#include "multipatch_type.hpp"

struct NullExtrapolationRule
{
    template <class... Dim, template <typename P> typename SplinesOnPatch, class... Patches>
    KOKKOS_FUNCTION double operator()(
            Coord<Dim...> const& coord_extrap,
            MultipatchType<SplinesOnPatch, Patches...> const& patches_splines,
            int const out_of_bounds_idx) const
    {
        assert(out_of_bounds_idx < 0);
        return 0.0;
    }
};
```


