

# File multipatch\_lagrange\_interpolation.hpp

[**File List**](files.md) **>** [**interpolation**](dir_54175eb71bf3bdb1ca864f2a0d352472.md) **>** [**multipatch\_lagrange\_interpolation.hpp**](multipatch__lagrange__interpolation_8hpp.md)

[Go to the documentation of this file](multipatch__lagrange__interpolation_8hpp.md)


```C++
// SPDX-License-Identifier: MIT

#pragma once
#include "ddc_aliases.hpp"
#include "lagrange_interpolation.hpp"
#include "multipatch_interpolation_2d.hpp"

namespace detail {

template <
        class ExecSpace,
        class DataType,
        template <class>
        class Basis1OnPatch,
        template <class>
        class Basis2OnPatch>
struct LagrangeInterpolatorOnPatch
{
    template <class Patch>
    using type = ::LagrangeInterpolator<
            ExecSpace,
            DataType,
            IdxRange<Basis1OnPatch<Patch>, Basis2OnPatch<Patch>>,
            typename Patch::IdxRange12,
            PatchExtrapolationRules<typename Patch::Dim1>,
            PatchExtrapolationRules<typename Patch::Dim2>>;
};

} // namespace detail

template <
        class PatchLocator,
        class MultipatchExtrapolationRule,
        template <class>
        class Basis1OnPatch,
        template <class>
        class Basis2OnPatch,
        class DataType = double>
using MultipatchLagrangeInterpolator2D = MultipatchInterpolation2D<
        detail::LagrangeInterpolatorOnPatch<
                typename PatchLocator::exec_space,
                DataType,
                Basis1OnPatch,
                Basis2OnPatch>::template type,
        PatchLocator,
        MultipatchExtrapolationRule>;
```


