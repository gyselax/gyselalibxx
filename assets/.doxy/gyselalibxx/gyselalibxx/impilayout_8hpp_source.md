

# File impilayout.hpp

[**File List**](files.md) **>** [**mpi\_parallelisation**](dir_a35b8fd75f8fad0c2619b083ab571e51.md) **>** [**impilayout.hpp**](impilayout_8hpp.md)

[Go to the documentation of this file](impilayout_8hpp.md)


```C++
// SPDX-License-Identifier: MIT
#pragma once
#include <ddc/ddc.hpp>

#include "ddc_aliases.hpp"

template <class IdxRangeData, class... DistributedDim>
class IMPILayout
{
    static_assert(ddc::is_discrete_domain_v<IdxRangeData>);

public:
    using discrete_domain_type = IdxRangeData;
    using distributed_sub_idx_range = IdxRange<DistributedDim...>;
    using distributed_type_seq = ddc::detail::TypeSeq<DistributedDim...>;

private:
    using type_seq = typename ddc::to_type_seq_t<discrete_domain_type>;

    template <std::size_t... I>
    static constexpr bool check_distributed_idx_ranges_are_first(
            std::integer_sequence<std::size_t, I...>)
    {
        return (ddc::in_tags_v<ddc::type_seq_element_t<I, type_seq>, distributed_type_seq> && ...);
    }

public:
    static constexpr int n_distributed_dimensions = sizeof...(DistributedDim);
    static constexpr bool distributed_idx_ranges_are_first = check_distributed_idx_ranges_are_first(
            std::make_index_sequence<n_distributed_dimensions>());
};
```


