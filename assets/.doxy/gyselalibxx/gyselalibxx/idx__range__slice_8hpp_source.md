

# File idx\_range\_slice.hpp

[**File List**](files.md) **>** [**data\_types**](dir_eaa769653453aaefd8cc10e98e9bb3eb.md) **>** [**idx\_range\_slice.hpp**](idx__range__slice_8hpp.md)

[Go to the documentation of this file](idx__range__slice_8hpp.md)


```C++
// SPDX-License-Identifier: MIT
#pragma once

#include <ddc/ddc.hpp>

#include "ddc_aliases.hpp"


template <class T>
struct IdxRangeToSlice;

template <class... Dims>
struct IdxRangeToSlice<ddc::detail::TypeSeq<Dims...>>
{
    using value = IdxRangeSlice<Dims...>;
};

template <class T>
using to_subidx_range_collection = typename IdxRangeToSlice<T>::value;
```


