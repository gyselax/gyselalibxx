

# File moments.hpp

[**File List**](files.md) **>** [**speciesinfo**](dir_661be8452a62f1b4720eb6eb57123ae7.md) **>** [**moments.hpp**](moments_8hpp.md)

[Go to the documentation of this file](moments_8hpp.md)


```C++
// SPDX-License-Identifier: MIT

#pragma once

class Moments
{
public:
    using discrete_dimension_type = Moments;

public:
    template <class Grid1D, class MemorySpace>
    class Impl
    {
        template <class ODDim, class OMemorySpace>
        friend class Impl;

    public:
        using discrete_dimension_type = Moments;

        template <class OMemorySpace>
        explicit Impl(Impl<Grid1D, OMemorySpace> const& impl)
        {
        }

        Impl() {}
    };
};
```


