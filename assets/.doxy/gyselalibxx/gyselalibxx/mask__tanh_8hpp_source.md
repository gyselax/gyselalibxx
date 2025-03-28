

# File mask\_tanh.hpp

[**File List**](files.md) **>** [**geometryXVx**](dir_e51b496b46dd687775e46e0826614574.md) **>** [**rhs**](dir_53474cb30a3389ee74cb3186cae99ac0.md) **>** [**mask\_tanh.hpp**](mask__tanh_8hpp.md)

[Go to the documentation of this file](mask__tanh_8hpp.md)


```C++
// SPDX-License-Identifier: MIT
#pragma once
#include "ddc_aliases.hpp"
#include "geometry.hpp"

enum class MaskType { Normal, Inverted };

host_t<DFieldMemX> mask_tanh(
        IdxRangeX const& gridx,
        double extent,
        double stiffness,
        MaskType const type,
        bool normalised);
```


