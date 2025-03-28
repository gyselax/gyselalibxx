

# File ivlasovsolver.hpp

[**File List**](files.md) **>** [**geometryXYVxVy**](dir_e4674dab6493cf35bbeb1b23e7fbbddd.md) **>** [**vlasov**](dir_0a9688649b1824bbfb2c211b845ba732.md) **>** [**ivlasovsolver.hpp**](ivlasovsolver_8hpp.md)

[Go to the documentation of this file](ivlasovsolver_8hpp.md)


```C++
// SPDX-License-Identifier: MIT

#pragma once

#include "geometry.hpp"

class IVlasovSolver
{
public:
    virtual ~IVlasovSolver() = default;

    virtual DFieldSpVxVyXY operator()(
            DFieldSpVxVyXY allfdistribu,
            DConstFieldXY efield_x,
            DConstFieldXY efield_y,
            double dt) const = 0;
};
```


