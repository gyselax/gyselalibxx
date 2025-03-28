

# File chargedensitycalculator.hpp

[**File List**](files.md) **>** [**geometryXYVxVy**](dir_e4674dab6493cf35bbeb1b23e7fbbddd.md) **>** [**poisson**](dir_14c5eb4d397dfd4e1a4d5c7bede9e118.md) **>** [**chargedensitycalculator.hpp**](geometryXYVxVy_2poisson_2chargedensitycalculator_8hpp.md)

[Go to the documentation of this file](geometryXYVxVy_2poisson_2chargedensitycalculator_8hpp.md)


```C++
// SPDX-License-Identifier: MIT

#pragma once

#include <ddc/ddc.hpp>

#include "ddc_helper.hpp"
#include "geometry.hpp"
#include "ichargedensitycalculator.hpp"
#include "quadrature.hpp"

class ChargeDensityCalculator : public IChargeDensityCalculator
{
private:
    Quadrature<IdxRangeVxVy, IdxRangeXYVxVy> m_quadrature;

public:
    explicit ChargeDensityCalculator(DConstFieldVxVy coeffs);

    void operator()(DFieldXY rho, DConstFieldSpVxVyXY allfdistribu) const final;
};
```


