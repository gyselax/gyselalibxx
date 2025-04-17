

# File poisson\_like\_rhs\_function.hpp

[**File List**](files.md) **>** [**geometryRTheta**](dir_e9f169004bcfe9f3cb1f8a27ce024e59.md) **>** [**poisson**](dir_131fdd0509f46f459997bddabd4481b1.md) **>** [**poisson\_like\_rhs\_function.hpp**](poisson__like__rhs__function_8hpp.md)

[Go to the documentation of this file](poisson__like__rhs__function_8hpp.md)


```C++
// SPDX-License-Identifier: MIT
#pragma once

#include <ddc/ddc.hpp>

#include "geometry.hpp"



template <class RadialExtrapolationRule>
class PoissonLikeRHSFunction
{
public:
    using evaluator_type = ddc::SplineEvaluator2D<
            Kokkos::DefaultHostExecutionSpace,
            Kokkos::HostSpace,
            BSplinesR,
            BSplinesTheta,
            GridR,
            GridTheta,
            RadialExtrapolationRule,
            RadialExtrapolationRule,
            ddc::PeriodicExtrapolationRule<Theta>,
            ddc::PeriodicExtrapolationRule<Theta>>;

private:
    host_t<ConstSpline2D> const m_coefs;
    evaluator_type const& m_evaluator;

public:
    PoissonLikeRHSFunction(host_t<ConstSpline2D> coefs, evaluator_type const& evaluator)
        : m_coefs(coefs)
        , m_evaluator(evaluator)
    {
    }

    double operator()(CoordRTheta const& coord_rtheta) const
    {
        return m_evaluator(coord_rtheta, m_coefs);
    }
};
```


