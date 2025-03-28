

# File chargedensitycalculator.cpp

[**File List**](files.md) **>** [**geometryXVx**](dir_e51b496b46dd687775e46e0826614574.md) **>** [**poisson**](dir_d78fdb6d05340e24a2e187de33ea09a4.md) **>** [**chargedensitycalculator.cpp**](geometryXVx_2poisson_2chargedensitycalculator_8cpp.md)

[Go to the documentation of this file](geometryXVx_2poisson_2chargedensitycalculator_8cpp.md)


```C++
// SPDX-License-Identifier: MIT

#include <ddc/ddc.hpp>

#include "chargedensitycalculator.hpp"

ChargeDensityCalculator::ChargeDensityCalculator(DConstFieldVx coeffs) : m_quadrature(coeffs) {}

DFieldX ChargeDensityCalculator::operator()(DFieldX const rho, DConstFieldSpXVx const allfdistribu)
        const
{
    Kokkos::Profiling::pushRegion("ChargeDensityCalculator");

    IdxRangeSp const kin_species_idx_range = get_idx_range<Species>(allfdistribu);

    m_quadrature(
            Kokkos::DefaultExecutionSpace(),
            rho,
            KOKKOS_LAMBDA(IdxXVx ixvx) {
                double sum = 0.0;
                for (IdxSp isp : kin_species_idx_range) {
                    sum += charge(isp) * allfdistribu(isp, ixvx);
                }
                return sum;
            });

    IdxSp const last_kin_species = kin_species_idx_range.back();
    IdxSp const last_species = kin_species_idx_range.back();
    if (last_kin_species != last_species) {
        double const chargedens_adiabspecies = charge(last_species);

        ddc::parallel_for_each(
                Kokkos::DefaultExecutionSpace(),
                get_idx_range(rho),
                KOKKOS_LAMBDA(IdxX ix) { rho(ix) += chargedens_adiabspecies; });
    }

    Kokkos::Profiling::popRegion();

    return rho;
}
```


