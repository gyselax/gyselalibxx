#include <cmath>
#include <iostream>

#include "species_info.hpp"

using std::sqrt, std::exp;

void maxwellian_initialization(
        const double density,
        const double temperature,
        const double mean_velocity,
        DSpanVx fMaxwellian)
{
    const double inv_sqrt_2piT = 1. / sqrt(2. * M_PI * temperature);
    auto gridvx = fMaxwellian.domain<IDimVx>();
    for (IndexVx iv : gridvx) {
        const CoordVx v = gridvx.to_real(iv);
        fMaxwellian(iv) = density * inv_sqrt_2piT
                          * exp(-(v - mean_velocity) * (v - mean_velocity) / (2. * temperature));
    }
}
