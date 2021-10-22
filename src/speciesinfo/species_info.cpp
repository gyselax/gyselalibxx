#include <cmath>
#include <iostream>

#include "species_info.hpp"

using std::sqrt, std::exp;

void maxwellian_initialization(
        double const density,
        double const temperature,
        double const mean_velocity,
        DSpanVx fMaxwellian)
{
    double const inv_sqrt_2piT = 1. / sqrt(2. * M_PI * temperature);
    IDomainVx gridvx = fMaxwellian.domain<IDimVx>();
    for (IndexVx iv : gridvx) {
        CoordVx const v = gridvx.to_real(iv);
        fMaxwellian(iv) = density * inv_sqrt_2piT
                          * exp(-(v - mean_velocity) * (v - mean_velocity) / (2. * temperature));
    }
}
