// SPDX-License-Identifier: MIT

#include "nullpoissonsolver.hpp"

void NullPoissonSolver::operator()(
        device_t<DSpanX> const electrostatic_potential_device,
        device_t<DSpanX> const electric_field_device,
        device_t<DViewSpXVx> const allfdistribu_device) const
{
}
