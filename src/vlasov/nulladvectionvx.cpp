#include "nulladvectionvx.h"

DBlockViewXVx NullAdvectionVx::operator()(DBlockViewXVx fdistribu, double mass_ratio, double dt)
        const
{
    return fdistribu;
}
