#include "nulladvectionvx.h"

DBlockSpanXVx NullAdvectionVx::operator()(DBlockSpanXVx fdistribu, double mass_ratio, double dt)
        const
{
    return fdistribu;
}
