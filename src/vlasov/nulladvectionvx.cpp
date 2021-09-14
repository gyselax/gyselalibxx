#include "nulladvectionvx.h"

DSpanXVx NullAdvectionVx::operator()(DSpanXVx fdistribu, double mass_ratio, double dt) const
{
    return fdistribu;
}
