#include "nulladvectionvx.h"

DSpanSpXVx NullAdvectionVx::operator()(DSpanSpXVx fdistribu, DViewX efield, double dt) const
{
    return fdistribu;
}
