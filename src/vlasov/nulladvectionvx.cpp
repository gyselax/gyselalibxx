#include "nulladvectionvx.h"

DSpanSpXVx NullAdvectionVx::operator()(DSpanSpXVx fdistribu, DViewX electric_potential, double dt)
        const
{
    return fdistribu;
}
