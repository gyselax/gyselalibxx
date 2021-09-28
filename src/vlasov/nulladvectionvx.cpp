#include "nulladvectionvx.h"

DSpanXVx NullAdvectionVx::operator()(
        DSpanXVx fdistribu,
        DViewX efield,
        double sqrt_me_on_mspecies,
        double dt) const
{
    return fdistribu;
}
