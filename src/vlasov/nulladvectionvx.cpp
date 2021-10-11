#include "nulladvectionvx.h"

DSpanXVx NullAdvectionVx::operator()(
        DSpanXVx fdistribu,
        DViewX efield,
        int charge_species,
        double sqrt_me_on_mspecies,
        double dt) const
{
    return fdistribu;
}
