#include "nulladvectionvx.h"

DSpanSpXVx NullAdvectionVx::operator()(
        DSpanSpXVx allfdistribu,
        DViewX electric_potential,
        double dt) const
{
    return allfdistribu;
}
