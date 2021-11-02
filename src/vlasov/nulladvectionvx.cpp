#include "nulladvectionvx.hpp"

DSpanSpXVx NullAdvectionVx::operator()(
        DSpanSpXVx const allfdistribu,
        DViewX const electric_potential,
        double const dt) const
{
    return allfdistribu;
}
