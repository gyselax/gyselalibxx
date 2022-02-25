// SPDX-License-Identifier: MIT

#include "nulladvectionvx.hpp"

DSpanSpXVx NullAdvectionVx::operator()(DSpanSpXVx const allfdistribu, DViewX const, double const)
        const
{
    return allfdistribu;
}
