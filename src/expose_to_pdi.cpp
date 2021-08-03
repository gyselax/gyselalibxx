#include <array>

#include <pdi.h>

#include "expose_to_pdi.h"
#include "mcoord.h"
#include "product_mdomain.h"

void expose_to_pdi(std::string const& name, DBlockSpanXVx b)
{
    auto dom = b.domain();
    auto extents = dom.extents().array();
    PDI_expose((name + "_extents").c_str(), extents.data(), PDI_OUT);
    PDI_expose(name.c_str(), b.data(), PDI_OUT);
}
