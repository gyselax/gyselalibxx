#include <array>

#include <ddc/MCoord>
#include <ddc/ProductMDomain>

#include <pdi.h>

#include "expose_to_pdi.h"

void expose_to_pdi(std::string const& name, DSpanXVx b)
{
    auto dom = b.domain();
    auto extents = dom.extents().array();
    PDI_expose((name + "_extents").c_str(), extents.data(), PDI_OUT);
    PDI_expose(name.c_str(), b.data(), PDI_OUT);
}
