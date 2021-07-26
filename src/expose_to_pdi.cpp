#include <pdi.h>

#include "expose_to_pdi.h"

void expose_to_pdi(std::string const& name, DBlockSpanXVx b)
{
    auto dom = b.domain();
    auto extents = dom.extents().array();
    PDI_expose((name + "_extents").c_str(), extents.data(), PDI_OUT);
    PDI_expose(name.c_str(), b.raw_view().data(), PDI_OUT);
}
