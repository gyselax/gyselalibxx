#pragma once

#include <string>
#include <type_traits>

#include <pdi.h>

#include "ddc/block_span.hpp"

template <class SupportType, class ElementType, class LayoutStridedPolicy>
void expose_to_pdi(
        std::string const& name,
        BlockSpan<SupportType, ElementType, LayoutStridedPolicy> const& data)
{
    auto dom = data.domain();
    auto extents = dom.extents().array();
    PDI_expose((name + "_extents").c_str(), extents.data(), PDI_OUT);
    if constexpr (std::is_const_v<ElementType>) {
        PDI_multi_expose(
                name.c_str(),
                (name + "_extents").c_str(),
                extents.data(),
                PDI_OUT,
                name.c_str(),
                const_cast<typename std::remove_const<ElementType>::type*>(data.data()),
                PDI_OUT,
                NULL);
    } else {
        PDI_expose(name.c_str(), data.data(), PDI_INOUT);
    }
}
