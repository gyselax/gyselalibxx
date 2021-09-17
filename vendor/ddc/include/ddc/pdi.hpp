#pragma once

#include <string>
#include <type_traits>

#include <pdi.h>

#include "ddc/block_span.hpp"

template <class ElementType>
static inline constexpr PDI_inout_t default_access()
{
    if constexpr (std::is_const_v<ElementType>) {
        return PDI_OUT;
    } else {
        return PDI_INOUT;
    }
};

template <class SupportType, class ElementType, class LayoutStridedPolicy>
void expose_to_pdi(
        std::string const& name,
        BlockSpan<SupportType, ElementType, LayoutStridedPolicy> const& data,
        PDI_inout_t access = default_access<ElementType>())
{
    auto const& extents = data.domain().extents().array();
    size_t rank = extents.size();
    if constexpr (std::is_const_v<ElementType>) {
        assert(!(access & PDI_IN) && "Can not use PDI to input constant data");
    }
    PDI_multi_expose(
            name.c_str(),
            (name + "_rank").c_str(),
            &rank,
            PDI_OUT,
            (name + "_extents").c_str(),
            const_cast<size_t*>(extents.data()),
            PDI_OUT,
            name.c_str(),
            data.data(),
            access,
            NULL);
}
