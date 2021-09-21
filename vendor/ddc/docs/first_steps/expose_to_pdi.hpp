#pragma once

#include <string>

#include <ddc/Block>
#include <ddc/pdi.hpp>

#include <pdi.h>

constexpr char const* const PDI_CFG = R"PDI_CFG(
metadata:
  iter : int

data:
  temperature_extents: { type: array, subtype: int64, size: 2 }
  temperature_subextents: { type: array, subtype: int64, size: 2 }
  temperature:
    type: array
    subtype: double
    size: [ '$temperature_extents[0]', '$temperature_extents[1]' ]
    subsize: [ '$temperature_subextents[0]', '$temperature_subextents[1]' ]

plugins:
  decl_hdf5:
    - file: 'temperature_${iter:04}.h5'
      on_event: temperature
      collision_policy: replace_and_warn
      write: [temperature]
  trace: ~
)PDI_CFG";

template <class SupportType, class ElementType>
void my_expose_to_pdi(
        std::string const& name,
        Block<ElementType, SupportType> const& data,
        SupportType const& subdomain,
        PDI_inout_t access = default_access_v<ElementType>)
{
    auto extents = data.extents().array();
    auto subextents = subdomain.extents().array();
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
            extents.data(),
            PDI_OUT,
            (name + "_subextents").c_str(),
            subextents.data(),
            PDI_OUT,
            name.c_str(),
            data.data(),
            access,
            NULL);
}
