// SPDX-License-Identifier: MIT

#pragma once

#include <ddc/ddc.hpp>

template <class... DDims>
using NDTag = ddc::detail::TypeSeq<DDims...>;
