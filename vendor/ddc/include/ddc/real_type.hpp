// Copyright (C) The DDC development team, see COPYRIGHT.md file
//
// SPDX-License-Identifier: MIT

#pragma once

namespace ddc {

#ifdef DDC_BUILD_DOUBLE_PRECISION

using Real = double;

#else

using Real = float;

#endif

} // namespace ddc
