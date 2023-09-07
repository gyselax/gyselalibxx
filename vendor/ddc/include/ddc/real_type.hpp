// SPDX-License-Identifier: MIT

#pragma once

namespace ddc {

#ifdef DDC_ENABLE_DOUBLE

using Real = double;

#else

using Real = float;

#endif

} // namespace ddc
