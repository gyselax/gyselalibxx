// SPDX-License-Identifier: MIT

#pragma once

#include "preprocessor.hpp"

////////////////////////////////////////
/// GSLX_DEBUG_BREAK
///
/// This macro will stop the program.
///
/// Example usage:
///
///     GSLX_DEBUG_BREAK();
///
////////////////////////////////////////

#include <cstdlib>
// NOTE: We could also use __debugbreak or int3
#define GSLX_DEBUG_BREAK() std::exit(EXIT_FAILURE)


////////////////////////////////////////
/// GSLX_ASSERT
///
/// Use this macro to assert a condition. The assertions trigger only if
/// GSLX_ASSERT_ENABLED is defined. Use that only for cheap assertion. The cost
/// should be negligible compared to the computation.
///
/// Example usage:
///
///     GSLX_ASSERT(x > 10);
///
////////////////////////////////////////

namespace gslx {
namespace error {
void AssertionFailure(const char* the_message) noexcept;
}
} // namespace gslx

// NOTE: For now we always enable this kind of assertion. We could use an option
// and a CMake configure_file.
#define GSLX_ASSERT_ENABLED 1

#if defined(GSLX_ASSERT_ENABLED)
#define GSLX_ASSERT(an_expression)                                                                 \
    do {                                                                                           \
        GSLX_UNUSED(                                                                               \
                (an_expression)                                                                    \
                || (gslx::error::AssertionFailure(                                                 \
                            #an_expression " in " __FILE__ ":" GSLX_UTILITY_STRINGIFY(__LINE__)),  \
                    false));                                                                       \
    } while (false)
#else
#define GSLX_ASSERT(an_expression)
#endif
