

# File assert.hpp

[**File List**](files.md) **>** [**src**](dir_68267d1309a1af8e8297ef4c3efbcdba.md) **>** [**utils**](dir_313caf1132e152dd9b58bea13a4052ca.md) **>** [**assert.hpp**](assert_8hpp.md)

[Go to the documentation of this file](assert_8hpp.md)


```C++
// SPDX-License-Identifier: MIT

#pragma once

#include "preprocessor.hpp"


#include <cstdlib>
// NOTE: We could also use __debugbreak or int3
#define GSLX_DEBUG_BREAK() std::exit(EXIT_FAILURE)



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
```


