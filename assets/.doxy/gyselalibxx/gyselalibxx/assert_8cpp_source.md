

# File assert.cpp

[**File List**](files.md) **>** [**src**](dir_68267d1309a1af8e8297ef4c3efbcdba.md) **>** [**utils**](dir_313caf1132e152dd9b58bea13a4052ca.md) **>** [**assert.cpp**](assert_8cpp.md)

[Go to the documentation of this file](assert_8cpp.md)


```C++
#include <cstdio>

#include "assert.hpp"

namespace gslx {
namespace error {

void PresentErrorExplanation(const char* a_category, const char* an_explanation) noexcept
{
    std::fprintf(stderr, "[GSLX][%s] %s\n", a_category, an_explanation);
}

void AssertionFailure(const char* the_message) noexcept
{
    PresentErrorExplanation("ASSERTION_FAILURE", the_message);

    std::fflush(stdout);
    std::fflush(stderr);

    GSLX_DEBUG_BREAK();
}

} // namespace error
} // namespace gslx
```


