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
