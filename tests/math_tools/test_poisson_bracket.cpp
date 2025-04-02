#include "../mapping/geometry_mapping_tests.hpp"

namespace {

struct Phi_cov;

struct Phi
{
    static bool constexpr PERIODIC = true;
    static bool constexpr IS_COVARIANT = false;
    static bool constexpr IS_CONTRAVARIANT = true;
    /// The corresponding type in the dual space.
    using Dual = Phi_cov;
};

struct Phi_cov
{
    static bool constexpr PERIODIC = true;
    static bool constexpr IS_COVARIANT = true;
    static bool constexpr IS_CONTRAVARIANT = false;
    /// The corresponding type in the dual space.
    using Dual = Phi;
};

} // namespace


TEST(GyrokineticPoissonBracket, Anticommutativity) {}

TEST(GyrokineticPoissonBracket, Bilinearity) {}

TEST(GyrokineticPoissonBracket, Leibniz) {}

TEST(GyrokineticPoissonBracket, Jacobi) {}
