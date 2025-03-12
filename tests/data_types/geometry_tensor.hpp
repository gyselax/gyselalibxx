// SPDX-License-Identifier: MIT
#pragma once

namespace {

struct R_cov;
struct Theta_cov;

struct R
{
    static bool constexpr IS_COVARIANT = false;
    static bool constexpr IS_CONTRAVARIANT = true;
    using Dual = R_cov;
};
struct Theta
{
    static bool constexpr IS_COVARIANT = false;
    static bool constexpr IS_CONTRAVARIANT = true;
    using Dual = Theta_cov;
};

struct R_cov
{
    static bool constexpr IS_COVARIANT = true;
    static bool constexpr IS_CONTRAVARIANT = false;
    using Dual = R;
};
struct Theta_cov
{
    static bool constexpr IS_COVARIANT = true;
    static bool constexpr IS_CONTRAVARIANT = false;
    using Dual = Theta;
};

struct X
{
    static bool constexpr IS_COVARIANT = true;
    static bool constexpr IS_CONTRAVARIANT = true;
    using Dual = X;
};
struct Y
{
    static bool constexpr IS_COVARIANT = true;
    static bool constexpr IS_CONTRAVARIANT = true;
    using Dual = Y;
};

} // namespace
