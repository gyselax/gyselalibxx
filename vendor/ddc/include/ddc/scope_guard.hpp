// Copyright (C) The DDC development team, see COPYRIGHT.md file
//
// SPDX-License-Identifier: MIT

#pragma once

#include <functional>
#include <map>
#include <optional>
#include <string>

#include "discrete_space.hpp"

namespace ddc {

class ScopeGuard
{
    void discretization_store_initialization() const
    {
        detail::g_discretization_store
                = std::make_optional<std::map<std::string, std::function<void()>>>();
    }

public:
    ScopeGuard()
    {
        discretization_store_initialization();
    }

    ScopeGuard([[maybe_unused]] int argc, [[maybe_unused]] char**& argv)
    {
        discretization_store_initialization();
    }

    ScopeGuard(ScopeGuard const& x) = delete;

    ScopeGuard(ScopeGuard&& x) noexcept = delete;

    ~ScopeGuard() noexcept
    {
        for (auto const& [name, fn] : *detail::g_discretization_store) {
            fn();
        }
        detail::g_discretization_store.reset();
    }

    ScopeGuard& operator=(ScopeGuard const& x) = delete;

    ScopeGuard& operator=(ScopeGuard&& x) noexcept = delete;
};

} // namespace ddc
