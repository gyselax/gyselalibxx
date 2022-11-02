#pragma once

#include <map>
#include <memory>
#include <optional>

#include <Kokkos_Core.hpp>

#include "discrete_space.hpp"

namespace ddc {

class ScopeGuard
{
    Kokkos::ScopeGuard m_kokkos_scope_guard;

    void discretization_store_initialization() const
    {
        ddc_detail::g_discretization_store
                = std::make_optional<std::map<std::string, std::function<void()>>>();
    }

public:
    ScopeGuard()
    {
        discretization_store_initialization();
    }

    ScopeGuard(int argc, char**& argv) : m_kokkos_scope_guard(argc, argv)
    {
        discretization_store_initialization();
    }

    ScopeGuard(ScopeGuard const& x) = delete;

    ScopeGuard(ScopeGuard&& x) noexcept = delete;

    ~ScopeGuard() noexcept
    {
        for (auto const& [name, fn] : *ddc_detail::g_discretization_store) {
            fn();
        }
        ddc_detail::g_discretization_store.reset();
    }

    ScopeGuard& operator=(ScopeGuard const& x) = delete;

    ScopeGuard& operator=(ScopeGuard&& x) noexcept = delete;
};

} // namespace ddc
