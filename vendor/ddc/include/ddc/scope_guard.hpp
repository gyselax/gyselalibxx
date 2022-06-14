#pragma once

#include <Kokkos_Core.hpp>

class ScopeGuard
{
    Kokkos::ScopeGuard m_kokkos_scope_guard;

public:
    ScopeGuard() = default;

    ScopeGuard(int argc, char**& argv) : m_kokkos_scope_guard(argc, argv) {}

    ScopeGuard(ScopeGuard const& x) = delete;

    ScopeGuard(ScopeGuard&& x) noexcept = delete;

    ~ScopeGuard() noexcept = default;

    ScopeGuard& operator=(ScopeGuard const& x) = delete;

    ScopeGuard& operator=(ScopeGuard&& x) noexcept = delete;
};
