// SPDX-License-Identifier: MIT

#pragma once

/**
 * @class MpiScopeGuard
 * @brief RAII wrapper for MPI initialisation and finalisation.
 *
 * This class manages the lifetime of the MPI environment using the
 * Resource Acquisition Is Initialisation (RAII) idiom. Constructing an
 * instance initialises MPI, and destroying the instance finalises MPI.
 *
 * The class is non-copyable and non-movable to ensure unique ownership
 * of the MPI runtime state within a process.
 */
class MpiScopeGuard
{
public:
    /**
     * @brief Initialises the MPI environment.
     *
     * Calls `MPI_Init` with no command-line arguments.
     *
     * @note This constructor assumes that MPI has not already been
     * initialised elsewhere.
     */
    MpiScopeGuard() noexcept;

    /**
     * @brief Initialises the MPI environment using command-line arguments.
     *
     * Calls `MPI_Init`, forwarding the provided command-line arguments
     * to the MPI implementation.
     *
     * @param argc Number of command-line arguments.
     * @param argv Array of command-line argument strings.
     *
     * @note The MPI implementation may modify `argc` and `argv`.
     * @note This constructor assumes that MPI has not already been
     * initialised elsewhere.
     */
    MpiScopeGuard(int& argc, char**& argv) noexcept;

    /**
     * @brief Deleted copy constructor.
     *
     * Copying an `MpiScopeGuard` is disallowed because MPI lifetime
     * ownership must remain unique.
     */
    MpiScopeGuard(MpiScopeGuard const& rhs) = delete;

    /**
     * @brief Deleted move constructor.
     *
     * Moving an `MpiScopeGuard` is disallowed because MPI lifetime
     * ownership cannot be transferred safely.
     */
    MpiScopeGuard(MpiScopeGuard&& rhs) noexcept = delete;

    /**
     * @brief Finalises the MPI environment.
     *
     * Calls `MPI_Finalize` if MPI was successfully initialised by
     * this guard instance.
     */
    ~MpiScopeGuard() noexcept;

    /**
     * @brief Deleted copy assignment operator.
     */
    MpiScopeGuard& operator=(MpiScopeGuard const& rhs) = delete;

    /**
     * @brief Deleted move assignment operator.
     */
    MpiScopeGuard& operator=(MpiScopeGuard&& rhs) noexcept = delete;
};
