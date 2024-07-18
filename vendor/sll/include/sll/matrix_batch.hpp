#pragma once

#include <Kokkos_Core.hpp>

#include "sll/view.hpp"


/**
 * @brief MatrixBatch superclass for managing a collection of linear systems. 
 * The main assumption is that all matrices have the same size.
 * It is also assumed that each matrix is used to solve one equation.
 * 
 * Classes inheriting from this class must manage other aspects:
 * Sparsity: kind of storage (Dense, Ell, Csr, etc.)
 * Kind of solver (direct, iterative)
 * Preconditioners and factorizations
 *
 * @tparam ExecSpace Execution space,needed by Kokkos for allocations and parallelism.
 * The simplest choice is to follow Kokkos, for that: specify Kokkos::DefaultExecutionSpace
 */

template <typename ExecSpace>
class MatrixBatch
{
public:
    /**
     * @brief Alias for 2D double Kokkos views, LayoutRight is specified.
     */
    using DKokkosView2D
            = Kokkos::View<double**, Kokkos::LayoutRight, typename ExecSpace::memory_space>;

    /**
     * @brief The constructor for MatrixBatch class.
     * @param[in] batch_size Number of linear systems to solve.
     * @param[in] mat_size Common matrix size for all the systems.
     */
    MatrixBatch(int batch_size, int mat_size) : m_batch_size(batch_size), m_matrix_size(mat_size) {}

    virtual ~MatrixBatch() = default;
    /**
     * @brief Generic function which could be used by children classes to execute
     * class specific implementation details.
     */
    virtual void factorize() = 0;

    /**
     * @brief A function which solves the collection of linear problems.
     * @param[inout] b 2d Kokkos view which stores the right hand side, 
     * @return  The computation result, stored in b.
     */
    virtual DKokkosView2D solve_inplace(DKokkosView2D b) const = 0;

    /**
     * @brief Getter function for matrix size.
     * @return  Value of common matrix size.
     */
    int get_size() const
    {
        return m_matrix_size;
    }

    /**
     * @brief Getter function for batch size, ie the number of linear systems to solve.
     * @return Value of batch size.
     */
    int get_batch_size() const
    {
        return m_batch_size;
    }

private:
    int const m_batch_size;
    int const m_matrix_size;
};
