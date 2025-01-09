// SPDX-License-Identifier: MIT
#pragma once

#include <sll/view.hpp>

#include <Kokkos_Core.hpp>


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
    /// @brief The type of a Kokkos::View storing batched right-hand sides. Second dimenion is batch dimension.
    using BatchedRHS = Kokkos::View<double**, Kokkos::LayoutRight, ExecSpace>;


private:
    std::size_t m_size;
    std::size_t m_batch_size;

protected:
    /**
     * @brief The constructor for MatrixBatch class.
     * @param[in] batch_size Number of linear systems to solve.
     * @param[in] mat_size Common matrix size for all the systems.
     */
    explicit MatrixBatch(const std::size_t batch_size, const std::size_t mat_size)
        : m_size(mat_size)
        , m_batch_size(batch_size)
    {
    }

public:
    /// @brief Destruct
    virtual ~MatrixBatch() = default;

    /**
     * @brief Perform a pre-process operation on the solver. Must be called after filling the matrix.
     */
    virtual void setup_solver() = 0;

    /**
     * @brief Solve the multiple right-hand sides linear problem Ax=b inplace.
     *
     * @param[in, out] b A 2D Kokkos::View storing the batched right-hand sides of the problem and receiving the corresponding solution.
     */
    virtual void solve(BatchedRHS b) const = 0;

    /**
     * @brief Get the size of the square matrix corresponding to a single batch in one of its dimensions.
     *
     * @return The size of the matrix in one of its dimensions.
     */
    std::size_t size() const
    {
        return m_size;
    }

    /**
     * @brief Get the batch size of the linear problem.
     *
     * @return The batch size of the linear problem.
     */
    std::size_t batch_size() const
    {
        return m_batch_size;
    }
};
