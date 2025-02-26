// SPDX-License-Identifier: MIT
#pragma once

#include <ginkgo/extensions/kokkos.hpp>
#include <ginkgo/ginkgo.hpp>

#include <Kokkos_Core.hpp>

#include "matrix_batch.hpp"
#include "matrix_utils.hpp"

/**
 * @brief  Matrix class which is able to manage and solve a batch of sparse linear systems. Executes on either CPU or GPU.
 * It takes advantage of the sparse structure, and the only batched solver available in Ginkgo : Stabilised Bicg.
 * The sparsity pattern is assumed to be the same for all matrices. ie the non-zero components are located at the same places for all matrices.
 * This class uses the ELL storage format which needs two 1D arrays, one stores values the other column indices.
 * The class returns these arrays (as Kokkos views) with the get_batch_idx_and_vals function, it is then possible to fill them outside the class.
 * Tolerance and maximal number of iterations, which are parameters for the iterative solver, are set in the constructor.
 * It is possible to get convergence information by activating the logger at constructor call.
 * @tparam ExecSpace Execution space,needed by Kokkos for allocations and parallelism.
 * The simplest choice is to follow Kokkos, for that: specify Kokkos::DefaultExecutionSpace
 */
template <class ExecSpace>
class MatrixBatchEll : public MatrixBatch<ExecSpace>
{
public:
    using typename MatrixBatch<ExecSpace>::BatchedRHS;
    using MatrixBatch<ExecSpace>::size;
    using MatrixBatch<ExecSpace>::batch_size;

private:
    using batch_sparse_type = gko::batch::matrix::Ell<double, int>;
    using solver_type = gko::batch::solver::Bicgstab<double>;

    std::shared_ptr<batch_sparse_type> m_batch_matrix_ell;
    std::shared_ptr<solver_type> m_solver;
    int m_max_iter;
    double m_tol;
    bool m_with_logger;


public:
    /**
     * @brief The constructor for MatrixBatchEll class.
     *
     * @param[in] batch_size Number of linear systems to solve.
     * @param[in] mat_size Common matrix size for all the systems.
     * @param[in] non_zeros_per_row number of non zero components per line.
     * @param[in] max_iter maximal number of iterations for the solver
     * @param[in] res_tol residual tolerance parameter, to ensure convergence. Be careful! the relative residual 
     * provided here, will be used as "implicit residual" in ginkgo solver.
     * @param[in] logger boolean parameter for saving log information such residual and interations count.
    */
    explicit MatrixBatchEll(
            const int batch_size,
            const int mat_size,
            const int non_zeros_per_row,
            std::optional<int> max_iter = std::nullopt,
            std::optional<double> res_tol = std::nullopt,
            std::optional<bool> logger = std::nullopt)
        : MatrixBatch<ExecSpace>(batch_size, mat_size)
        , m_max_iter(max_iter.value_or(500))
        , m_tol(res_tol.value_or(1e-15))
        , m_with_logger(logger.value_or(false))
    {
        std::shared_ptr const gko_exec = gko::ext::kokkos::create_executor(ExecSpace());
        m_batch_matrix_ell = gko::share(
                batch_sparse_type::
                        create(gko_exec,
                               gko::batch_dim<2>(batch_size, gko::dim<2>(mat_size, mat_size)),
                               non_zeros_per_row));
    }

    /**
     * @brief Constructor for MatrixBatchEll class.
     *
     * @param[in] cols_idx  A Kokkos view which stores the column indices of non-zero components.
     * @param[in] batch_values A Kokkos view which stores the values of non-zero elements.
     * @param[in] max_iter maximal number of iterations for the solver, default 500.
     * @param[in] res_tol residual tolerance parameter, to ensure convergence. Be careful! The residual 
     * provided here, set as relative residual, will be used as "implicit residual" in ginkgo solver.
     * Default value is set to 1e-15.
     * @param[in] logger boolean parameter to save logger information. Default value false.
     */
    explicit MatrixBatchEll(
            Kokkos::View<int**, Kokkos::LayoutLeft, ExecSpace> cols_idx,
            Kokkos::View<double***, Kokkos::LayoutStride, ExecSpace> batch_values,
            std::optional<int> max_iter = std::nullopt,
            std::optional<double> res_tol = std::nullopt,
            std::optional<bool> logger = std::nullopt)
        : MatrixBatch<ExecSpace>(batch_values.extent(0), batch_values.extent(1))
        , m_max_iter(max_iter.value_or(500))
        , m_tol(res_tol.value_or(1e-15))
        , m_with_logger(logger.value_or(false))


    {
        std::shared_ptr const gko_exec = gko::ext::kokkos::create_executor(ExecSpace());
        m_batch_matrix_ell = gko::share(
                gko::batch::matrix::Ell<double>::
                        create(gko_exec,
                               gko::batch_dim<2>(batch_size(), gko::dim<2>(size(), size())),
                               batch_values.extent(2),
                               gko::array<double>::
                                       view(gko_exec, batch_values.span(), batch_values.data()),
                               gko::array<int>::view(gko_exec, cols_idx.span(), cols_idx.data())));
    }

    /**
    * @brief A function to get information about values and indices for the whole batch.
    * Data is managed by two Kokkos Views stored on the host.
    * @return idx_view   Column indices for the non-zero values.
    * @return vals_view  The non-zero values.
    */
    std::pair<
            Kokkos::View<int**, Kokkos::LayoutLeft, ExecSpace>,
            Kokkos::View<double***, Kokkos::LayoutStride, ExecSpace>>
    get_batch_ell()
    {
        int* idx_buffer = m_batch_matrix_ell->get_col_idxs();
        double* vals_buffer = m_batch_matrix_ell->get_values();
        Kokkos::LayoutStride values_layout(
                batch_size(),
                m_batch_matrix_ell->get_num_stored_elements_per_row() * size(),
                size(),
                1,
                m_batch_matrix_ell->get_num_stored_elements_per_row(),
                size());
        Kokkos::View<int**, Kokkos::LayoutLeft, ExecSpace>
                idx_view(idx_buffer, size(), m_batch_matrix_ell->get_num_stored_elements_per_row());
        Kokkos::View<double***, Kokkos::LayoutStride, ExecSpace>
                vals_view(vals_buffer, values_layout);
        return {idx_view, vals_view};
    }

    /**
    * @brief A getter function for a value located at a specified place.
    *  @param[in] batch_idx Index in the batch.
    *  @param[in] line_idx Line index inside the matrix.
    *  @param[in] non_zero_col_idx Non-zero index element in the line.
    * @return  value of the component.
    */
    double get_ell_element(int batch_idx, int line_idx, int non_zero_col_idx) const
    {
        // Checks index of column according to the sparsity pattern.
        assert(non_zero_col_idx >= 0
               && non_zero_col_idx <= m_batch_matrix_ell->get_num_stored_elements_per_row());
        return m_batch_matrix_ell->get_const_values()
                [batch_idx * size() * m_batch_matrix_ell->get_num_stored_elements_per_row()
                 + non_zero_col_idx * size() + line_idx];
    }

    /**
    * @brief A setter function to modify a value located at a specified place.
    *  @param[in] batch_idx Index in the batch.
    *  @param[in] line_idx Line index inside the matrix.
    *  @param[in] non_zero_col_idx Non-zero index element in the line.
    *  @param[in] aij New value.
    */
    void set_ell_element(int batch_idx, int line_idx, int non_zero_col_idx, double aij)
    {
        // Checks index of column according to the sparsity pattern.
        assert(non_zero_col_idx >= 0
               && non_zero_col_idx <= m_batch_matrix_ell->get_num_stored_elements_per_row());
        m_batch_matrix_ell->get_values()
                [batch_idx * size() * m_batch_matrix_ell->get_num_stored_elements_per_row()
                 + non_zero_col_idx * size() + line_idx]
                = aij;
    }

    /**
     * @brief Perform a pre-process operation on the solver. Must be called after filling the matrix.
     *
     * It uses parameters like maximum number of iterations and tolerance are used to instantiate a Ginkgo solver.
     *
     * The stopping criterion is a reduction factor ||Axe-b||/||b||<tol with max_iter maximum iterations.
     */
    void setup_solver() final
    {
        std::shared_ptr const gko_exec = m_batch_matrix_ell->get_executor();
        gko::batch::stop::tolerance_type tol_type = gko::batch::stop::tolerance_type::relative;

        std::shared_ptr solver_factory = solver_type::build()
                                                 .with_max_iterations(m_max_iter)
                                                 .with_tolerance(m_tol)
                                                 .with_tolerance_type(tol_type)
                                                 .on(gko_exec);
        m_solver = solver_factory->generate(m_batch_matrix_ell);
        gko_exec->synchronise();
    }

    /**
     * @brief Solve the batched linear problem Axe=b.
     *
     * @param[in, out] b A 2D Kokkos::View storing the batched right-hand sides of the problem and receiving the corresponding solutions.
     */
    void solve(BatchedRHS const b) const final
    {
        std::shared_ptr const gko_exec = m_solver->get_executor();
        BatchedRHS x_view("x_view", batch_size(), size());

        // Create a logger to obtain the iteration counts and "implicit" residual norms for every system after the solve.
        std::shared_ptr<const gko::batch::log::BatchConvergence<double>> logger
                = gko::batch::log::BatchConvergence<double>::create();
        m_solver->add_logger(logger);
        gko_exec->synchronise();

        Kokkos::deep_copy(x_view, b);
        m_solver->apply(to_gko_multivector(gko_exec, b), to_gko_multivector(gko_exec, x_view));
        m_solver->remove_logger(logger);
        // save logger data
        if (m_with_logger) {
            std::fstream log_file("ell_log.txt", std::ios::out | std::ios::app);
            save_logger(log_file, m_batch_matrix_ell, x_view, b, logger, m_tol);
            log_file.close();
        }
        //check convergence
        check_conv(batch_size(), m_tol, gko_exec, logger);

        Kokkos::deep_copy(b, x_view);
    }

    /**
    * @brief A function returns the norm of a matrix located at batch_idx.
    * @param[in] batch_idx integer, index of the matrix in the batch. 
    * @return  value of the matrix infinite-norm.
    */
    double norm(int batch_idx) const
    {
        int const tmp_mat_size = size();
        int const tmp_batch_size = batch_size();
        int const non_zeros = m_batch_matrix_ell->get_num_stored_elements_per_row();
        double* vals_proxy = m_batch_matrix_ell->get_values();
        Kokkos::LayoutStride values_layout(
                tmp_batch_size,
                non_zeros * tmp_mat_size,
                tmp_mat_size,
                1,
                non_zeros,
                tmp_mat_size);
        Kokkos::View<double***, Kokkos::LayoutStride, typename ExecSpace::memory_space>
                vals_view(vals_proxy, values_layout);

        double result = 0;
        Kokkos::parallel_reduce(
                "L-infinitty norm",
                Kokkos::RangePolicy<ExecSpace>(0, tmp_mat_size),
                KOKKOS_LAMBDA(int i, double& res) {
                    double row_sum = 0.;
                    for (int k = 0; k < non_zeros; k++) {
                        row_sum += Kokkos::abs(vals_view(batch_idx, i, k));
                    }
                    if (row_sum > res) {
                        res = row_sum;
                    }
                },
                Kokkos::Max<double>(result));
        return result;
    }
};
