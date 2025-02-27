// SPDX-License-Identifier: MIT
#pragma once

#include <ginkgo/extensions/kokkos.hpp>
#include <ginkgo/ginkgo.hpp>

#include <Kokkos_Core.hpp>

#include "matrix_batch.hpp"
#include "matrix_utils.hpp"

/**
* @brief A tag to choose between the batched iterative solvers provided by Ginkgo.
*
* BICGSTAB BiConjugate Gradient Stabilise method. BICGSTAB is able to solve general sparse matrices problems.
* CG Conjugate Gradient method (For symmetric positive definite matrices).
* If the matrices structures verify CG requirements, we encourage the use of this method for its lower computation cost.
*/
enum class MatrixBatchCsrSolver { CG, BICGSTAB, BATCH_CG, BATCH_BICGSTAB };

/**
 * @brief  Matrix class which is able to manage and solve a batch of sparse linear systems. Executes on either CPU or GPU.
 * It takes advantage of the sparse structure, and the only batched solver available in Ginkgo : Stabilised Bicg.
 * This class uses the CSR storage format which needs three arrays, one stores values, the other column indices.
 * The third array contains the count of non-zero inside the matrix lines.(eg:for a given line index i nn_per_row[i]= sum of non-zeros until line i)
 * The class returns these arrays (as Kokkos views) with the get_csr_views function, it is then possible to fill them outside the class.
 * The sparsity pattern is the same for all matrices, hence column indices are stored only for one system.
 * Tolerance and maximal number of iterations, which are parameters for the iterative solver, are set in the constructor.
 * It is possible to get convergence information by activating the logger at constructor call.
 * @tparam ExecSpace Execution space,needed by Kokkos for allocations and parallelism.
 * The simplest choice is to follow Kokkos, for that: specify Kokkos::DefaultExecutionSpace
 * @tparam Solver Refers to the solver type, default value is the Bicgstab which is more general.
 * The use of a CG solver is also possible, in this case, please make sure that matrices structure fulfils CG requirements.
 */
template <class ExecSpace, MatrixBatchCsrSolver Solver = MatrixBatchCsrSolver::BICGSTAB>
class MatrixBatchCsr : public MatrixBatch<ExecSpace>
{
public:
    using typename MatrixBatch<ExecSpace>::BatchedRHS;
    using MatrixBatch<ExecSpace>::size;
    using MatrixBatch<ExecSpace>::batch_size;

private:
    using batch_sparse_type = gko::batch::matrix::Csr<double, int>;
    using solver_type = std::conditional_t<
            Solver == MatrixBatchCsrSolver::CG,
            gko::solver::Cg<double>,
            std::conditional_t<
                    Solver == MatrixBatchCsrSolver::BICGSTAB,
                    gko::solver::Bicgstab<double>,
                    std::conditional_t<
                            Solver == MatrixBatchCsrSolver::BATCH_CG,
                            gko::batch::solver::Cg<double>,
                            gko::batch::solver::Bicgstab<double>>>>;

    std::shared_ptr<batch_sparse_type> m_batch_matrix_csr;
    std::conditional_t<
            Solver == MatrixBatchCsrSolver::CG || Solver == MatrixBatchCsrSolver::BICGSTAB,
            std::vector<std::shared_ptr<solver_type>>,
            std::shared_ptr<solver_type>>
            m_solver;
    int m_max_iter;
    double m_tol;
    bool m_with_logger;
    unsigned int m_preconditionner_max_block_size; // Maximum size of Jacobi-block preconditionner

public:
    /**
     * @brief The constructor for MatrixBatchCsr class.
     *
     * @param[in] batch_size Number of linear systems to solve.
     * @param[in] mat_size Common matrix size for all the systems.
     * @param[in] nnz_per_system Number of non-zero components per matrix.
     * @param[in] max_iter maximal number of iterations for the solver, default 1000.
     * @param[in] res_tol residual tolerance parameter, to ensure convergence. Be careful! the relative residual 
     * provided here, will be used as "implicit residual" in ginkgo solver.
     * @param[in] logger boolean parameter for saving log information such residual and interactions count.
     * @param[in] preconditionner_max_block_size An optional parameter used to define the maximum size of a block
     */
    explicit MatrixBatchCsr(
            const int batch_size,
            const int mat_size,
            const int nnz_per_system,
            std::optional<int> max_iter = std::nullopt,
            std::optional<double> res_tol = std::nullopt,
            std::optional<bool> logger = std::nullopt,
            std::optional<int> preconditionner_max_block_size = 1u)
        : MatrixBatch<ExecSpace>(batch_size, mat_size)
        , m_max_iter(max_iter.value_or(1000))
        , m_tol(res_tol.value_or(1e-15))
        , m_with_logger(logger.value_or(false))
        , m_preconditionner_max_block_size(preconditionner_max_block_size.value_or(
                  default_preconditionner_max_block_size<ExecSpace>()))
    {
        std::shared_ptr const gko_exec = gko::ext::kokkos::create_executor(ExecSpace());
        m_batch_matrix_csr = gko::share(
                batch_sparse_type::
                        create(gko_exec,
                               gko::batch_dim<2>(batch_size, gko::dim<2>(mat_size, mat_size)),
                               nnz_per_system));
    }

    /**
     * @brief Constructor for MatrixBatchCsr class.
     *
     * @param[in] batch_values A 2D Kokkos view which stores the values of non-zero elements for the whole batch.
     * @param[in] cols_idx  A 1D Kokkos view which stores the column indices for each non-zero component.(only for one matrix)
     * @param[in] nnz_per_row  A 1D Kokkos view of length matrix_size+1 which stores the count of the non-zeros along the lines of the matrix.
     * It is defined as:
     * nnz_per_row[0] = 0.
     * nnz_per_row[matrix_size] = total_number_of_nonzero.
     * To get the number of non-zero for a line i,one have to compute : n_non_zeros_at_line_in = nnz_per_row[i+1]-nnz_per_row[i].
     * @param[in] max_iter maximal number of iterations for the solver, default 1000.
     * @param[in] res_tol residual tolerance parameter, to ensure convergence. Be careful! The residual 
     * provided here, set as relative residual, will be used as "implicit residual" in ginkgo solver.
     * Default value is set to 1e-15.
     * @param[in] logger boolean parameter to save logger information. Default value false.
     * @param[in] preconditionner_max_block_size An optional parameter used to define the maximum size of a block
     */
    explicit MatrixBatchCsr(
            Kokkos::View<double**, Kokkos::LayoutRight, ExecSpace> batch_values,
            Kokkos::View<int*, Kokkos::LayoutRight, ExecSpace> cols_idx,
            Kokkos::View<int*, Kokkos::LayoutRight, ExecSpace> nnz_per_row,
            std::optional<int> max_iter = std::nullopt,
            std::optional<double> res_tol = std::nullopt,
            std::optional<bool> logger = std::nullopt,
            std::optional<int> preconditionner_max_block_size = 1u)
        : MatrixBatch<ExecSpace>(batch_values.extent(0), nnz_per_row.size() - 1)
        , m_max_iter(max_iter.value_or(1000))
        , m_tol(res_tol.value_or(1e-15))
        , m_with_logger(logger.value_or(false))
        , m_preconditionner_max_block_size(preconditionner_max_block_size.value_or(
                  default_preconditionner_max_block_size<ExecSpace>()))
    {
        std::shared_ptr const gko_exec = gko::ext::kokkos::create_executor(ExecSpace());
        m_batch_matrix_csr = gko::share(
                batch_sparse_type::
                        create(gko_exec,
                               gko::batch_dim<2>(batch_size(), gko::dim<2>(size(), size())),
                               std::move(gko::array<double>::
                                                 view(gko_exec,
                                                      batch_values.span(),
                                                      batch_values.data())),
                               std::move(gko::array<
                                         int>::view(gko_exec, cols_idx.span(), cols_idx.data())),
                               std::move(gko::array<int>::
                                                 view(gko_exec,
                                                      nnz_per_row.span(),
                                                      nnz_per_row.data()))));
    }

    /**
     * @brief A function to update information about values,indices and the number of non-zero per row for the whole batch.
     *        Data is managed by Kokkos Views stored on the host.
     * @return vals_view A 2D Kokkos view which stores the values of non-zero elements for the whole batch.
     * @return col_idx_view  A 1D Kokkos view which stores the column indices for each non-zero component.(only for one matrix)
     * @return nnz_per_row_view A 1D Kokkos view which stores the count of non-zero per line, in an additive way.
     *         see nnz_per_row parameter in constructor.
     */
    std::tuple<
            Kokkos::View<double**, Kokkos::LayoutRight, ExecSpace>,
            Kokkos::View<int*, Kokkos::LayoutRight, ExecSpace>,
            Kokkos::View<int*, Kokkos::LayoutRight, ExecSpace>>
    get_batch_csr()
    {
        int* idx_buffer = m_batch_matrix_csr->get_col_idxs();
        int* nnz_per_row_buffer = m_batch_matrix_csr->get_row_ptrs();
        double* vals_buffer = m_batch_matrix_csr->get_values();

        Kokkos::View<int*, Kokkos::LayoutRight, ExecSpace>
                col_idx_view_buffer(idx_buffer, m_batch_matrix_csr->get_num_elements_per_item());
        Kokkos::View<int*, Kokkos::LayoutRight, ExecSpace>
                nnz_per_row_view_buffer(nnz_per_row_buffer, size() + 1);
        Kokkos::View<double**, Kokkos::LayoutRight, ExecSpace> vals_view_buffer(
                vals_buffer,
                batch_size(),
                m_batch_matrix_csr->get_num_elements_per_item());

        return {vals_view_buffer, col_idx_view_buffer, nnz_per_row_view_buffer};
    }

    /**
     * @brief Perform a pre-process operation on the solver. Must be called after filling the matrix.
     *
     * It uses the batch of matrices to generate a batched Jacobi preconditioner. Other parameters like maximum number
     * of iterations and tolerance are also used to instantiate a Ginkgo solver.
     *
     * The stopping criterion is a reduction factor ||Axe-b||/||b||<tol with max_iter maximum iterations.
     */
    void setup_solver() final
    {
        std::shared_ptr const gko_exec = m_batch_matrix_csr->get_executor();
        //Check if the indices array is sorted, and sort it if necessary.
        //The values array corresponding to the indices is also reordered.
        for (size_t i = 0; i < batch_size(); i++) {
            std::unique_ptr<gko::matrix::Csr<double, int>> tmp_matrix
                    = m_batch_matrix_csr->create_view_for_item(i);
            if (!tmp_matrix->is_sorted_by_column_index()) {
                tmp_matrix->sort_by_column_index();
            }
        }
        if constexpr (
                Solver == MatrixBatchCsrSolver::CG || Solver == MatrixBatchCsrSolver::BICGSTAB) {
            // Create the solver factory
            std::shared_ptr const residual_criterion
                    = gko::stop::ResidualNorm<double>::build().with_reduction_factor(m_tol).on(
                            gko_exec);

            std::shared_ptr const iterations_criterion
                    = gko::stop::Iteration::build().with_max_iters(m_max_iter).on(gko_exec);

            std::shared_ptr const preconditioner
                    = gko::preconditioner::Jacobi<double>::build()
                              .with_max_block_size(m_preconditionner_max_block_size)
                              .on(gko_exec);

            std::unique_ptr const solver_factory
                    = solver_type::build()
                              .with_preconditioner(preconditioner)
                              .with_criteria(residual_criterion, iterations_criterion)
                              .on(gko_exec);

            // Create the solvers
            for (size_t i = 0; i < batch_size(); i++) {
                m_solver.emplace_back(solver_factory->generate(
                        m_batch_matrix_csr->create_const_view_for_item(i)));
            }
        } else {
            // Create the solver factory
            std::shared_ptr const preconditioner
                    = gko::batch::preconditioner::Jacobi<double, int>::build()
                              .with_max_block_size(m_preconditionner_max_block_size)
                              .on(gko_exec);

            std::shared_ptr solver_factory = solver_type::build()
                                                     .with_max_iterations(m_max_iter)
                                                     .with_tolerance(m_tol)
                                                     .with_preconditioner(preconditioner)
                                                     .on(gko_exec);

            // Create the solver
            m_solver = solver_factory->generate(m_batch_matrix_csr);
        }
        gko_exec->synchronize();
    }

    /**
     * @brief Solve the batched linear problem Axe=b.
     *
     * @param[in, out] b A 2D Kokkos::View storing the batched right-hand sides of the problem and receiving the corresponding solutions.
     */
    void solve(BatchedRHS const b) const final
    {
        BatchedRHS x_view("x_view", batch_size(), size());
        Kokkos::deep_copy(x_view, b);
        solve(x_view, b);
        Kokkos::deep_copy(b, x_view);
    }

    /**
     * @brief Solve the batched linear problem Axe=b.
     *
     * @param[in] x A 2D Kokkos::View storing the batched  initial guests (useful for iterative solver) of the problems, and receiving the corresponding solutions.
     *
     * @param[in] b A 2D Kokkos::View storing the batched right-hand side of the problems.
     */
    void solve(BatchedRHS const x, BatchedRHS const b) const
    {
        if constexpr (
                Solver == MatrixBatchCsrSolver::CG || Solver == MatrixBatchCsrSolver::BICGSTAB) {
            for (size_t i = 0; i < batch_size(); i++) {
                std::shared_ptr const gko_exec = m_solver[i]->get_executor();

                // Create a logger to obtain the iteration counts and "implicit" residual norms for every system after the solve.
                std::shared_ptr const logger = gko::log::Convergence<double>::create();

                // Solve & log
                m_solver[i]->add_logger(logger);
                m_solver[i]
                        ->apply(to_gko_multivector(gko_exec, b)->create_const_view_for_item(i),
                                to_gko_multivector(gko_exec, x)->create_view_for_item(i));
                m_solver[i]->remove_logger(logger);
                // save logger data
                if (m_with_logger) {
                    std::fstream log_file("csr_log.txt", std::ios::out | std::ios::app);
                    save_logger(
                            log_file,
                            i,
                            m_batch_matrix_csr->create_const_view_for_item(i),
                            Kokkos::subview(x, i, Kokkos::ALL),
                            Kokkos::subview(b, i, Kokkos::ALL),
                            logger,
                            m_tol);
                    log_file.close();
                }
                // Check convergency
                if (!logger->has_converged()) {
                    throw std::runtime_error("Ginkgo did not converge in MatrixBatchCsr");
                }
            }
        } else {
            std::shared_ptr const gko_exec = m_solver->get_executor();

            // Create a logger to obtain the iteration counts and "implicit" residual norms for every system after the solve.
            std::shared_ptr const logger = gko::batch::log::BatchConvergence<double>::create();

            // Solve & log
            m_solver->add_logger(logger);
            m_solver->apply(to_gko_multivector(gko_exec, b), to_gko_multivector(gko_exec, x));
            m_solver->remove_logger(logger);

            // Save logger data
            if (m_with_logger) {
                std::fstream log_file("csr_log.txt", std::ios::out | std::ios::app);
                save_logger(log_file, m_batch_matrix_csr, x, b, logger, m_tol);
                log_file.close();
            }
            // Check convergence
            check_conv(batch_size(), m_tol, gko_exec, logger);
        }
    }

    /**
     * @brief A function returning the norm of a matrix located at batch_idx.
     * @param[in] batch_idx The index of the matrix in the batch. 
     * @return  The value of the matrix infinite-norm.
     */
    double norm(int batch_idx) const
    {
        int const tmp_mat_size = size();
        int const tmp_batch_size = batch_size();
        double result = 0.;
        double* vals_proxy = m_batch_matrix_csr->get_values();
        int* row_ptr_proxy = m_batch_matrix_csr->get_row_ptrs();
        Kokkos::View<double**, Kokkos::LayoutRight, typename ExecSpace::memory_space> vals_view(
                vals_proxy,
                tmp_batch_size,
                m_batch_matrix_csr->get_num_elements_per_item());
        Kokkos::View<int*, Kokkos::LayoutRight, typename ExecSpace::memory_space>
                row_ptr_view(row_ptr_proxy, tmp_mat_size + 1);

        Kokkos::parallel_reduce(
                "L-infinity norm",
                Kokkos::RangePolicy<ExecSpace>(0, tmp_mat_size),
                KOKKOS_LAMBDA(int const i, double& res) {
                    double row_sum = 0.;
                    for (int k = row_ptr_view[i]; k < row_ptr_view[i + 1]; k++) {
                        row_sum += Kokkos::abs(vals_view(batch_idx, k));
                    }
                    if (row_sum > res) {
                        res = row_sum;
                    }
                },
                Kokkos::Max<double>(result));

        return result;
    }
};

/**
     * @brief A function which converts COO data storage into CSR.
     * The results of the conversion are directly stored inside MatrixBatchCsr instance and used for computations.
     * @param[in] matrix A MatrixBatchCsr instance, provides Kokkos views to fill.
     * All input data is stored in 1D host allocated Kokkos View.
     * @param[in] vals_coo_host Values.
     * @param[in] row_coo_host Rows indices.
     * @param[in] col_coo_host Columns indices 
     */
template <MatrixBatchCsrSolver Solver = MatrixBatchCsrSolver::BICGSTAB>
void convert_coo_to_csr(
        std::unique_ptr<MatrixBatchCsr<Kokkos::DefaultExecutionSpace, Solver>>& matrix,
        Kokkos::View<double*, Kokkos::LayoutRight, Kokkos::DefaultHostExecutionSpace> vals_coo_host,
        Kokkos::View<int*, Kokkos::LayoutRight, Kokkos::DefaultHostExecutionSpace> row_coo_host,
        Kokkos::View<int*, Kokkos::LayoutRight, Kokkos::DefaultHostExecutionSpace> col_coo_host)
{
    int const mat_size = matrix->size();
    int const batch_size = matrix->batch_size();
    int const non_zero_per_system = vals_coo_host.extent(0) / batch_size;

    Kokkos::View<double**, Kokkos::LayoutRight, Kokkos::DefaultHostExecutionSpace>
            values_view_host("", batch_size, non_zero_per_system);
    Kokkos::View<int*, Kokkos::LayoutRight, Kokkos::DefaultHostExecutionSpace>
            idx_view_host("", non_zero_per_system);
    Kokkos::View<int*, Kokkos::LayoutRight, Kokkos::DefaultHostExecutionSpace>
            nnz_per_row_view_host("", mat_size + 1);
    nnz_per_row_view_host(0) = 0;
    nnz_per_row_view_host(mat_size) = non_zero_per_system;

    for (int i = 0; i < mat_size; i++) {
        int n_non_zero_in_row = 0;
        for (int k = 0; k < non_zero_per_system; k++) {
            if (row_coo_host(k) == i) {
                int const idx = nnz_per_row_view_host(i) + n_non_zero_in_row;
                idx_view_host(idx) = col_coo_host(k);
                values_view_host(0, idx) = vals_coo_host(k);
                n_non_zero_in_row++;
            }
        }
        nnz_per_row_view_host(i + 1) = nnz_per_row_view_host(i) + n_non_zero_in_row;
    }
    for (int batch_idx = 1; batch_idx < batch_size; batch_idx++) {
        for (int idx = 0; idx < non_zero_per_system; idx++) {
            values_view_host(batch_idx, idx) = vals_coo_host(batch_idx * non_zero_per_system + idx);
        }
    }
    auto [values, col_idx, nnz_per_row] = matrix->get_batch_csr();
    Kokkos::deep_copy(values, values_view_host);
    Kokkos::deep_copy(col_idx, idx_view_host);
    Kokkos::deep_copy(nnz_per_row, nnz_per_row_view_host);
}
