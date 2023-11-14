#pragma once
#include <algorithm>
#include <cassert>
#include <iomanip>
#include <iostream>
#include <memory>
#include <type_traits>
#include <utility>

#include <ginkgo/ginkgo.hpp>

#include <Kokkos_Core.hpp>

#include "ginkgo/core/matrix/dense.hpp"

#include "matrix.hpp"
#include "view.hpp"

namespace ddc::detail {

// Matrix class for Csr storage and iterative solve
template <class ExecSpace>
class Matrix_Sparse : public Matrix
{
private:
    const int m_m;
    const int m_n;
    Kokkos::View<int*, Kokkos::HostSpace> m_rows;
    Kokkos::View<int*, Kokkos::HostSpace> m_cols;
    Kokkos::View<double*, Kokkos::HostSpace> m_data;

    std::unique_ptr<
            gko::solver::Bicgstab<gko::default_precision>::Factory,
            std::default_delete<gko::solver::Bicgstab<gko::default_precision>::Factory>>
            m_solver_factory;

    int m_cols_per_par_chunk; // Maximum number of columns of B to be passed to a Ginkgo solver
    int m_par_chunks_per_seq_chunk; // Maximum number of teams to be executed in parallel
    int m_preconditionner_max_block_size; // Maximum size of Jacobi-block preconditionner

public:
    // Constructor
    explicit Matrix_Sparse(
            const int mat_size,
            std::optional<int> cols_per_par_chunk = std::nullopt,
            std::optional<int> par_chunks_per_seq_chunk = std::nullopt,
            std::optional<unsigned int> preconditionner_max_block_size = std::nullopt)
        : Matrix(mat_size)
        , m_m(mat_size)
        , m_n(mat_size)
        , m_rows("rows", mat_size + 1)
        , m_cols("cols", mat_size * mat_size)
        , m_data("data", mat_size * mat_size)
    {
        // Fill the csr indexes as a dense matrix and initialize with zeros (zeros will be removed once non-zeros elements will be set)
        for (int i = 0; i < m_m * m_n; i++) {
            if (i < m_m + 1) {
                m_rows(i) = i * m_n; //CSR
            }
            m_cols(i) = i % m_n;
            m_data(i) = 0;
        }

        if (cols_per_par_chunk.has_value()) {
            m_cols_per_par_chunk = cols_per_par_chunk.value();
        } else {
#ifdef KOKKOS_ENABLE_SERIAL
            if (std::is_same_v<ExecSpace, Kokkos::Serial>) {
                m_cols_per_par_chunk = 512;
            }
#endif
#ifdef KOKKOS_ENABLE_OPENMP
            if (std::is_same_v<ExecSpace, Kokkos::OpenMP>) {
                m_cols_per_par_chunk = 512;
            }
#endif
#ifdef KOKKOS_ENABLE_CUDA
            if (std::is_same_v<ExecSpace, Kokkos::Cuda>) {
                m_cols_per_par_chunk = 65535;
            }
#endif
#ifdef KOKKOS_ENABLE_HIP
            if (std::is_same_v<ExecSpace, Kokkos::HIP>) {
                m_cols_per_par_chunk = 65535;
            }
#endif
        }

        if (par_chunks_per_seq_chunk.has_value()) {
            m_par_chunks_per_seq_chunk = par_chunks_per_seq_chunk.value();
        } else {
#ifdef KOKKOS_ENABLE_SERIAL
            if (std::is_same_v<ExecSpace, Kokkos::Serial>) {
                m_par_chunks_per_seq_chunk = 1;
            }
#endif
#ifdef KOKKOS_ENABLE_OPENMP
            if (std::is_same_v<ExecSpace, Kokkos::OpenMP>) {
                m_par_chunks_per_seq_chunk = Kokkos::DefaultHostExecutionSpace().concurrency();
            }
#endif
#ifdef KOKKOS_ENABLE_CUDA
            if (std::is_same_v<ExecSpace, Kokkos::Cuda>) {
                m_par_chunks_per_seq_chunk = 1;
            }
#endif
#ifdef KOKKOS_ENABLE_HIP
            if (std::is_same_v<ExecSpace, Kokkos::HIP>) {
                m_par_chunks_per_seq_chunk = 1;
            }
#endif
        }

        if (preconditionner_max_block_size.has_value()) {
            m_preconditionner_max_block_size = preconditionner_max_block_size.value();
        } else {
#ifdef KOKKOS_ENABLE_SERIAL
            if (std::is_same_v<ExecSpace, Kokkos::Serial>) {
                m_preconditionner_max_block_size = 8u;
            }
#endif
#ifdef KOKKOS_ENABLE_OPENMP
            if (std::is_same_v<ExecSpace, Kokkos::OpenMP>) {
                m_preconditionner_max_block_size = 8u;
            }
#endif
#ifdef KOKKOS_ENABLE_CUDA
            if (std::is_same_v<ExecSpace, Kokkos::Cuda>) {
                m_preconditionner_max_block_size = 1u;
            }
#endif
#ifdef KOKKOS_ENABLE_HIP
            if (std::is_same_v<ExecSpace, Kokkos::HIP>) {
                m_preconditionner_max_block_size = 1u;
            }
#endif
        }

        // Create the solver factory
        std::shared_ptr<gko::Executor> gko_exec;
        if (false) {
        }
#ifdef KOKKOS_ENABLE_OPENMP
        else if (std::is_same_v<ExecSpace, Kokkos::OpenMP>) {
            gko_exec = create_gko_exec<Kokkos::Serial>();
        }
#endif
        else {
            gko_exec = create_gko_exec<ExecSpace>();
        }
        std::shared_ptr<gko::stop::ResidualNorm<>::Factory> residual_criterion
                = gko::stop::ResidualNorm<>::build().with_reduction_factor(1e-20).on(gko_exec);
        m_solver_factory
                = gko::solver::Bicgstab<>::build()
                          .with_preconditioner(
                                  gko::preconditioner::Jacobi<>::build()
                                          .with_max_block_size(m_preconditionner_max_block_size)
                                          .on(gko_exec))
                          .with_criteria(
                                  residual_criterion,
                                  gko::stop::Iteration::build().with_max_iters(1000u).on(gko_exec))
                          .on(gko_exec);
    }

    std::unique_ptr<gko::matrix::Dense<>, std::default_delete<gko::matrix::Dense<>>> to_gko_vec(
            double* vec_ptr,
            size_t n,
            size_t n_equations,
            std::shared_ptr<gko::Executor> gko_exec) const
    {
        auto v = gko::matrix::Dense<>::
                create(gko_exec,
                       gko::dim<2>(n, n_equations),
                       gko::array<double>::view(gko_exec, n * n_equations, vec_ptr),
                       n_equations);
        return v;
    }

    std::unique_ptr<gko::matrix::Csr<>, std::default_delete<gko::matrix::Csr<>>> to_gko_mat(
            double* mat_ptr,
            size_t n_nonzero_rows,
            size_t n_nonzeros,
            std::shared_ptr<gko::Executor> gko_exec) const
    {
        auto M = gko::matrix::Csr<>::
                create(gko_exec,
                       gko::dim<2>(m_m, m_n),
                       gko::array<double>::view(gko_exec, n_nonzeros, mat_ptr),
                       gko::array<int>::view(gko_exec, n_nonzeros, m_cols.data()),
                       gko::array<int>::view(gko_exec, n_nonzero_rows + 1, m_rows.data()));
        return M;
    }

    virtual double get_element(int i, int j) const override
    {
        throw std::runtime_error("MatrixSparse::get_element() is not implemented because no API is "
                                 "provided by Ginkgo");
    }

    virtual void set_element(int i, int j, double aij) override
    {
        m_data(i * m_n + j) = aij;
    }

    int factorize_method() override
    {
        std::shared_ptr<gko::Executor> gko_exec = create_gko_exec<ExecSpace>();
        // Remove zeros
        auto data_mat
                = gko::share(to_gko_mat(m_data.data(), m_m, m_m * m_n, gko_exec->get_master()));
        auto data_mat_ = gko::matrix_data<>(gko::dim<2>(m_m, m_n));
        data_mat->write(data_mat_);
        data_mat_.remove_zeros();
        data_mat->read(data_mat_);

        // Realloc Kokkos::Views without zeros
        Kokkos::realloc(Kokkos::WithoutInitializing, m_cols, data_mat_.nonzeros.size());
        Kokkos::realloc(Kokkos::WithoutInitializing, m_data, data_mat_.nonzeros.size());
        Kokkos::deep_copy(
                m_rows,
                Kokkos::View<
                        int*,
                        Kokkos::HostSpace,
                        Kokkos::MemoryTraits<
                                Kokkos::Unmanaged>>(data_mat->get_row_ptrs(), m_m + 1));
        Kokkos::deep_copy(
                m_cols,
                Kokkos::View<int*, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>(
                        data_mat->get_col_idxs(),
                        data_mat->get_num_stored_elements()));
        Kokkos::deep_copy(
                m_data,
                Kokkos::View<double*, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>(
                        data_mat->get_values(),
                        data_mat->get_num_stored_elements()));

        return 0;
    }

    virtual int solve_inplace_method(double* b, char transpose, int n_equations) const override
    {
        std::shared_ptr<gko::Executor> gko_exec = create_gko_exec<ExecSpace>();
        auto data_mat = gko::share(to_gko_mat(
                m_data.data(),
                m_rows.size() - 1,
                m_cols.size(),
                gko_exec->get_master()));
        auto data_mat_device = gko::share(gko::clone(gko_exec, data_mat));
        Kokkos::View<
                double**,
                Kokkos::LayoutRight,
                ExecSpace,
                Kokkos::MemoryTraits<Kokkos::Unmanaged>>
                b_view(b, m_n, n_equations);

        const int n_seq_chunks
                = n_equations / m_cols_per_par_chunk / m_par_chunks_per_seq_chunk + 1;
        const int par_chunks_per_last_seq_chunk
                = (n_equations % (m_cols_per_par_chunk * m_par_chunks_per_seq_chunk))
                          / m_cols_per_par_chunk
                  + 1;
        const int cols_per_last_par_chunk
                = (n_equations % (m_cols_per_par_chunk * m_par_chunks_per_seq_chunk * n_seq_chunks))
                  % m_cols_per_par_chunk;

        Kokkos::View<double***, Kokkos::LayoutRight, ExecSpace> b_buffer(
                "b_buffer",
                std::min(n_equations / m_cols_per_par_chunk, m_par_chunks_per_seq_chunk),
                m_n,
                std::min(n_equations, m_cols_per_par_chunk));
        // Last par_chunk of last seq_chunk do not have same number of columns than the others. To get proper layout (because we passe the pointers to Ginkgo), we need a dedicated allocation
        Kokkos::View<double**, Kokkos::LayoutRight, ExecSpace>
                b_last_buffer("b_last_buffer", m_n, cols_per_last_par_chunk);

        for (int i = 0; i < n_seq_chunks; i++) {
            int n_par_chunks_in_seq_chunk = i < n_seq_chunks - 1 ? m_par_chunks_per_seq_chunk
                                                                 : par_chunks_per_last_seq_chunk;
            Kokkos::parallel_for(
                    Kokkos::RangePolicy<
                            Kokkos::DefaultHostExecutionSpace>(0, n_par_chunks_in_seq_chunk),
                    [&](int const j) {
                        int n_equations_in_par_chunk
                                = (i < n_seq_chunks - 1 || j < n_par_chunks_in_seq_chunk - 1)
                                          ? m_cols_per_par_chunk
                                          : cols_per_last_par_chunk;
                        if (n_equations_in_par_chunk != 0) {
                            auto solver = m_solver_factory->generate(data_mat_device);
                            std::pair<int, int> par_chunk_window(
                                    (i * m_par_chunks_per_seq_chunk + j) * m_cols_per_par_chunk,
                                    (i * m_par_chunks_per_seq_chunk + j) * m_cols_per_par_chunk
                                            + n_equations_in_par_chunk);
                            Kokkos::View<double**, Kokkos::LayoutRight, ExecSpace> b_par_chunk;
                            if (i < n_seq_chunks - 1 || j < n_par_chunks_in_seq_chunk - 1) {
                                b_par_chunk = Kokkos::
                                        subview(b_buffer,
                                                j,
                                                Kokkos::ALL,
                                                std::pair<int, int>(0, n_equations_in_par_chunk));
                            } else {
                                b_par_chunk = Kokkos::
                                        subview(b_last_buffer,
                                                Kokkos::ALL,
                                                std::pair<int, int>(0, n_equations_in_par_chunk));
                            }
                            Kokkos::deep_copy(
                                    b_par_chunk,
                                    Kokkos::subview(b_view, Kokkos::ALL, par_chunk_window));
                            auto b_vec_batch = to_gko_vec(
                                    b_par_chunk.data(),
                                    m_n,
                                    n_equations_in_par_chunk,
                                    gko_exec);

                            solver->apply(b_vec_batch, b_vec_batch); // inplace solve
                            Kokkos::deep_copy(
                                    Kokkos::subview(b_view, Kokkos::ALL, par_chunk_window),
                                    b_par_chunk);
                        }
                    });
        }
        return 1;
    }
};

} // namespace ddc::detail
