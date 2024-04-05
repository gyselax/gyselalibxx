// Copyright (C) The DDC development team, see COPYRIGHT.md file
//
// SPDX-License-Identifier: MIT

#pragma once

#include <algorithm>
#include <memory>
#include <optional>
#include <type_traits>
#include <utility>

#include <ginkgo/ginkgo.hpp>

#include <Kokkos_Core.hpp>

#include "ginkgo_executors.hpp"
#include "matrix.hpp"

namespace ddc::detail {

/**
 * @param gko_exec[in] A Ginkgo executor that has access to the Kokkos::View memory space
 * @param view[in] A 2-D Kokkos::View with unit stride in the second dimension
 * @return A Ginkgo Dense matrix view over the Kokkos::View data
 */
template <class KokkosViewType>
auto to_gko_dense(std::shared_ptr<const gko::Executor> const& gko_exec, KokkosViewType const& view)
{
    static_assert((Kokkos::is_view_v<KokkosViewType> && KokkosViewType::rank == 2));
    using value_type = typename KokkosViewType::traits::value_type;

    if (view.stride_1() != 1) {
        throw std::runtime_error("The view needs to be contiguous in the second dimension");
    }

    return gko::matrix::Dense<value_type>::
            create(gko_exec,
                   gko::dim<2>(view.extent(0), view.extent(1)),
                   gko::array<value_type>::view(gko_exec, view.span(), view.data()),
                   view.stride_0());
}

template <class ExecSpace>
int default_cols_per_chunk() noexcept
{
#ifdef KOKKOS_ENABLE_SERIAL
    if (std::is_same_v<ExecSpace, Kokkos::Serial>) {
        return 8192;
    }
#endif
#ifdef KOKKOS_ENABLE_OPENMP
    if (std::is_same_v<ExecSpace, Kokkos::OpenMP>) {
        return 8192;
    }
#endif
#ifdef KOKKOS_ENABLE_CUDA
    if (std::is_same_v<ExecSpace, Kokkos::Cuda>) {
        return 65535;
    }
#endif
#ifdef KOKKOS_ENABLE_HIP
    if (std::is_same_v<ExecSpace, Kokkos::HIP>) {
        return 65535;
    }
#endif
    return 1;
}

template <class ExecSpace>
unsigned int default_preconditionner_max_block_size() noexcept
{
#ifdef KOKKOS_ENABLE_SERIAL
    if (std::is_same_v<ExecSpace, Kokkos::Serial>) {
        return 32u;
    }
#endif
#ifdef KOKKOS_ENABLE_OPENMP
    if (std::is_same_v<ExecSpace, Kokkos::OpenMP>) {
        return 1u;
    }
#endif
#ifdef KOKKOS_ENABLE_CUDA
    if (std::is_same_v<ExecSpace, Kokkos::Cuda>) {
        return 1u;
    }
#endif
#ifdef KOKKOS_ENABLE_HIP
    if (std::is_same_v<ExecSpace, Kokkos::HIP>) {
        return 1u;
    }
#endif
    return 1u;
}

// Matrix class for sparse storage and iterative solve
template <class ExecSpace>
class Matrix_Sparse : public Matrix
{
    using matrix_sparse_type = gko::matrix::Csr<double, int>;
#ifdef KOKKOS_ENABLE_OPENMP
    using solver_type = std::conditional_t<
            std::is_same_v<ExecSpace, Kokkos::OpenMP>,
            gko::solver::Gmres<double>,
            gko::solver::Bicgstab<double>>;
#else
    using solver_type = gko::solver::Bicgstab<double>;
#endif


private:
    std::unique_ptr<gko::matrix::Dense<double>> m_matrix_dense;

    std::shared_ptr<matrix_sparse_type> m_matrix_sparse;

    std::shared_ptr<solver_type> m_solver;
    std::shared_ptr<gko::LinOp> m_solver_tr;

    int m_cols_per_chunk; // Maximum number of columns of B to be passed to a Ginkgo solver

    unsigned int m_preconditionner_max_block_size; // Maximum size of Jacobi-block preconditionner

public:
    // Constructor
    explicit Matrix_Sparse(
            const int mat_size,
            std::optional<int> cols_per_chunk = std::nullopt,
            std::optional<unsigned int> preconditionner_max_block_size = std::nullopt)
        : Matrix(mat_size)
        , m_cols_per_chunk(cols_per_chunk.value_or(default_cols_per_chunk<ExecSpace>()))
        , m_preconditionner_max_block_size(preconditionner_max_block_size.value_or(
                  default_preconditionner_max_block_size<ExecSpace>()))
    {
        std::shared_ptr const gko_exec = create_gko_exec<ExecSpace>();
        m_matrix_dense = gko::matrix::Dense<
                double>::create(gko_exec->get_master(), gko::dim<2>(mat_size, mat_size));
        m_matrix_dense->fill(0);
        m_matrix_sparse = matrix_sparse_type::create(gko_exec, gko::dim<2>(mat_size, mat_size));
    }

    virtual double get_element([[maybe_unused]] int i, [[maybe_unused]] int j) const override
    {
        throw std::runtime_error("MatrixSparse::get_element() is not implemented because no API is "
                                 "provided by Ginkgo");
    }

    virtual void set_element(int i, int j, double aij) override
    {
        m_matrix_dense->at(i, j) = aij;
    }

    int factorize_method() override
    {
        // Remove zeros
        gko::matrix_data<double> matrix_data(gko::dim<2>(get_size(), get_size()));
        m_matrix_dense->write(matrix_data);
        m_matrix_dense.reset();
        matrix_data.remove_zeros();
        m_matrix_sparse->read(matrix_data);
        std::shared_ptr const gko_exec = m_matrix_sparse->get_executor();

        // Create the solver factory
        std::shared_ptr const residual_criterion
                = gko::stop::ResidualNorm<double>::build().with_reduction_factor(1e-15).on(
                        gko_exec);

        std::shared_ptr const iterations_criterion
                = gko::stop::Iteration::build().with_max_iters(1000u).on(gko_exec);

        std::shared_ptr const preconditioner
                = gko::preconditioner::Jacobi<double>::build()
                          .with_max_block_size(m_preconditionner_max_block_size)
                          .on(gko_exec);

        std::unique_ptr const solver_factory
                = solver_type::build()
                          .with_preconditioner(preconditioner)
                          .with_criteria(residual_criterion, iterations_criterion)
                          .on(gko_exec);

        m_solver = solver_factory->generate(m_matrix_sparse);
        m_solver_tr = m_solver->transpose();
        gko_exec->synchronize();

        return 0;
    }

    virtual int solve_inplace_method(double* b, char transpose, int n_equations) const override
    {
        std::shared_ptr const gko_exec = m_solver->get_executor();

        int const main_chunk_size = std::min(m_cols_per_chunk, n_equations);

        Kokkos::View<double**, Kokkos::LayoutRight, ExecSpace> const
                b_view(b, get_size(), n_equations);
        Kokkos::View<double**, Kokkos::LayoutRight, ExecSpace> const
                x_view("", get_size(), main_chunk_size);

        int const iend = (n_equations + main_chunk_size - 1) / main_chunk_size;
        for (int i = 0; i < iend; ++i) {
            int const subview_begin = i * main_chunk_size;
            int const subview_end
                    = (i + 1 == iend) ? n_equations : (subview_begin + main_chunk_size);

            auto const b_subview = Kokkos::
                    subview(b_view, Kokkos::ALL, Kokkos::pair(subview_begin, subview_end));
            auto const x_subview = Kokkos::
                    subview(x_view, Kokkos::ALL, Kokkos::pair(0, subview_end - subview_begin));

            Kokkos::deep_copy(x_subview, b_subview);

            if (transpose == 'N') {
                m_solver
                        ->apply(to_gko_dense(gko_exec, b_subview),
                                to_gko_dense(gko_exec, x_subview));
            } else if (transpose == 'T') {
                m_solver_tr
                        ->apply(to_gko_dense(gko_exec, b_subview),
                                to_gko_dense(gko_exec, x_subview));
            } else {
                throw std::domain_error("transpose option not recognized");
            }


            Kokkos::deep_copy(b_subview, x_subview);
        }

        return 1;
    }
};

} // namespace ddc::detail
