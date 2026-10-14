

# File matrix\_utils.hpp

[**File List**](files.md) **>** [**matrix\_tools**](dir_8cedd1260cc2f2819c8df2fc66ad98b5.md) **>** [**matrix\_utils.hpp**](matrix__utils_8hpp.md)

[Go to the documentation of this file](matrix__utils_8hpp.md)


```C++
// SPDX-License-Identifier: MIT
#pragma once
#include <source_location>

#include <ginkgo/extensions/kokkos/spaces.hpp>
#include <ginkgo/extensions/kokkos/types.hpp>
#include <ginkgo/ginkgo.hpp>

#include <Kokkos_Core.hpp>


enum class GkoMatrixResidualNorm { TWO, INF };

template <class KokkosViewType>
auto to_gko_multivector(
        std::shared_ptr<const gko::Executor> const& gko_exec,
        KokkosViewType const& view)
{
    static_assert((Kokkos::is_view_v<KokkosViewType> && KokkosViewType::rank == 2));
    using value_type = typename KokkosViewType::traits::value_type;

    assert(view.stride(1) == 1);
    return gko::share(
            gko::batch::MultiVector<value_type>::
                    create(gko_exec,
                           gko::batch_dim<2>(view.extent(0), gko::dim<2>(view.extent(1), 1)),
                           gko::array<value_type>::view(gko_exec, view.span(), view.data())));
}

template <class ExecSpace = Kokkos::DefaultExecutionSpace>
double inf_norm(gko::matrix::Dense<double> const* vec)
{
    // Ensure there is only 1 column
    assert(vec->get_size()[1] == 1);
    Kokkos::View<double const**, Kokkos::LayoutRight> const vec_view
            = gko::ext::kokkos::map_data(vec);
    double result = 0.;
    const std::source_location location = std::source_location::current();
    Kokkos::parallel_reduce(
            location.function_name(),
            Kokkos::RangePolicy<ExecSpace>(0, vec_view.extent(0)),
            KOKKOS_LAMBDA(int const i, double& local_max) {
                local_max = Kokkos::max(local_max, Kokkos::abs(vec_view(i, 0)));
            },
            Kokkos::Max<double>(result));
    return result;
}

inline double vector_norm(gko::matrix::Dense<double> const* vec, GkoMatrixResidualNorm const norm)
{
    if (norm == GkoMatrixResidualNorm::INF) {
        return inf_norm(vec);
    } else {
        // norm == GkoMatrixResidualNorm::TWO
        auto norm_host = gko::matrix::Dense<
                double>::create(vec->get_executor()->get_master(), gko::dim<2>(1, 1));
        vec->compute_norm2(norm_host);
        return norm_host->at(0, 0);
    }
}

inline void check_conv(
        int const batch_size,
        double const tol,
        std::shared_ptr<const gko::Executor> gko_exec,
        std::shared_ptr<const gko::batch::log::BatchConvergence<double>> logger)
{
    auto logger_residual_host
            = gko::make_temporary_clone(gko_exec->get_master(), &logger->get_residual_norm());
    bool has_converged = false;
    Kokkos::parallel_reduce(
            "convergence",
            Kokkos::RangePolicy<Kokkos::DefaultHostExecutionSpace>(0, batch_size),
            [&](int batch_idx, bool& check_tol) {
                check_tol = check_tol && (logger_residual_host->get_const_data()[batch_idx] <= tol);
            },
            Kokkos::LAnd<bool>(has_converged));
    if (!has_converged) {
        throw ::std::runtime_error("Ginkgo did not converge");
    }
}

inline void write_log(
        std::fstream& log_file,
        int const batch_index,
        int const num_iterations,
        double const implicit_res_norm,
        double const true_res_norm,
        double const b_norm,
        double const tol)
{
    log_file << " System no. " << batch_index << ":" << std::endl;
    log_file << " Number of iterations = " << num_iterations << std::endl;
    log_file << " Implicit residual norm = " << implicit_res_norm << std::endl;
    log_file << " True (Ax-b) residual norm = " << true_res_norm << std::endl;
    log_file << " Right-hand side (b) norm = " << b_norm << std::endl;
    if (!(true_res_norm <= tol)) {
        log_file << " --- System " << batch_index << " did not converge! ---" << std::endl;
    }
    log_file << "------------------------------------------------" << std::endl;
}

template <class sparse_type>
void save_logger(
        std::fstream& log_file,
        int const batch_index,
        std::unique_ptr<sparse_type> matrix,
        Kokkos::View<double*, Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace> const x_view,
        Kokkos::View<double*, Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace> const b_view,
        std::shared_ptr<const gko::log::Convergence<double>> logger,
        double const tol,
        GkoMatrixResidualNorm const norm)
{
    std::shared_ptr const gko_exec = matrix->get_executor();

    auto x = gko::matrix::Dense<double>::
            create(gko_exec,
                   gko::dim<2>(x_view.extent(0), 1),
                   gko::array<double>::view(gko_exec, x_view.span(), x_view.data()),
                   x_view.stride(0));
    auto b = gko::matrix::Dense<double>::
            create(gko_exec,
                   gko::dim<2>(b_view.extent(0), 1),
                   gko::array<double>::view(gko_exec, b_view.span(), b_view.data()),
                   b_view.stride(0));

    // allocate the residual
    auto res = gko::matrix::Dense<double>::create(gko_exec, gko::dim<2>(b_view.extent(0), 1));
    res->copy_from(b);

    gko::dim<2> norm_dim(1, 1);
    // we need constants on the device
    auto one = gko::matrix::Dense<double>::create(gko_exec, norm_dim);
    one->fill(1.0);
    auto neg_one = gko::matrix::Dense<double>::create(gko_exec, norm_dim);
    neg_one->fill(-1.0);
    //to estimate the "true" residual, the apply function below computes Ax-res, and stores the result in res.
    matrix->apply(one, x, neg_one, res);

    int const log_iters_host = logger->get_num_iterations();
    // The norm of the residual used by the solver.
    double const log_resid_host
            = vector_norm(gko::as<gko::matrix::Dense<double>>(logger->get_residual()), norm);

    write_log(
            log_file,
            batch_index,
            log_iters_host,
            log_resid_host,
            vector_norm(res.get(), norm),
            vector_norm(b.get(), norm),
            tol);
}

template <class batch_sparse_type>
void save_logger(
        std::fstream& log_file,
        std::shared_ptr<batch_sparse_type> batch_matrix,
        Kokkos::View<double**, Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace> const x_view,
        Kokkos::View<double**, Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace> const b_view,
        std::shared_ptr<const gko::batch::log::BatchConvergence<double>> logger,
        double const tol)
{
    std::shared_ptr const gko_exec = batch_matrix->get_executor();
    int const batch_size = x_view.extent(0);

    auto x = to_gko_multivector(gko_exec, x_view);
    auto b = to_gko_multivector(gko_exec, b_view);

    // allocate the residual
    auto res = gko::batch::MultiVector<double>::
            create(gko_exec, gko::batch_dim<2>(batch_size, gko::dim<2>(b_view.extent(1), 1)));
    res->copy_from(b);

    gko::batch_dim<2> norm_dim(batch_size, gko::dim<2>(1, 1));
    // allocate rhs norm on host.
    auto b_norm_host = gko::batch::MultiVector<double>::create(gko_exec->get_master(), norm_dim);
    b_norm_host->fill(0.0);
    // allocate the residual norm on host.
    auto res_norm_host = gko::batch::MultiVector<double>::create(gko_exec->get_master(), norm_dim);
    res_norm_host->fill(0.0);
    // compute rhs norm.
    b->compute_norm2(b_norm_host);
    // we need constants on the device
    auto one = gko::batch::MultiVector<double>::create(gko_exec, norm_dim);
    one->fill(1.0);
    auto neg_one = gko::batch::MultiVector<double>::create(gko_exec, norm_dim);
    neg_one->fill(-1.0);
    //to estimate the "true" residual, the apply function below computes Ax-res, and stores the result in res.
    batch_matrix->apply(one, x, neg_one, res);
    //compute residual norm.
    res->compute_norm2(res_norm_host);

    auto log_iters_host
            = gko::make_temporary_clone(gko_exec->get_master(), &logger->get_num_iterations());
    auto log_resid_host
            = gko::make_temporary_clone(gko_exec->get_master(), &logger->get_residual_norm());

    // "unbatch" converts a batch object into a vector
    // of objects of the corresponding single type.
    for (int i = 0; i < batch_size; ++i) {
        write_log(
                log_file,
                i,
                log_iters_host->get_const_data()[i],
                log_resid_host->get_const_data()[i],
                res_norm_host->create_const_view_for_item(i)->at(0, 0),
                b_norm_host->create_const_view_for_item(i)->at(0, 0),
                tol);
    }
}

template <class ExecSpace>
unsigned int default_preconditioner_max_block_size() noexcept
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
```


