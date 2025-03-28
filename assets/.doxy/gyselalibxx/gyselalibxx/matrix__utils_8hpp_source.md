

# File matrix\_utils.hpp

[**File List**](files.md) **>** [**matrix\_tools**](dir_8cedd1260cc2f2819c8df2fc66ad98b5.md) **>** [**matrix\_utils.hpp**](matrix__utils_8hpp.md)

[Go to the documentation of this file](matrix__utils_8hpp.md)


```C++
// SPDX-License-Identifier: MIT
#pragma once
#include <ginkgo/ginkgo.hpp>

#include <Kokkos_Core.hpp>

template <class KokkosViewType>
auto to_gko_multivector(
        std::shared_ptr<const gko::Executor> const& gko_exec,
        KokkosViewType const& view)
{
    static_assert((Kokkos::is_view_v<KokkosViewType> && KokkosViewType::rank == 2));
    using value_type = typename KokkosViewType::traits::value_type;

    assert(view.stride_1() == 1);
    return gko::share(
            gko::batch::MultiVector<value_type>::
                    create(gko_exec,
                           gko::batch_dim<2>(view.extent(0), gko::dim<2>(view.extent(1), 1)),
                           gko::array<value_type>::view(gko_exec, view.span(), view.data())));
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
        double const tol)
{
    std::shared_ptr const gko_exec = matrix->get_executor();

    auto x = gko::matrix::Dense<double>::
            create(gko_exec,
                   gko::dim<2>(x_view.extent(0), 1),
                   gko::array<double>::view(gko_exec, x_view.span(), x_view.data()),
                   x_view.stride_0());
    auto b = gko::matrix::Dense<double>::
            create(gko_exec,
                   gko::dim<2>(b_view.extent(0), 1),
                   gko::array<double>::view(gko_exec, b_view.span(), b_view.data()),
                   b_view.stride_0());

    // allocate the residual
    auto res = gko::matrix::Dense<double>::create(gko_exec, gko::dim<2>(b_view.extent(0), 1));
    res->copy_from(b);

    gko::dim<2> norm_dim(1, 1);
    // allocate rhs norm on host.
    auto b_norm_host = gko::matrix::Dense<double>::create(gko_exec->get_master(), norm_dim);
    b_norm_host->fill(0.0);
    // allocate the residual norm on host.
    auto res_norm_host = gko::matrix::Dense<double>::create(gko_exec->get_master(), norm_dim);
    res_norm_host->fill(0.0);
    // compute rhs norm.
    b->compute_norm2(b_norm_host);
    // we need constants on the device
    auto one = gko::matrix::Dense<double>::create(gko_exec, norm_dim);
    one->fill(1.0);
    auto neg_one = gko::matrix::Dense<double>::create(gko_exec, norm_dim);
    neg_one->fill(-1.0);
    //to estimate the "true" residual, the apply function below computes Ax-res, and stores the result in res.
    matrix->apply(one, x, neg_one, res);
    //compute residual norm.
    res->compute_norm2(res_norm_host);

    int const log_iters_host = logger->get_num_iterations();
    double const log_resid_host
            = gko::make_temporary_clone(
                      gko_exec->get_master(),
                      gko::as<gko::matrix::Dense<double>>(logger->get_residual_norm()))
                      ->at(0, 0);

    write_log(
            log_file,
            batch_index,
            log_iters_host,
            log_resid_host,
            res_norm_host->at(0, 0),
            b_norm_host->at(0, 0),
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
```


