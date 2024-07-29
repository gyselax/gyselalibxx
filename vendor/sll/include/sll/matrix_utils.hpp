#pragma once
#include <ginkgo/ginkgo.hpp>

#include <Kokkos_Core.hpp>

#include "ddc/kernels/splines/ginkgo_executors.hpp"


/**
 * @brief A function to convert a 2D Kokkos view into a ginkgo multivector structure.
 * @param gko_exec[in] A Ginkgo executor that has access to the Kokkos::View memory space
 * @param view[in] A 2-D Kokkos::View with unit stride in the second dimension
 * @return A Ginkgo Multivector view over the Kokkos::View data
 */

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

/**
 * @brief A function extracted from ginkgo to unbatch ginkgo objects.
 * @param[in] batch_object A batch of ginkgo objects ( eg. gko::batch::matrix::Dense
 *            converted to std::vector<gko::matrix::Dense> ), works also for multivectors, loggers.
 * @return a vector of ginkgo structures.
 */
template <typename InputType>
std::vector<std::unique_ptr<typename InputType::unbatch_type>> unbatch(
        const InputType* batch_object)
{
    std::vector<std::unique_ptr<typename InputType::unbatch_type>> unbatched_mats;
    for (int b = 0; b < batch_object->get_num_batch_items(); ++b) {
        unbatched_mats.emplace_back(batch_object->create_const_view_for_item(b)->clone());
    }
    return unbatched_mats;
}

/**
 * @brief A function for checking convergence. It loops over the batch and checks 
 *        the if residual is lower or equal to the prescribed tolerance.
 * @param[in] batch_size the size of the batch , ie number of linears problems.
 * @param[in] logger Ginkgo logger which contains residual and numbers of iterations
 *                   For the whole batch.
 * @param[in] gko_exec Ginkgo executor, refers to the execution space.
 * @param[in] logger Ginkgo convergence object which stores iterations number and residual for the whole batch.
 */
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
        throw ::std::runtime_error("Residual tolerance is not reached");
    }
}

/**
 * @brief A function to save convergence data using the logger.
 * @param[in] batch_matrix Batch of matrices.
 * @param[in] x_view 2d Kokkos view containing the batch of computed solutions.
 * @param[in] b_view 2d Kokkos view containing the batch of rhs. 
 * @param[in] logger Ginkgo logger which stores residual and numbers of iterations for the whole batch.
 * @param[in] tol tolerance 
 * @param[in] filename prefix used to generate the logger file.
 */
template <class batch_sparse_type>
void save_logger(
        std::shared_ptr<batch_sparse_type> batch_matrix,
        Kokkos::View<double**, Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace> const x_view,
        Kokkos::View<double**, Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace> const b_view,
        std::shared_ptr<const gko::batch::log::BatchConvergence<double>> logger,
        double const tol,
        std::string const& filename)
{
    int const batch_size = x_view.extent(0);
    std::shared_ptr const gko_exec = batch_matrix->get_executor();
    // allocate the residual
    auto res = gko::batch::MultiVector<double>::
            create(gko_exec, gko::batch_dim<2>(batch_size, gko::dim<2>(x_view.extent(1), 1)));
    res->copy_from(to_gko_multivector(gko_exec, b_view));

    auto norm_dim = gko::batch_dim<2>(batch_size, gko::dim<2>(1, 1));
    // allocate rhs norm on host.
    auto b_norm_host = gko::batch::MultiVector<double>::create(gko_exec->get_master(), norm_dim);
    b_norm_host->fill(0.0);
    // allocate the residual norm on host.
    auto res_norm_host = gko::batch::MultiVector<double>::create(gko_exec->get_master(), norm_dim);
    res_norm_host->fill(0.0);
    // compute rhs norm.
    to_gko_multivector(gko_exec, b_view)->compute_norm2(b_norm_host);
    // we need constants on the device
    auto one = gko::batch::MultiVector<double>::create(gko_exec, norm_dim);
    one->fill(1.0);
    auto neg_one = gko::batch::MultiVector<double>::create(gko_exec, norm_dim);
    neg_one->fill(-1.0);
    //to estimate the "true" residual, the apply function below computes Ax-res, and stores the result in res.
    batch_matrix->apply(one, to_gko_multivector(gko_exec, x_view), neg_one, res);
    //compute residual norm.
    res->compute_norm2(res_norm_host);

    auto log_resid_host
            = gko::make_temporary_clone(gko_exec->get_master(), &logger->get_residual_norm());
    auto log_iters_host
            = gko::make_temporary_clone(gko_exec->get_master(), &logger->get_num_iterations());

    std::fstream log_file(filename + "_logger.txt", std::ios::out | std::ios::app);
    // "unbatch" converts a batch object into a vector
    // of objects of the corresponding single type.
    auto unb_res_norm = unbatch(res_norm_host.get());
    auto unb_bnorm = unbatch(b_norm_host.get());
    for (int i = 0; i < batch_size; ++i) {
        // Logger  output
        log_file << " System no. " << i << ": Ax-b residual norm = " << unb_res_norm[i]->at(0, 0)
                 << ", implicit residual norm = " << log_resid_host->get_const_data()[i]
                 << ", iterations = " << log_iters_host->get_const_data()[i] << std::endl;
        log_file << " unbatched bnorm at(i,0)" << unb_bnorm[i]->at(0, 0) << std::endl;
        log_file << " unbatched residual norm at(i,0)" << unb_res_norm[i]->at(0, 0) << std::endl;
        if (!(unb_res_norm[i]->at(0, 0) <= tol)) {
            log_file << "System " << i << " converged only to " << unb_res_norm[i]->at(0, 0)
                     << " relative residual." << std::endl;
        }
        log_file << "------------------------------------------------" << std::endl;
    }
    log_file.close();
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