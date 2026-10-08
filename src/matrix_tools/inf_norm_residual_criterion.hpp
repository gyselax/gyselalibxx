// SPDX-License-Identifier: MIT
#pragma once

#include <source_location>

#include <ginkgo/ginkgo.hpp>

#include <Kokkos_Core.hpp>

#include "matrix_utils.hpp"

/**
 * @brief A Ginkgo stopping criterion based on the infinite norm of the residual.
 *
 * The criterion is satisfied when ||r||_inf <= reduction_factor * ||b||_inf, where r is the
 * residual vector provided by the solver at each iteration (for CG/BiCGSTAB this is the residual
 * computed by the recurrence of the method, not an explicitly re-computed b - Ax).
 * Only single right-hand sides are supported.
 *
 * @tparam ExecSpace The Kokkos execution space which has access to the Ginkgo executor memory.
 */
template <class ExecSpace = Kokkos::DefaultExecutionSpace>
class InfNormResidual
    : public gko::EnablePolymorphicObject<InfNormResidual<ExecSpace>, gko::stop::Criterion>
{
    friend class gko::EnablePolymorphicObject<InfNormResidual<ExecSpace>, gko::stop::Criterion>;

public:
    /// The parameters which can be passed to the factory.
    GKO_CREATE_FACTORY_PARAMETERS(parameters, Factory)
    {
        // Residual norm reduction factor, relative to the infinite norm of the right-hand side.
        // This defaults to 1e-8 but can be modified with .with_reduction_factor(m_tol)
        double GKO_FACTORY_PARAMETER_SCALAR(reduction_factor, 1e-8);
    };
    GKO_ENABLE_CRITERION_FACTORY(InfNormResidual<ExecSpace>, parameters, Factory);
    GKO_ENABLE_BUILD_METHOD(Factory);

    /**
     * @brief Check if the infinite norm of the residual is below the threshold.
     *
     * Note: This function should be private but it is public due to CUDA restrictions.
     *
     * @param[in] stopping_id The id of the stopping criterion, saved in the stopping status.
     * @param[in] set_finalized Controls if the current version should count as finalized.
     * @param[inout] stop_status The status of each right-hand side.
     * @param[out] one_changed Indicates if the status of a right-hand side was modified.
     * @param[in] updater The object containing the solver's current state.
     *
     * @return True if the criterion is satisfied.
     */
    bool check_impl(
            gko::uint8 stopping_id,
            bool set_finalized,
            gko::array<gko::stopping_status>* stop_status,
            bool* one_changed,
            gko::stop::Criterion::Updater const& updater) override
    {
        if (updater.ignore_residual_check_) {
            return false;
        }
        if (updater.residual_ == nullptr) {
            throw std::runtime_error(
                    "InfNormResidual requires a solver which provides the residual");
        }
        double const res_norm
                = inf_norm<ExecSpace>(gko::as<gko::matrix::Dense<double>>(updater.residual_));
        if (res_norm <= m_threshold) {
            gko::stopping_status* const status = stop_status->get_data();
            const std::source_location location = std::source_location::current();
            Kokkos::parallel_for(
                    location.function_name(),
                    Kokkos::RangePolicy<ExecSpace>(0, stop_status->get_size()),
                    KOKKOS_LAMBDA(int const i) { status[i].converge(stopping_id, set_finalized); });
            Kokkos::fence();
            *one_changed = true;
            return true;
        }
        return false;
    }

protected:
    /**
     * @brief Constructor used by gko::EnablePolymorphicObject.
     * @param[in] exec The Ginkgo executor.
     */
    explicit InfNormResidual(std::shared_ptr<gko::Executor const> exec)
        : gko::EnablePolymorphicObject<InfNormResidual<ExecSpace>, gko::stop::Criterion>(
                std::move(exec))
        , m_threshold(0.)
    {
    }

    /**
     * @brief Constructor used by the factory. Computes the convergence threshold from the right-hand side.
     * @param[in] factory The factory holding the parameters.
     * @param[in] args The arguments of the solver (system matrix, right-hand side, ...).
     */
    explicit InfNormResidual(Factory const* factory, gko::stop::CriterionArgs const& args)
        : gko::EnablePolymorphicObject<InfNormResidual<ExecSpace>, gko::stop::Criterion>(
                factory->get_executor())
        , parameters_ {factory->get_parameters()}
        , m_threshold(
                  parameters_.reduction_factor
                  * inf_norm<ExecSpace>(gko::as<gko::matrix::Dense<double>>(args.b.get())))
    {
    }

private:
    // reduction_factor * ||b||_inf
    double m_threshold;
};
