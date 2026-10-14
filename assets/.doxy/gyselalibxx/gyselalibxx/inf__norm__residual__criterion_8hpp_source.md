

# File inf\_norm\_residual\_criterion.hpp

[**File List**](files.md) **>** [**matrix\_tools**](dir_8cedd1260cc2f2819c8df2fc66ad98b5.md) **>** [**inf\_norm\_residual\_criterion.hpp**](inf__norm__residual__criterion_8hpp.md)

[Go to the documentation of this file](inf__norm__residual__criterion_8hpp.md)


```C++
// SPDX-License-Identifier: MIT
#pragma once

#include <source_location>

#include <ginkgo/ginkgo.hpp>

#include <Kokkos_Core.hpp>

#include "matrix_utils.hpp"

template <class ExecSpace = Kokkos::DefaultExecutionSpace>
class InfNormResidual
    : public gko::EnablePolymorphicObject<InfNormResidual<ExecSpace>, gko::stop::Criterion>
{
    friend class gko::EnablePolymorphicObject<InfNormResidual<ExecSpace>, gko::stop::Criterion>;

public:
    GKO_CREATE_FACTORY_PARAMETERS(parameters, Factory)
    {
        // Residual norm reduction factor, relative to the infinite norm of the right-hand side.
        // This defaults to 1e-8 but can be modified with .with_reduction_factor(m_tol)
        double GKO_FACTORY_PARAMETER_SCALAR( // cppcheck-suppress unreadVariable
                reduction_factor,
                1e-8);
    };
    GKO_ENABLE_CRITERION_FACTORY(InfNormResidual<ExecSpace>, parameters, Factory);
    GKO_ENABLE_BUILD_METHOD(Factory);

    bool check_impl(
            gko::uint8 stopping_id,
            bool set_finalised,
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
                    KOKKOS_LAMBDA(int const i) { status[i].converge(stopping_id, set_finalised); });
            Kokkos::fence();
            *one_changed = true;
            return true;
        }
        return false;
    }

protected:
    explicit InfNormResidual(std::shared_ptr<gko::Executor const> exec)
        : gko::EnablePolymorphicObject<InfNormResidual<ExecSpace>, gko::stop::Criterion>(
                std::move(exec))
        , m_threshold(0.)
    {
    }

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
```


