// SPDX-License-Identifier: MIT

#pragma once

#include <ddc/kernels/splines.hpp>

#include "iinterpolator.hpp"

/**
 * @brief A class for interpolating a function using splines.
 *
 * The class is parametrised by multiple template parameters. Please note that CTAD will deduce all these
 * template parameters from the Builder and Evaluator passed as constructor arguments.
 *
 * @tparam DDimI The dimension of interest.
 * @tparam BSplines The BSplines along the dimension of interest.
 * @tparam BcMin The boundary condition at the lower boundary.
 * @tparam BcMax The boundary condition at the upper boundary.
 * @tparam Grid1D... All the dimensions of the interpolation problem (batched + interpolated).
 */
template <
        class DDimI,
        class BSplines,
        ddc::BoundCond BcMin,
        ddc::BoundCond BcMax,
        class LeftExtrapolationRule,
        class RightExtrapolationRule,
        ddc::SplineSolver Solver,
        class... Grid1D>
class SplineInterpolator : public IInterpolator<DDimI, Grid1D...>
{
    using BuilderType = ddc::SplineBuilder<
            Kokkos::DefaultExecutionSpace,
            Kokkos::DefaultExecutionSpace::memory_space,
            BSplines,
            DDimI,
            BcMin,
            BcMax,
            Solver,
            Grid1D...>;
    using EvaluatorType = ddc::SplineEvaluator<
            Kokkos::DefaultExecutionSpace,
            Kokkos::DefaultExecutionSpace::memory_space,
            BSplines,
            DDimI,
            LeftExtrapolationRule,
            RightExtrapolationRule,
            Grid1D...>;
    using deriv_type = typename IInterpolator<DDimI, Grid1D...>::deriv_type;
    using batched_derivs_domain_type =
            typename IInterpolator<DDimI, Grid1D...>::batched_derivs_domain_type;

private:
    BuilderType const& m_builder;

    EvaluatorType const& m_evaluator;

    mutable ddc::Chunk<
            double,
            typename BuilderType::batched_spline_domain_type,
            ddc::DeviceAllocator<double>>
            m_coefs;


public:
    /**
     * @brief Create a spline interpolator object.
     * @param[in] builder An operator which builds spline coefficients from the values of a function at known interpolation points.
     * @param[in] evaluator An operator which evaluates the value of a spline at requested coordinates.
     */
    SplineInterpolator(BuilderType const& builder, EvaluatorType const& evaluator)
        : m_builder(builder)
        , m_evaluator(evaluator)
        , m_coefs(builder.batched_spline_domain())
    {
    }

    ~SplineInterpolator() override = default;

    batched_derivs_domain_type batched_derivs_domain_xmin(
            ddc::DiscreteDomain<Grid1D...> dom) const override
    {
        return ddc::replace_dim_of<DDimI, deriv_type>(
                dom,
                ddc::DiscreteDomain<deriv_type>(
                        ddc::DiscreteElement<deriv_type>(1),
                        ddc::DiscreteVector<deriv_type>(BuilderType::s_nbc_xmin)));
    }

    batched_derivs_domain_type batched_derivs_domain_xmax(
            ddc::DiscreteDomain<Grid1D...> dom) const override
    {
        return ddc::replace_dim_of<DDimI, deriv_type>(
                dom,
                ddc::DiscreteDomain<deriv_type>(
                        ddc::DiscreteElement<deriv_type>(1),
                        ddc::DiscreteVector<deriv_type>(BuilderType::s_nbc_xmax)));
    }

    /**
     * @brief Approximate the value of a function at a set of coordinates using the
     * current values at a known set of interpolation points.
     *
     * @param[in, out] inout_data On input: an array containing the value of the function at the interpolation points.
     *           On output: an array containing the value of the function at the coordinates.
     * @param[in] coordinates The coordinates where the function should be evaluated.
     * @param[in] derivs_xmin The values of the derivatives at the lower boundary
     * (used only with ddc::BoundCond::HERMITE lower boundary condition).
     * @param[in] derivs_xmax The values of the derivatives at the upper boundary
     * (used only with ddc::BoundCond::HERMITE upper boundary condition).
     *
     * @return A reference to the inout_data array containing the value of the function at the coordinates.
     */
    device_t<ddc::ChunkSpan<double, ddc::DiscreteDomain<Grid1D...>>> operator()(
            device_t<ddc::ChunkSpan<double, ddc::DiscreteDomain<Grid1D...>>> const inout_data,
            device_t<ddc::ChunkSpan<
                    const ddc::Coordinate<typename DDimI::continuous_dimension_type>,
                    ddc::DiscreteDomain<Grid1D...>>> const coordinates,
            std::optional<device_t<
                    ddc::ChunkSpan<double const, typename BuilderType::batched_derivs_domain_type>>>
                    derivs_xmin
            = std::nullopt,
            std::optional<device_t<
                    ddc::ChunkSpan<double const, typename BuilderType::batched_derivs_domain_type>>>
                    derivs_xmax
            = std::nullopt) const override
    {
        m_builder(m_coefs.span_view(), inout_data.span_cview(), derivs_xmin, derivs_xmax);
        m_evaluator(inout_data, coordinates, m_coefs.span_cview());
        return inout_data;
    }
};

/**
 * @brief A class which stores information necessary to create an instance of the SplineInterpolator class.
 *
 * This class allows an instance of the SplineInterpolator class where necessary. This allows the
 * memory allocated in the private members of the SplineInterpolator to be freed when the object is not in use.
 * These objects are: m_coefs, m_derivs_min_alloc, m_derivs_max_alloc.
 */
template <
        class DDimI,
        class BSplines,
        ddc::BoundCond BcMin,
        ddc::BoundCond BcMax,
        class LeftExtrapolationRule,
        class RightExtrapolationRule,
        ddc::SplineSolver Solver,
        class... Grid1D>
class PreallocatableSplineInterpolator : public IPreallocatableInterpolator<DDimI, Grid1D...>
{
    using BuilderType = ddc::SplineBuilder<
            Kokkos::DefaultExecutionSpace,
            Kokkos::DefaultExecutionSpace::memory_space,
            BSplines,
            DDimI,
            BcMin,
            BcMax,
            Solver,
            Grid1D...>;
    using EvaluatorType = ddc::SplineEvaluator<
            Kokkos::DefaultExecutionSpace,
            Kokkos::DefaultExecutionSpace::memory_space,
            BSplines,
            DDimI,
            LeftExtrapolationRule,
            RightExtrapolationRule,
            Grid1D...>;

    BuilderType const& m_builder;

    EvaluatorType const& m_evaluator;

public:
    /**
     * @brief Create an object capable of creating SplineInterpolator objects.
     * @param[in] builder An operator which builds spline coefficients from the values of a function at known interpolation points.
     * @param[in] evaluator An operator which evaluates the value of a spline at requested coordinates.
     */
    PreallocatableSplineInterpolator(BuilderType const& builder, EvaluatorType const& evaluator)
        : m_builder(builder)
        , m_evaluator(evaluator)
    {
    }

    ~PreallocatableSplineInterpolator() override = default;

    /**
     * Create an instance of the SplineInterpolator class.
     *
     * @return A unique pointer to an instance of the SplineInterpolator class.
     */
    std::unique_ptr<IInterpolator<DDimI, Grid1D...>> preallocate() const override
    {
        return std::make_unique<SplineInterpolator<
                DDimI,
                BSplines,
                BcMin,
                BcMax,
                LeftExtrapolationRule,
                RightExtrapolationRule,
                Solver,
                Grid1D...>>(m_builder, m_evaluator);
    }
};
