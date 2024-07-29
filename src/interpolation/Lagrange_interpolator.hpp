// SPDX-License-Identifier: MIT
#pragma once
#include <ddc/ddc.hpp>

#include "Lagrange.hpp"
#include "iinterpolator.hpp"

/**
 * @brief A class for interpolating a function using Lagrange polynomials.
 * It is designed to work with both uniform and non-uniform mesh, and have the advantage to be local.
 *
 */
template <class DDimI, BCond BcMin, BCond BcMax, class... Grid1D>
class LagrangeInterpolator : public IInterpolator<DDimI, Grid1D...>
{
    using CDim = typename DDimI::continuous_dimension_type;
    using deriv_type = typename IInterpolator<DDimI, Grid1D...>::deriv_type;
    using batched_derivs_domain_type =
            typename IInterpolator<DDimI, Grid1D...>::batched_derivs_domain_type;

private:
    int m_degree;
    ddc::DiscreteVector<DDimI> m_ghost;

public:
    /**
     * @brief Create a  Lagrange interpolator object.
     * @param[in] degree Degree of polynomials
     * @param[in] ghost  Discrete vector which gives the number of ghost points. By default choose 2.
    */
    LagrangeInterpolator(int degree, ddc::DiscreteVector<DDimI> ghost)
        : m_degree(degree)
        , m_ghost(ghost)
    {
    }

    ~LagrangeInterpolator() override = default;

    batched_derivs_domain_type batched_derivs_domain_xmin(
            ddc::DiscreteDomain<Grid1D...> dom) const override
    {
        return ddc::replace_dim_of<DDimI, deriv_type>(
                dom,
                ddc::DiscreteDomain<deriv_type>(
                        ddc::DiscreteElement<deriv_type>(1),
                        ddc::DiscreteVector<deriv_type>(0)));
    }

    batched_derivs_domain_type batched_derivs_domain_xmax(
            ddc::DiscreteDomain<Grid1D...> dom) const override
    {
        return ddc::replace_dim_of<DDimI, deriv_type>(
                dom,
                ddc::DiscreteDomain<deriv_type>(
                        ddc::DiscreteElement<deriv_type>(1),
                        ddc::DiscreteVector<deriv_type>(0)));
    }

    /**
     * @brief Approximate the value of a function at a set of coordinates using the
     * current values at a known set of interpolation points.
     *
     * @param[in, out] inout_data On input: an array containing the value of the function at the interpolation points.
     *                   On output: an array containing the value of the function at the coordinates.
     * @param[in] coordinates The coordinates where the function should be evaluated.
     * @param[in] derivs_xmin The values of the derivatives at the lower boundary (unused in this class).
     * @param[in] derivs_xmax The values of the derivatives at the upper boundary (unused in this class).
     *
     * @return A reference to the inout_data array containing the value of the function at the coordinates.
     */
    device_t<ddc::ChunkSpan<double, ddc::DiscreteDomain<Grid1D...>>> operator()(
            device_t<ddc::ChunkSpan<double, ddc::DiscreteDomain<Grid1D...>>> inout_data,
            device_t<ddc::ChunkSpan<
                    const ddc::Coordinate<typename DDimI::continuous_dimension_type>,
                    ddc::DiscreteDomain<Grid1D...>>> coordinates,
            [[maybe_unused]] std::optional<device_t<ddc::ChunkSpan<
                    double const,
                    typename IInterpolator<DDimI, Grid1D...>::batched_derivs_domain_type>>>
                    derivs_xmin
            = std::nullopt,
            [[maybe_unused]] std::optional<device_t<ddc::ChunkSpan<
                    double const,
                    typename IInterpolator<DDimI, Grid1D...>::batched_derivs_domain_type>>>
                    derivs_xmax
            = std::nullopt) const override
    {
        static_assert(
                BcMin != BCond::PERIODIC,
                "PERIODIC Boundary condition is not supported yet in LagrangeInterpolator.");

        int const deg = m_degree;
        auto const ghost = m_ghost;
        auto inout_data_tmp_alloc
                = ddc::create_mirror_and_copy(Kokkos::DefaultExecutionSpace(), inout_data);
        auto inout_data_tmp = inout_data_tmp_alloc.span_view();
        auto batch_domain = ddc::remove_dims_of<DDimI>(inout_data.domain());
        ddc::parallel_for_each(
                Kokkos::DefaultExecutionSpace(),
                batch_domain,
                KOKKOS_LAMBDA(typename decltype(batch_domain)::discrete_element_type const i) {
                    Lagrange<Kokkos::DefaultExecutionSpace, DDimI, BcMin, BcMax> evaluator(
                            deg,
                            inout_data_tmp[i],
                            inout_data.template domain<DDimI>(),
                            ghost);
                    for (ddc::DiscreteElement<DDimI> j : inout_data.template domain<DDimI>()) {
                        inout_data(i, j) = evaluator.evaluate(coordinates(i, j));
                    }
                });
        return inout_data;
    }
};

/**
 * @brief A class which stores information necessary to create an instance of the LagrangeInterpolator class.
 *
 * This class allows an instance of the LagrangeInterpolator class where necessary. This allows the
 * memory allocated in the private members of the Interpolator to be freed when the object is not in use.
 */
template <class DDimI, BCond BcMin, BCond BcMax, class... Grid1D>
class PreallocatableLagrangeInterpolator : public IPreallocatableInterpolator<DDimI, Grid1D...>
{
    LagrangeInterpolator<DDimI, BcMin, BcMax, Grid1D...> const& m_evaluator;

public:
    /**
     * @brief Create an object capable of creating LagrangeInterpolator objects.
     * @param[in] evaluator An operator which evaluates the value of the interpolation polynomial at requested coordinates.
     */
    explicit PreallocatableLagrangeInterpolator(
            LagrangeInterpolator<DDimI, BcMin, BcMax, Grid1D...> const& evaluator)
        : m_evaluator(evaluator)
    {
    }

    ~PreallocatableLagrangeInterpolator() override = default;

    /**
     * Create an instance of the LagrangeInterpolator class.
     *
     * @return A unique pointer to an instance of the LagrangeInterpolator class.
     */
    std::unique_ptr<IInterpolator<DDimI, Grid1D...>> preallocate() const override
    {
        return std::make_unique<LagrangeInterpolator<DDimI, BcMin, BcMax, Grid1D...>>(m_evaluator);
    }
};
