// SPDX-License-Identifier: MIT
#pragma once
#include <ddc/ddc.hpp>

#include "Lagrange.hpp"
#include "i_interpolator_batched.hpp"

/**
 * @brief A class for interpolating a function using Lagrange polynomials.
 * It is designed to work with both uniform and non-uniform mesh, and have the advantage to be local.
 *
 */
template <class DDimI, BCond BcMin, BCond BcMax, class... DDim>
class LagrangeInterpolatorBatched : public IInterpolatorBatched<DDimI, DDim...>
{
    using CDim = typename DDimI::continuous_dimension_type;

private:
    int m_degree;
    ddc::DiscreteVector<DDimI> m_ghost;

public:
    /**
     * @brief Create a Batched Lagrange interpolator object.
     * @param[in] degree Degree of polynomials
     * @param[in] ghost  Discrete vector which gives the number of ghost points. By default choose 2.
    */
    LagrangeInterpolatorBatched(int degree, ddc::DiscreteVector<DDimI> ghost)
        : m_degree(degree)
        , m_ghost(ghost)
    {
    }

    ~LagrangeInterpolatorBatched() override = default;

    /**
     * @brief Approximate the value of a function at a set of coordinates using the
     * current values at a known set of interpolation points.
     *
     * @param[in, out] inout_data On input: an array containing the value of the function at the interpolation points.
     *                   On output: an array containing the value of the function at the coordinates.
     * @param[in] coordinates The coordinates where the function should be evaluated.
     *
     * @return A reference to the inout_data array containing the value of the function at the coordinates.
     */
    device_t<ddc::ChunkSpan<double, ddc::DiscreteDomain<DDim...>>> operator()(
            device_t<ddc::ChunkSpan<double, ddc::DiscreteDomain<DDim...>>> inout_data,
            device_t<ddc::ChunkSpan<
                    const ddc::Coordinate<typename DDimI::continuous_dimension_type>,
                    ddc::DiscreteDomain<DDim...>>> coordinates) const override
    {
        static_assert(
                BcMin != BCond::PERIODIC,
                "PERIODIC Boundary condition is not supported yet in LagrangeInterpolatorBatched.");

        int const deg = m_degree;
        auto const ghost = m_ghost;
        auto inout_data_tmp_alloc
                = ddc::create_mirror_and_copy(Kokkos::DefaultExecutionSpace(), inout_data);
        auto inout_data_tmp = inout_data_tmp_alloc.span_view();
        auto batch_domain
                = ddc::remove_dims_of(inout_data.domain(), inout_data.template domain<DDimI>());
        ddc::for_each(
                ddc::policies::parallel_device,
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
 * @brief A class which stores information necessary to create an instance of the LagrangeInterpolatorBatched class.
 *
 * This class allows an instance of the LagrangeInterpolatorBatched class where necessary. This allows the
 * memory allocated in the private members of the InterpolatorBatched to be freed when the object is not in use.
 */
template <class DDimI, BCond BcMin, BCond BcMax, class... DDim>
class PreallocatableLagrangeInterpolatorBatched
    : public IPreallocatableInterpolatorBatched<DDimI, DDim...>
{
    LagrangeInterpolatorBatched<DDimI, BcMin, BcMax, DDim...> const& m_evaluator;

public:
    /**
     * @brief Create an object capable of creating LagrangeInterpolatorBatched objects.
     * @param[in] evaluator An operator which evaluates the value of the interpolation polynomial at requested coordinates.
     */
    explicit PreallocatableLagrangeInterpolatorBatched(
            LagrangeInterpolatorBatched<DDimI, BcMin, BcMax, DDim...> const& evaluator)
        : m_evaluator(evaluator)
    {
    }

    ~PreallocatableLagrangeInterpolatorBatched() override = default;

    /**
     * Create an instance of the LagrangeInterpolatorBatched class.
     *
     * @return A unique pointer to an instance of the LagrangeInterpolatorBatched class.
     */
    std::unique_ptr<IInterpolatorBatched<DDimI, DDim...>> preallocate() const override
    {
        return std::make_unique<LagrangeInterpolatorBatched<DDimI, BcMin, BcMax, DDim...>>(
                m_evaluator);
    }
};
