// SPDX-License-Identifier: MIT
#pragma once
#include <ddc/ddc.hpp>

#include "Lagrange.hpp"
#include "i_interpolator.hpp"

/**
 * @brief A class for interpolating a function using Lagrange polynomials.
 * It is designed to work with both uniform and non-uniform mesh, and have the advantage to be local.
 *
 */
template <class DDim, BCond BcMin, BCond BcMax>
class LagrangeInterpolator : public IInterpolator<DDim>
{
    using CDim = typename DDim::continuous_dimension_type;

private:
    int m_degree;
    ddc::DiscreteDomain<DDim> m_domain;
    ddc::DiscreteVector<DDim> m_ghost;

public:
    /**
     * @brief Create a Lagrange interpolator object.
     * @param[in] degree Degree of polynomials
     * @param[in] domain Discrete domain related to direction of interest for computations
     * @param[in] ghost  Discrete vector which gives the number of ghost points. By default choose 2.
    */
    LagrangeInterpolator(
            int degree,
            ddc::DiscreteDomain<DDim> domain,
            ddc::DiscreteVector<DDim> ghost)
        : m_degree(degree)
        , m_domain(domain)
        , m_ghost(ghost)
    {
    }

    ~LagrangeInterpolator() override = default;

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
    ddc::ChunkSpan<double, ddc::DiscreteDomain<DDim>> operator()(
            ddc::ChunkSpan<double, ddc::DiscreteDomain<DDim>> const inout_data,
            ddc::ChunkSpan<const ddc::Coordinate<CDim>, ddc::DiscreteDomain<DDim>> const
                    coordinates) const override
    {
        Lagrange<Kokkos::DefaultHostExecutionSpace, DDim, BcMin, BcMax>
                evaluator(m_degree, inout_data, m_domain, m_ghost);

        ddc::for_each(
                ddc::policies::serial_host,
                coordinates.domain(),
                DDC_LAMBDA(ddc::DiscreteElement<DDim> const ix) {
                    inout_data(ix) = evaluator.evaluate(coordinates(ix));
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
template <class DDim, BCond BcMin, BCond BcMax>
class PreallocatableLagrangeInterpolator : public IPreallocatableInterpolator<DDim>
{
    LagrangeInterpolator<DDim, BcMin, BcMax> const& evaluator;

public:
    /**
     * @brief Create an object capable of creating LagrangeInterpolator objects.
     * @param[in] evaluator An operator which evaluates the value of the interpolation polynomial at requested coordinates.
     */
    explicit PreallocatableLagrangeInterpolator(
            LagrangeInterpolator<DDim, BcMin, BcMax> const& evaluator)
        : evaluator(evaluator)
    {
    }

    ~PreallocatableLagrangeInterpolator() override = default;

    /**
     * Create an instance of the LagrangeInterpolator class.
     *
     * @return A unique pointer to an instance of the LagrangeInterpolator class.
     */
    std::unique_ptr<IInterpolator<DDim>> preallocate() const override
    {
        return std::make_unique<LagrangeInterpolator<DDim, BcMin, BcMax>>(evaluator);
    }
};
