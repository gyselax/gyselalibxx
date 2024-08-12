// SPDX-License-Identifier: MIT
#pragma once
#include <ddc/ddc.hpp>

#include "Lagrange.hpp"
#include "ddc_alias_inline_functions.hpp"
#include "ddc_aliases.hpp"
#include "iinterpolator.hpp"

/**
 * @brief A class for interpolating a function using Lagrange polynomials.
 * It is designed to work with both uniform and non-uniform mesh, and have the advantage to be local.
 *
 */
template <class GridInterp, BCond BcMin, BCond BcMax, class... Grid1D>
class LagrangeInterpolator : public IInterpolator<GridInterp, Grid1D...>
{
    using InterpDim = typename GridInterp::continuous_dimension_type;
    using deriv_type = typename IInterpolator<GridInterp, Grid1D...>::deriv_type;
    using batched_derivs_idx_range_type =
            typename IInterpolator<GridInterp, Grid1D...>::batched_derivs_idx_range_type;

private:
    int m_degree;
    IdxStep<GridInterp> m_ghost;

public:
    /**
     * @brief Create a  Lagrange interpolator object.
     * @param[in] degree Degree of polynomials
     * @param[in] ghost  Discrete vector which gives the number of ghost points. By default choose 2.
    */
    LagrangeInterpolator(int degree, IdxStep<GridInterp> ghost) : m_degree(degree), m_ghost(ghost)
    {
    }

    ~LagrangeInterpolator() override = default;

    batched_derivs_idx_range_type batched_derivs_idx_range_xmin(
            IdxRange<Grid1D...> dom) const override
    {
        return ddc::replace_dim_of<
                GridInterp,
                deriv_type>(dom, IdxRange<deriv_type>(Idx<deriv_type>(1), IdxStep<deriv_type>(0)));
    }

    batched_derivs_idx_range_type batched_derivs_idx_range_xmax(
            IdxRange<Grid1D...> dom) const override
    {
        return ddc::replace_dim_of<
                GridInterp,
                deriv_type>(dom, IdxRange<deriv_type>(Idx<deriv_type>(1), IdxStep<deriv_type>(0)));
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
    Field<double, IdxRange<Grid1D...>> operator()(
            Field<double, IdxRange<Grid1D...>> inout_data,
            Field<const Coord<typename GridInterp::continuous_dimension_type>, IdxRange<Grid1D...>>
                    coordinates,
            [[maybe_unused]] std::optional<Field<
                    double const,
                    typename IInterpolator<GridInterp, Grid1D...>::batched_derivs_idx_range_type>>
                    derivs_xmin
            = std::nullopt,
            [[maybe_unused]] std::optional<Field<
                    double const,
                    typename IInterpolator<GridInterp, Grid1D...>::batched_derivs_idx_range_type>>
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
        auto inout_data_tmp = get_field(inout_data_tmp_alloc);
        auto batch_idx_range = ddc::remove_dims_of<GridInterp>(get_idx_range(inout_data));
        ddc::parallel_for_each(
                Kokkos::DefaultExecutionSpace(),
                batch_idx_range,
                KOKKOS_LAMBDA(typename decltype(batch_idx_range)::discrete_element_type const i) {
                    Lagrange<Kokkos::DefaultExecutionSpace, GridInterp, BcMin, BcMax> evaluator(
                            deg,
                            inout_data_tmp[i],
                            get_idx_range<GridInterp>(inout_data),
                            ghost);
                    for (Idx<GridInterp> j : get_idx_range<GridInterp>(inout_data)) {
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
template <class GridInterp, BCond BcMin, BCond BcMax, class... Grid1D>
class PreallocatableLagrangeInterpolator : public IPreallocatableInterpolator<GridInterp, Grid1D...>
{
    LagrangeInterpolator<GridInterp, BcMin, BcMax, Grid1D...> const& m_evaluator;

public:
    /**
     * @brief Create an object capable of creating LagrangeInterpolator objects.
     * @param[in] evaluator An operator which evaluates the value of the interpolation polynomial at requested coordinates.
     */
    explicit PreallocatableLagrangeInterpolator(
            LagrangeInterpolator<GridInterp, BcMin, BcMax, Grid1D...> const& evaluator)
        : m_evaluator(evaluator)
    {
    }

    ~PreallocatableLagrangeInterpolator() override = default;

    /**
     * Create an instance of the LagrangeInterpolator class.
     *
     * @return A unique pointer to an instance of the LagrangeInterpolator class.
     */
    std::unique_ptr<IInterpolator<GridInterp, Grid1D...>> preallocate() const override
    {
        return std::make_unique<LagrangeInterpolator<GridInterp, BcMin, BcMax, Grid1D...>>(
                m_evaluator);
    }
};
