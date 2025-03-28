

# File Lagrange\_interpolator.hpp

[**File List**](files.md) **>** [**interpolation**](dir_264890e5c091f8c8d7fe1f842870c25e.md) **>** [**Lagrange\_interpolator.hpp**](Lagrange__interpolator_8hpp.md)

[Go to the documentation of this file](Lagrange__interpolator_8hpp.md)


```C++
// SPDX-License-Identifier: MIT
#pragma once
#include <ddc/ddc.hpp>

#include "Lagrange.hpp"
#include "ddc_alias_inline_functions.hpp"
#include "ddc_aliases.hpp"
#include "iinterpolator.hpp"

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
    LagrangeInterpolator(int degree, IdxStep<GridInterp> ghost) : m_degree(degree), m_ghost(ghost)
    {
    }

    ~LagrangeInterpolator() override = default;

    batched_derivs_idx_range_type batched_derivs_idx_range_xmin(
            IdxRange<Grid1D...> idx_range) const override
    {
        return ddc::replace_dim_of<GridInterp, deriv_type>(
                idx_range,
                IdxRange<deriv_type>(Idx<deriv_type>(1), IdxStep<deriv_type>(0)));
    }

    batched_derivs_idx_range_type batched_derivs_idx_range_xmax(
            IdxRange<Grid1D...> idx_range) const override
    {
        return ddc::replace_dim_of<GridInterp, deriv_type>(
                idx_range,
                IdxRange<deriv_type>(Idx<deriv_type>(1), IdxStep<deriv_type>(0)));
    }

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
        auto const interp_range = get_idx_range<GridInterp>(inout_data);
        ddc::parallel_for_each(
                Kokkos::DefaultExecutionSpace(),
                batch_idx_range,
                KOKKOS_LAMBDA(typename decltype(batch_idx_range)::discrete_element_type const i) {
                    Lagrange<Kokkos::DefaultExecutionSpace, GridInterp, BcMin, BcMax>
                            evaluator(deg, inout_data_tmp[i], interp_range, ghost);
                    for (Idx<GridInterp> j : interp_range) {
                        inout_data(i, j) = evaluator.evaluate(coordinates(i, j));
                    }
                });
        return inout_data;
    }
};

template <class GridInterp, BCond BcMin, BCond BcMax, class... Grid1D>
class PreallocatableLagrangeInterpolator : public IPreallocatableInterpolator<GridInterp, Grid1D...>
{
    LagrangeInterpolator<GridInterp, BcMin, BcMax, Grid1D...> const& m_evaluator;

public:
    explicit PreallocatableLagrangeInterpolator(
            LagrangeInterpolator<GridInterp, BcMin, BcMax, Grid1D...> const& evaluator)
        : m_evaluator(evaluator)
    {
    }

    ~PreallocatableLagrangeInterpolator() override = default;

    std::unique_ptr<IInterpolator<GridInterp, Grid1D...>> preallocate() const override
    {
        return std::make_unique<LagrangeInterpolator<GridInterp, BcMin, BcMax, Grid1D...>>(
                m_evaluator);
    }
};
```


