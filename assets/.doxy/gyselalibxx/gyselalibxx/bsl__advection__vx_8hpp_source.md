

# File bsl\_advection\_vx.hpp

[**File List**](files.md) **>** [**advection**](dir_b90fde0f10c67a9aef841a6e6700f1f6.md) **>** [**bsl\_advection\_vx.hpp**](bsl__advection__vx_8hpp.md)

[Go to the documentation of this file](bsl__advection__vx_8hpp.md)


```C++
// SPDX-License-Identifier: MIT
#pragma once
#include <type_traits>

#include <ddc/ddc.hpp>

#include "ddc_alias_inline_functions.hpp"
#include "ddc_aliases.hpp"
#include "ddc_helper.hpp"
#include "i_interpolation.hpp"
#include "iadvectionvx.hpp"
#include "species_info.hpp"

template <class Geometry, concepts::Interpolation1D FunctionInterpolator, class DataType = double>
class BslAdvectionVelocity
    : public IAdvectionVelocity<
              Geometry,
              interpolation_grid_t<typename FunctionInterpolator::BuilderType>,
              DataType>
{
    static_assert(std::is_floating_point_v<DataType>);

    using FunctionBuilder = typename FunctionInterpolator::BuilderType;
    using FunctionEvaluator = typename FunctionInterpolator::EvaluatorType;

    using GridV = interpolation_grid_t<typename FunctionInterpolator::BuilderType>;

    using IdxRangeFdistribu = typename Geometry::IdxRangeFdistribu;
    using IdxRangeSpatial = typename Geometry::IdxRangeSpatial;
    using IdxSpatial = typename IdxRangeSpatial::discrete_element_type;
    using IdxV = Idx<GridV>;
    using DimV = typename GridV::continuous_dimension_type;
    using IdxRangeSpaceVelocity
            = ddc::remove_dims_of_t<typename Geometry::IdxRangeFdistribu, Species>;

private:
    using IdxRangeFunctionBasis = typename InterpolationBuilderTraits<
            FunctionBuilder>::template batched_basis_idx_range_type<IdxRangeSpaceVelocity>;

    FunctionBuilder const& m_function_builder;
    FunctionEvaluator const& m_function_evaluator;

public:
    [[deprecated]] explicit BslAdvectionVelocity(
            FunctionBuilder const& function_builder,
            FunctionEvaluator const& function_evaluator)
        : m_function_builder(function_builder)
        , m_function_evaluator(function_evaluator)
    {
    }

    explicit BslAdvectionVelocity(FunctionInterpolator const& function_interpolator)
        : m_function_builder(function_interpolator.get_builder())
        , m_function_evaluator(function_interpolator.get_evaluator())
    {
    }

    ~BslAdvectionVelocity() override = default;

    Field<DataType, IdxRangeFdistribu> operator()(
            Field<DataType, IdxRangeFdistribu> const allfdistribu,
            ConstField<DataType, IdxRangeSpatial> const electric_field,
            DataType const dt) const override
    {
        using IdxRangeBatch = ddc::remove_dims_of_t<IdxRangeFdistribu, Species, GridV>;
        using IdxBatch = typename IdxRangeBatch::discrete_element_type;


        Kokkos::Profiling::pushRegion("BslAdvectionVelocity");
        IdxRangeFdistribu const idx_range = get_idx_range(allfdistribu);
        IdxRange<GridV> const idx_range_v = ddc::select<GridV>(idx_range);
        IdxRange<Species> const idx_range_sp = ddc::select<Species>(idx_range);

        IdxRangeSpaceVelocity batched_feet_idx_range(idx_range);

        // pre-allocate some memory to prevent allocation later in loop
        FieldMem<Coord<DimV>, IdxRangeSpaceVelocity> feet_coords_alloc(
                "feet_coords (BslAdvectionVelocity::operator())",
                batched_feet_idx_range);
        Field<Coord<DimV>, IdxRangeSpaceVelocity> feet_coords(get_field(feet_coords_alloc));
        DFieldMem<IdxRangeFunctionBasis> function_coefs_alloc(
                "function_coefs (BslAdvectionVelocity::operator())",
                batched_basis_idx_range(m_function_builder, batched_feet_idx_range));

        IdxRangeBatch batch_idx_range(idx_range);

        ddc::host_for_each(idx_range_sp, [&](IdxSp const isp) {
            DataType const charge_proxy
                    = charge(isp); // TODO: consider proper way to access charge from device
            DataType const sqrt_me_on_mspecies = std::sqrt(mass(ielec()) / mass(isp));
            const std::source_location location = std::source_location::current();
            ddc::parallel_for_each(
                    location.function_name(),
                    Kokkos::DefaultExecutionSpace(),
                    batch_idx_range,
                    KOKKOS_LAMBDA(IdxBatch const ib) {
                        IdxSpatial const ix(ib);
                        // compute the displacement
                        DataType const dvx
                                = charge_proxy * sqrt_me_on_mspecies * dt * electric_field(ix);

                        // compute the coordinates of the feet
                        for (IdxV const iv : idx_range_v) {
                            feet_coords(iv, ib) = Coord<DimV>(ddc::coordinate(iv) - dvx);
                        }
                    });
            m_function_builder(get_field(function_coefs_alloc), get_const_field(allfdistribu[isp]));
            m_function_evaluator(
                    allfdistribu[isp],
                    get_const_field(feet_coords),
                    get_const_field(function_coefs_alloc));
        });

        Kokkos::Profiling::popRegion();
        return allfdistribu;
    }
};
```


