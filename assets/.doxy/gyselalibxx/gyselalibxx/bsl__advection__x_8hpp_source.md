

# File bsl\_advection\_x.hpp

[**File List**](files.md) **>** [**advection**](dir_b90fde0f10c67a9aef841a6e6700f1f6.md) **>** [**bsl\_advection\_x.hpp**](bsl__advection__x_8hpp.md)

[Go to the documentation of this file](bsl__advection__x_8hpp.md)


```C++
// SPDX-License-Identifier: MIT
#pragma once
#include <type_traits>

#include <ddc/ddc.hpp>

#include "ddc_alias_inline_functions.hpp"
#include "ddc_aliases.hpp"
#include "i_interpolation.hpp"
#include "iadvectionx.hpp"
#include "species_info.hpp"

template <class Geometry, concepts::Interpolation FunctionInterpolator, class DataType = double>
class BslAdvectionSpatial
    : public IAdvectionSpatial<
              Geometry,
              typename InterpolationBuilderTraits<
                      typename FunctionInterpolator::BuilderType>::interpolation_grid_type,
              DataType>
{
    static_assert(std::is_floating_point_v<DataType>);

    using FunctionBuilder = typename FunctionInterpolator::BuilderType;
    using FunctionEvaluator = typename FunctionInterpolator::EvaluatorType;

    using GridX = typename InterpolationBuilderTraits<FunctionBuilder>::interpolation_grid_type;

    using GridV = typename Geometry::template velocity_dim_for<GridX>;
    using IdxRangeFdistrib = typename Geometry::IdxRangeFdistribu;
    using IdxX = Idx<GridX>;
    using IdxV = Idx<GridV>;
    using DimX = typename GridX::continuous_dimension_type;
    using DimV = typename GridV::continuous_dimension_type;
    using IdxRangeSpaceVelocity
            = ddc::remove_dims_of_t<typename Geometry::IdxRangeFdistribu, Species>;

private:
    using IdxRangeFunctionBasis = typename InterpolationBuilderTraits<
            FunctionBuilder>::template batched_basis_idx_range_type<IdxRangeSpaceVelocity>;

    FunctionBuilder const& m_function_builder;
    FunctionEvaluator const& m_function_evaluator;

public:
    [[deprecated]] explicit BslAdvectionSpatial(
            FunctionBuilder const& function_builder,
            FunctionEvaluator const& function_evaluator)
        : m_function_builder(function_builder)
        , m_function_evaluator(function_evaluator)
    {
    }

    explicit BslAdvectionSpatial(FunctionInterpolator const& function_interpolator)
        : m_function_builder(function_interpolator.get_builder())
        , m_function_evaluator(function_interpolator.get_evaluator())
    {
    }

    ~BslAdvectionSpatial() override = default;

    Field<DataType, IdxRangeFdistrib> operator()(
            Field<DataType, IdxRangeFdistrib> const allfdistribu,
            DataType const dt) const override
    {
        using IdxRangeBatch = ddc::remove_dims_of_t<IdxRangeFdistrib, Species, GridX>;
        using IdxBatch = typename IdxRangeBatch::discrete_element_type;

        Kokkos::Profiling::pushRegion("BslAdvectionSpatial");
        IdxRangeFdistrib const idx_range = get_idx_range(allfdistribu);
        IdxRange<GridX> const x_idx_range = ddc::select<GridX>(idx_range);
        IdxRange<Species> const sp_idx_range = ddc::select<Species>(idx_range);

        // pre-allocate some memory to prevent allocation later in loop
        IdxRangeSpaceVelocity batched_feet_idx_range(idx_range);
        FieldMem<Coord<DimX>, IdxRangeSpaceVelocity> feet_coords_alloc(
                "feet_coords (BslAdvectionVelocity::operator())",
                batched_feet_idx_range);
        Field<Coord<DimX>, IdxRangeSpaceVelocity> feet_coords(get_field(feet_coords_alloc));
        DFieldMem<IdxRangeFunctionBasis> function_coefs_alloc(
                "function_coefs (BslAdvectionVelocity::operator())",
                batched_basis_idx_range(m_function_builder, batched_feet_idx_range));

        IdxRangeBatch batch_idx_range(idx_range);

        for (IdxSp const isp : sp_idx_range) {
            DataType const sqrt_me_on_mspecies = std::sqrt(mass(ielec()) / mass(isp));
            const std::source_location location = std::source_location::current();
            ddc::parallel_for_each(
                    location.function_name(),
                    Kokkos::DefaultExecutionSpace(),
                    batch_idx_range,
                    KOKKOS_LAMBDA(IdxBatch const ib) {
                        // compute the displacement
                        IdxV const iv(ib);
                        Coord<DimV> const coord_iv = ddc::coordinate(iv);
                        DataType const dx = sqrt_me_on_mspecies * dt * coord_iv;

                        // compute the coordinates of the feet
                        for (IdxX const ix : x_idx_range) {
                            feet_coords(ix, ib) = Coord<DimX>(ddc::coordinate(ix) - dx);
                        }
                    });
            m_function_builder(get_field(function_coefs_alloc), get_const_field(allfdistribu[isp]));
            m_function_evaluator(
                    allfdistribu[isp],
                    get_const_field(feet_coords),
                    get_const_field(function_coefs_alloc));
        }

        Kokkos::Profiling::popRegion();
        return allfdistribu;
    }
};
```


