// SPDX-License-Identifier: MIT
#pragma once
#include <ddc/ddc.hpp>

#include "ddc_alias_inline_functions.hpp"
#include "ddc_aliases.hpp"
#include "i_interpolation_builder.hpp"
#include "i_interpolation_evaluator.hpp"
#include "iadvectionx.hpp"
#include "species_info.hpp"

/**
 * @brief A class which computes the spatial advection along the dimension of interest GridX. Working for every Cartesian geometry.
 */
template <class Geometry, class GridX, class FunctionBuilder, class FunctionEvaluator>
class BslAdvectionSpatial : public IAdvectionSpatial<Geometry, GridX>
{
    static_assert(concepts::InterpolationBuilder<FunctionBuilder>);
    static_assert(concepts::InterpolationEvaluator<FunctionEvaluator>);

    using GridV = typename Geometry::template velocity_dim_for<GridX>;
    using IdxRangeFdistrib = typename Geometry::IdxRangeFdistribu;
    using IdxX = Idx<GridX>;
    using IdxV = Idx<GridV>;
    using DimX = typename GridX::continuous_dimension_type;
    using DimV = typename GridV::continuous_dimension_type;
    using IdxRangeSpaceVelocity
            = ddc::remove_dims_of_t<typename Geometry::IdxRangeFdistribu, Species>;

private:
    using FunctionBasisIdxRange = typename InterpolationBuilderTraits<
            FunctionBuilder>::template batched_basis_idx_range_type<IdxRangeSpaceVelocity>;

    FunctionBuilder const& m_function_builder;
    FunctionEvaluator const& m_function_evaluator;

public:
    /**
     * @brief Constructor
     * @param[in] function_builder Builder along the GridX direction used to build
     *          the interpolation representation of the advected function.
     * @param[in] function_evaluator Evaluator along the GridX direction used to evaluate
     *          the advected function at the characteristic feet.
     */
    explicit BslAdvectionSpatial(
            FunctionBuilder const& function_builder,
            FunctionEvaluator const& function_evaluator)
        : m_function_builder(function_builder)
        , m_function_evaluator(function_evaluator)
    {
    }

    ~BslAdvectionSpatial() override = default;

    /**
     * @brief Advects fdistribu along GridX for a duration dt.
     * @param[in, out] allfdistribu Reference to the whole distribution function for one species, allocated on the device (ie it lets the choice of the location depend on the build configuration).
     * @param[in] dt Time step
     * @return A reference to the allfdistribu array containing the value of the function at the coordinates.
     */
    Field<double, IdxRangeFdistrib> operator()(
            Field<double, IdxRangeFdistrib> const allfdistribu,
            double const dt) const override
    {
        using IdxRangeBatch = ddc::remove_dims_of_t<IdxRangeFdistrib, Species, GridX>;
        using IdxBatch = typename IdxRangeBatch::discrete_element_type;

        Kokkos::Profiling::pushRegion("BslAdvectionSpatial");
        IdxRangeFdistrib const idx_range = get_idx_range(allfdistribu);
        IdxRange<GridX> const x_idx_range = ddc::select<GridX>(idx_range);
        IdxRange<Species> const sp_idx_range = ddc::select<Species>(idx_range);

        // pre-allocate some memory to prevent allocation later in loop
        IdxRangeSpaceVelocity batched_feet_idx_range(idx_range);
        FieldMem<Coord<DimX>, IdxRangeSpaceVelocity> feet_coords_alloc(batched_feet_idx_range);
        Field<Coord<DimX>, IdxRangeSpaceVelocity> feet_coords(get_field(feet_coords_alloc));
        DFieldMem<FunctionBasisIdxRange> function_coefs_alloc(
                batched_basis_idx_range(m_function_builder, batched_feet_idx_range));

        IdxRangeBatch batch_idx_range(idx_range);

        for (IdxSp const isp : sp_idx_range) {
            double const sqrt_me_on_mspecies = std::sqrt(mass(ielec()) / mass(isp));
            ddc::parallel_for_each(
                    Kokkos::DefaultExecutionSpace(),
                    batch_idx_range,
                    KOKKOS_LAMBDA(IdxBatch const ib) {
                        // compute the displacement
                        IdxV const iv(ib);
                        Coord<DimV> const coord_iv = ddc::coordinate(iv);
                        double const dx = sqrt_me_on_mspecies * dt * coord_iv;

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
