// SPDX-License-Identifier: MIT
#pragma once
#include <type_traits>

#include <ddc/ddc.hpp>

#include "ddc_alias_inline_functions.hpp"
#include "ddc_aliases.hpp"
#include "ddc_helper.hpp"
#include "i_interpolation_builder.hpp"
#include "i_interpolation_evaluator.hpp"
#include "iadvectionvx.hpp"
#include "species_info.hpp"

/**
 * @brief A class which computes the velocity advection along the dimension of interest GridV. Working for every Cartesian geometry.
 */
template <
        class Geometry,
        class GridV,
        concepts::InterpolationBuilder FunctionBuilder,
        concepts::InterpolationEvaluator FunctionEvaluator,
        class DataType = double>
class BslAdvectionVelocity : public IAdvectionVelocity<Geometry, GridV, DataType>
{
    static_assert(std::is_floating_point_v<DataType>);

    using IdxRangeFdistribu = typename Geometry::IdxRangeFdistribu;
    using IdxRangeSpatial = typename Geometry::IdxRangeSpatial;
    using IdxSpatial = typename IdxRangeSpatial::discrete_element_type;
    using IdxV = Idx<GridV>;
    using DimV = typename GridV::continuous_dimension_type;
    using IdxRangeSpaceVelocity
            = ddc::remove_dims_of_t<typename Geometry::IdxRangeFdistribu, Species>;

private:
    using IdxRangeFunctionDeriv = typename InterpolationBuilderTraits<
            FunctionBuilder>::template batched_derivs_idx_range_type<IdxRangeSpaceVelocity>;
    using IdxRangeFunctionBasis = typename InterpolationBuilderTraits<
            FunctionBuilder>::template batched_basis_idx_range_type<IdxRangeSpaceVelocity>;

    FunctionBuilder const& m_function_builder;
    FunctionEvaluator const& m_function_evaluator;

public:
    /**
     * @brief Constructor
     * @param[in] function_builder Builder along the GridV direction used to build
     *          the interpolation representation of the advected function.
     * @param[in] function_evaluator Evaluator along the GridV direction used to evaluate
     *          the advected function at the characteristic feet.
     */
    explicit BslAdvectionVelocity(
            FunctionBuilder const& function_builder,
            FunctionEvaluator const& function_evaluator)
        : m_function_builder(function_builder)
        , m_function_evaluator(function_evaluator)
    {
    }

    ~BslAdvectionVelocity() override = default;

    /**
     * @brief Advects fdistribu along GridV for a duration dt.
     * @param[in, out] allfdistribu Reference to the whole distribution function for one species, allocated on the device (ie it lets the choice of the location depend on the build configuration).
     * @param[in] electric_field Reference to the electric field which derives from electrostatic potential, allocated on the device.
     * @param[in] dt Time step
     * @return A reference to the allfdistribu array containing the value of the function at the coordinates.
     */
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

        FieldMem<DataType, IdxRangeFunctionDeriv> derivs_min(
                m_function_builder.batched_derivs_xmin_domain(batched_feet_idx_range));
        FieldMem<DataType, IdxRangeFunctionDeriv> derivs_max(
                m_function_builder.batched_derivs_xmax_domain(batched_feet_idx_range));
        ddc::parallel_fill(derivs_min, 0.);
        ddc::parallel_fill(derivs_max, 0.);

        // pre-allocate some memory to prevent allocation later in loop
        FieldMem<Coord<DimV>, IdxRangeSpaceVelocity> feet_coords_alloc(batched_feet_idx_range);
        Field<Coord<DimV>, IdxRangeSpaceVelocity> feet_coords(get_field(feet_coords_alloc));
        DFieldMem<IdxRangeFunctionBasis> function_coefs_alloc(
                batched_basis_idx_range(m_function_builder, batched_feet_idx_range));

        IdxRangeBatch batch_idx_range(idx_range);

        ddc::host_for_each(idx_range_sp, [&](IdxSp const isp) {
            DataType const charge_proxy
                    = charge(isp); // TODO: consider proper way to access charge from device
            DataType const sqrt_me_on_mspecies = std::sqrt(mass(ielec()) / mass(isp));
            ddc::parallel_for_each(
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
            m_function_builder(
                    get_field(function_coefs_alloc),
                    get_const_field(allfdistribu[isp]),
                    std::optional(get_const_field(derivs_min)),
                    std::optional(get_const_field(derivs_max)));
            m_function_evaluator(
                    allfdistribu[isp],
                    get_const_field(feet_coords),
                    get_const_field(function_coefs_alloc));
        });

        Kokkos::Profiling::popRegion();
        return allfdistribu;
    }
};
