

# File bsl\_advection\_polar.hpp

[**File List**](files.md) **>** [**advection**](dir_b90fde0f10c67a9aef841a6e6700f1f6.md) **>** [**bsl\_advection\_polar.hpp**](bsl__advection__polar_8hpp.md)

[Go to the documentation of this file](bsl__advection__polar_8hpp.md)


```C++
// SPDX-License-Identifier: MIT
#pragma once

#include "ddc_alias_inline_functions.hpp"
#include "ddc_aliases.hpp"
#include "i_interpolation.hpp"
#include "i_interpolation_builder.hpp"
#include "indexed_tensor.hpp"
#include "l_norm_tools.hpp"
#include "metric_tensor_evaluator.hpp"
#include "vector_field.hpp"
#include "vector_field_mem.hpp"



template <class FootFinder, class LogicalToPhysicalMapping, concepts::Interpolation Interpolator2D>
class BslAdvectionPolar
{
    using R = typename CoordWithOPoint<
            typename LogicalToPhysicalMapping::CoordArg>::curvilinear_tag_r;
    using Theta = typename CoordWithOPoint<
            typename LogicalToPhysicalMapping::CoordArg>::curvilinear_tag_theta;

    using CoordRTheta = typename LogicalToPhysicalMapping::CoordArg;
    using CoordXY = typename LogicalToPhysicalMapping::CoordResult;

    using CartesianBasis = ddc::to_type_seq_t<CoordXY>;
    using CurvilinearBasis = ddc::to_type_seq_t<CoordRTheta>;

    using IdxRangeBatched = typename FootFinder::IdxRangeOperator;

    using GridR = find_grid_t<R, ddc::to_type_seq_t<IdxRangeBatched>>;
    using GridTheta = find_grid_t<Theta, ddc::to_type_seq_t<IdxRangeBatched>>;

    using IdxRangeRTheta = IdxRange<GridR, GridTheta>;
    using IdxRangeR = IdxRange<GridR>;
    using IdxRangeTheta = IdxRange<GridTheta>;
    using IdxRangeBatch = ddc::remove_dims_of_t<IdxRangeBatched, GridR, GridTheta>;

    using IdxRTheta = typename IdxRangeRTheta::discrete_element_type;
    using IdxR = Idx<GridR>;
    using IdxTheta = Idx<GridTheta>;
    using IdxBatched = typename IdxRangeBatched::discrete_element_type;
    using IdxBatch = typename IdxRangeBatch::discrete_element_type;
    using IdxStepR = typename IdxRangeR::discrete_vector_type;

    using MemorySpace = typename FootFinder::memory_space;
    using ExecSpace = typename FootFinder::ExecSpace;

    using Builder2D = typename Interpolator2D::BuilderType;
    using Evaluator2D = typename Interpolator2D::EvaluatorType;

    using IdxRangeCoeffBatchedRTheta = typename InterpolationBuilderTraits<
            Builder2D>::template batched_basis_idx_range_type<IdxRangeBatched>;

    using DFieldFDistribu = DField<IdxRangeBatched, MemorySpace>;

    using CFieldMemFeetRTheta = FieldMem<CoordRTheta, IdxRangeBatched, MemorySpace>;

    using CFieldFeetRTheta = Field<CoordRTheta, IdxRangeBatched, MemorySpace>;

    using DVectorFieldMemAdvectionXY
            = DVectorFieldMem<IdxRangeBatched, CartesianBasis, MemorySpace>;
    using DVectorFieldAdvectionXY = DVectorField<IdxRangeBatched, CartesianBasis, MemorySpace>;
    using DVectorConstFieldAdvectionXY
            = DVectorConstField<IdxRangeBatched, CartesianBasis, MemorySpace>;

    using DVectorConstFieldAdvection = DVectorConstField<
            IdxRangeBatched,
            VectorIndexSet<typename FootFinder::AdvDim1, typename FootFinder::AdvDim2>,
            MemorySpace>;

    using DVectorFieldMemAdvectionXYOnBatch
            = DVectorFieldMem<IdxRangeBatch, CartesianBasis, MemorySpace>;
    using DVectorFieldAdvectionXYOnBatch = DVectorField<IdxRangeBatch, CartesianBasis, MemorySpace>;

    using DVectorFieldAdvectionRTheta
            = DVectorField<IdxRangeBatched, CurvilinearBasis, MemorySpace>;
    using DVectorConstFieldAdvectionRTheta
            = DVectorConstField<IdxRangeBatched, CurvilinearBasis, MemorySpace>;

private:
    Builder2D const& m_builder_2d;

    Evaluator2D const& m_evaluator_2d;

    FootFinder& m_find_feet_method;

    LogicalToPhysicalMapping const& m_logical_to_physical_mapping;

    std::optional<IdxRangeBatched> m_idx_range_advected_points;

public:
    BslAdvectionPolar(
            Interpolator2D const& interpolator_2d,
            FootFinder& foot_finder,
            LogicalToPhysicalMapping const& logical_to_physical_mapping,
            std::optional<IdxRangeBatched> idx_range_advected_points = std::nullopt)
        : m_builder_2d(interpolator_2d.get_builder())
        , m_evaluator_2d(interpolator_2d.get_evaluator())
        , m_find_feet_method(foot_finder)
        , m_logical_to_physical_mapping(logical_to_physical_mapping)
        , m_idx_range_advected_points(idx_range_advected_points)
    {
    }

    ~BslAdvectionPolar() = default;


    DFieldFDistribu operator()(
            DFieldFDistribu allfdistribu,
            DVectorConstFieldAdvection advection_field,
            double dt) const
    {
        // Pre-allocate coefficient storage
        DFieldMem<IdxRangeCoeffBatchedRTheta, MemorySpace> coefs_alloc(
                batched_basis_idx_range(m_builder_2d, get_idx_range(allfdistribu)));

        // Compute the feet of the characteristics at tn -----------------------------------------
        typename FootFinder::ElementwiseOperator find_foot_alloc
                = m_find_feet_method(get_const_field(advection_field));
        typename FootFinder::ElementwiseOperator::GPUCompat find_foot = find_foot_alloc(dt);

        // Interpolate the function on the feet of the characteristics. --------------------------
        Kokkos::Profiling::pushRegion("(GSLX) BslAdvectionPolar/Interpolation");
        m_builder_2d(get_field(coefs_alloc), get_const_field(allfdistribu));

        Evaluator2D const& evaluator_2d_proxy = m_evaluator_2d;

        DConstField<IdxRangeCoeffBatchedRTheta, MemorySpace> coefs = get_const_field(coefs_alloc);

        IdxRangeBatched idx_range_advected_points;
        if (m_idx_range_advected_points) {
            idx_range_advected_points = *m_idx_range_advected_points;
        } else {
            idx_range_advected_points = get_idx_range(allfdistribu);
        }

        const std::source_location location = std::source_location::current();
        ddc::parallel_for_each(
                location.function_name(),
                ExecSpace(),
                idx_range_advected_points,
                KOKKOS_LAMBDA(IdxBatched const idx) {
                    IdxBatch batch_idx(idx);
                    allfdistribu(idx) = evaluator_2d_proxy(find_foot(idx), coefs[batch_idx]);
                });
        Kokkos::Profiling::popRegion();

        unify_value_at_centre_pt(allfdistribu);

        return allfdistribu;
    }


    template <class T>
    static void unify_value_at_centre_pt(Field<T, IdxRangeBatched, MemorySpace> values)
    {
        IdxRangeBatched full_idx_range = get_idx_range(values);
        IdxRangeBatch const batched_idx_range(full_idx_range);
        IdxRangeR const r_idx_range(full_idx_range);
        IdxRangeTheta const theta_idx_range(full_idx_range);
        IdxR r0_idx = r_idx_range.front();
        IdxTheta theta0_idx = theta_idx_range.front();
        if (std::fabs(ddc::coordinate(r0_idx)) < 1e-15) {
            const std::source_location location = std::source_location::current();
            ddc::parallel_for_each(
                    location.function_name(),
                    ExecSpace(),
                    batched_idx_range,
                    KOKKOS_LAMBDA(const IdxBatch ib) {
                        for (IdxTheta itheta : theta_idx_range) {
                            values(ib, r0_idx, itheta) = values(ib, r0_idx, theta0_idx);
                        }
                    });
        }
    }
};
```


