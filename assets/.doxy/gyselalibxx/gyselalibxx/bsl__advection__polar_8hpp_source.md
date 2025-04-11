

# File bsl\_advection\_polar.hpp

[**File List**](files.md) **>** [**advection**](dir_b90fde0f10c67a9aef841a6e6700f1f6.md) **>** [**bsl\_advection\_polar.hpp**](bsl__advection__polar_8hpp.md)

[Go to the documentation of this file](bsl__advection__polar_8hpp.md)


```C++
// SPDX-License-Identifier: MIT
#pragma once

#include "ddc_alias_inline_functions.hpp"
#include "ddc_aliases.hpp"
#include "i_interpolator_2d.hpp"
#include "indexed_tensor.hpp"
#include "metric_tensor_evaluator.hpp"
#include "spline_interpolator_2d.hpp"
#include "spline_polar_foot_finder.hpp"
#include "vector_field.hpp"
#include "vector_field_mem.hpp"



template <class FootFinder, class LogicalToPhysicalMapping, class InterpolatorPolar>
class BslAdvectionPolar
{
    using R = typename LogicalToPhysicalMapping::curvilinear_tag_r;
    using Theta = typename LogicalToPhysicalMapping::curvilinear_tag_theta;

    using CoordRTheta = typename LogicalToPhysicalMapping::CoordArg;
    using CoordXY = typename LogicalToPhysicalMapping::CoordResult;

    using CartesianBasis = ddc::to_type_seq_t<CoordXY>;
    using CurvilinearBasis = ddc::to_type_seq_t<CoordRTheta>;

    using IdxRangeBatched = typename FootFinder::IdxRangeOperator;

    using GridR = find_grid_t<R, ddc::to_type_seq_t<IdxRangeBatched>>;
    using GridTheta = find_grid_t<Theta, ddc::to_type_seq_t<IdxRangeBatched>>;

    using IdxRangeRTheta = IdxRange<GridR, GridTheta>;
    using IdxRangeTheta = IdxRange<GridTheta>;

    using IdxRTheta = typename IdxRangeRTheta::discrete_element_type;
    using IdxStepRTheta = typename IdxRangeRTheta::discrete_vector_type;

    using MemorySpace = typename FootFinder::memory_space;

    using DFieldFDistribu = DField<IdxRangeBatched, MemorySpace>;

    using CFieldMemFeetRTheta = FieldMem<CoordRTheta, IdxRangeBatched, MemorySpace>;

    using CFieldFeetRTheta = Field<CoordRTheta, IdxRangeBatched, MemorySpace>;

    using DVectorFieldMemAdvectionXY
            = DVectorFieldMem<IdxRangeBatched, CartesianBasis, MemorySpace>;
    using DVectorFieldAdvectionXY = DVectorField<IdxRangeBatched, CartesianBasis, MemorySpace>;
    using DVectorConstFieldAdvectionXY
            = DVectorConstField<IdxRangeBatched, CartesianBasis, MemorySpace>;

    using DVectorFieldAdvectionRTheta
            = DVectorField<IdxRangeBatched, CurvilinearBasis, MemorySpace>;
    using DVectorConstFieldAdvectionRTheta
            = DVectorConstField<IdxRangeBatched, CurvilinearBasis, MemorySpace>;

private:
    InterpolatorPolar const& m_interpolator;

    FootFinder const& m_find_feet;

    LogicalToPhysicalMapping const& m_logical_to_physical_mapping;


public:
    BslAdvectionPolar(
            InterpolatorPolar const& function_interpolator,
            FootFinder const& foot_finder,
            LogicalToPhysicalMapping const& logical_to_physical_mapping)
        : m_interpolator(function_interpolator)
        , m_find_feet(foot_finder)
        , m_logical_to_physical_mapping(logical_to_physical_mapping)
    {
    }

    ~BslAdvectionPolar() = default;


    DFieldFDistribu operator()(
            DFieldFDistribu allfdistribu,
            DVectorConstFieldAdvectionXY advection_field_xy,
            double dt) const
    {
        // Pre-allocate some memory to prevent allocation later in loop
        std::unique_ptr<IInterpolator2D<IdxRangeRTheta, IdxRangeRTheta>> const interpolator_ptr
                = m_interpolator.preallocate();

        // Initialise the feet
        CFieldMemFeetRTheta feet_rtheta_alloc(get_idx_range(advection_field_xy));
        CFieldFeetRTheta feet_rtheta = get_field(feet_rtheta_alloc);
        ddc::parallel_for_each(
                Kokkos::DefaultExecutionSpace(),
                get_idx_range(advection_field_xy),
                KOKKOS_LAMBDA(IdxRTheta const irtheta) {
                    feet_rtheta(irtheta) = ddc::coordinate(irtheta);
                });

        // Compute the feet of the characteristics at tn -----------------------------------------
        m_find_feet(feet_rtheta, get_const_field(advection_field_xy), dt);

        // Interpolate the function on the feet of the characteristics. --------------------------
        (*interpolator_ptr)(get_field(allfdistribu), get_const_field(feet_rtheta));

        return allfdistribu;
    }


    DFieldFDistribu operator()(
            DFieldFDistribu allfdistribu,
            DVectorConstFieldAdvectionRTheta advection_field_rtheta,
            DTensor<CartesianBasis> const& advection_field_xy_centre,
            double dt) const
    {
        Kokkos::Profiling::pushRegion("PolarAdvection");
        IdxRangeRTheta grid(get_idx_range(allfdistribu));

        const int npoints_theta = IdxRangeTheta(grid).size();
        IdxRangeRTheta const grid_without_Opoint(grid.remove_first(IdxStepRTheta(1, 0)));
        IdxRangeRTheta const Opoint_grid(grid.take_first(IdxStepRTheta(1, npoints_theta)));


        // Convert advection field on RTheta to advection field on XY
        DVectorFieldMemAdvectionXY advection_field_xy_alloc(grid);
        DVectorFieldAdvectionXY advection_field_xy = get_field(advection_field_xy_alloc);

        LogicalToPhysicalMapping const& logical_to_physical_mapping_proxy
                = m_logical_to_physical_mapping;

        // (Ax, Ay) = J (Ar, Atheta)
        ddc::parallel_for_each(
                Kokkos::DefaultExecutionSpace(),
                grid_without_Opoint,
                KOKKOS_LAMBDA(IdxRTheta const irtheta) {
                    CoordRTheta const coord_rtheta(ddc::coordinate(irtheta));

                    Tensor J = logical_to_physical_mapping_proxy.jacobian_matrix(coord_rtheta);

                    ddcHelper::assign_vector_field_element(
                            advection_field_xy,
                            irtheta,
                            tensor_mul(
                                    index<'i', 'j'>(J),
                                    index<'j'>(advection_field_rtheta(irtheta))));
                });

        ddc::parallel_for_each(
                Kokkos::DefaultExecutionSpace(),
                Opoint_grid,
                KOKKOS_LAMBDA(IdxRTheta const irtheta) {
                    ddcHelper::assign_vector_field_element(
                            advection_field_xy,
                            irtheta,
                            advection_field_xy_centre);
                });

        (*this)(get_field(allfdistribu), get_const_field(advection_field_xy), dt);

        Kokkos::Profiling::popRegion();

        return allfdistribu;
    }
};
```


