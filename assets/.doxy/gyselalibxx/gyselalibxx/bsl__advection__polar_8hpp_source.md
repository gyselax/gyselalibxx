

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
#include "l_norm_tools.hpp"
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

    using DimX = typename LogicalToPhysicalMapping::cartesian_tag_x;
    using DimY = typename LogicalToPhysicalMapping::cartesian_tag_y;

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

    using DFieldFDistribu = DField<IdxRangeBatched, MemorySpace>;

    using CFieldMemFeetRTheta = FieldMem<CoordRTheta, IdxRangeBatched, MemorySpace>;

    using CFieldFeetRTheta = Field<CoordRTheta, IdxRangeBatched, MemorySpace>;

    using DVectorFieldMemAdvectionXY
            = DVectorFieldMem<IdxRangeBatched, CartesianBasis, MemorySpace>;
    using DVectorFieldAdvectionXY = DVectorField<IdxRangeBatched, CartesianBasis, MemorySpace>;
    using DVectorConstFieldAdvectionXY
            = DVectorConstField<IdxRangeBatched, CartesianBasis, MemorySpace>;

    using DVectorFieldMemAdvectionXYOnBatch
            = DVectorFieldMem<IdxRangeBatch, CartesianBasis, MemorySpace>;
    using DVectorFieldAdvectionXYOnBatch = DVectorField<IdxRangeBatch, CartesianBasis, MemorySpace>;

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
        std::unique_ptr<IInterpolator2D<IdxRangeRTheta, IdxRangeBatched>> const interpolator_ptr
                = m_interpolator.preallocate();

        // Initialise the feet
        CFieldMemFeetRTheta feet_rtheta_alloc(get_idx_range(advection_field_xy));
        CFieldFeetRTheta feet_rtheta = get_field(feet_rtheta_alloc);
        ddc::parallel_for_each(
                ExecSpace(),
                get_idx_range(advection_field_xy),
                KOKKOS_LAMBDA(IdxBatched const idx) {
                    IdxRTheta const irtheta(idx);
                    feet_rtheta(idx) = ddc::coordinate(irtheta);
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
        using IdxRangeBatchedWithoutR = ddc::remove_dims_of_t<IdxRangeBatched, GridR>;
        Kokkos::Profiling::pushRegion("PolarAdvection");
        IdxRangeBatched grid(get_idx_range(allfdistribu));
        IdxRangeR radial_grid(grid);
        IdxRangeBatchedWithoutR no_r_grid(grid);

        // Check the first points on R correspond to the O-point.
        assert(ddc::coordinate(radial_grid.front()) < 1e-13);

        IdxRangeBatched const
                grid_without_Opoint(radial_grid.remove_first(IdxStep<GridR>(1)), no_r_grid);
        IdxRangeBatched const Opoint_grid(radial_grid.take_first(IdxStepR(1)), no_r_grid);


        // Convert advection field on RTheta to advection field on XY
        DVectorFieldMemAdvectionXY advection_field_xy_alloc(grid);
        DVectorFieldAdvectionXY advection_field_xy = get_field(advection_field_xy_alloc);

        LogicalToPhysicalMapping const& logical_to_physical_mapping_proxy
                = m_logical_to_physical_mapping;

        // (Ax, Ay) = J (Ar, Atheta)
        copy_to_vector_space<CartesianBasis>(
                ExecSpace(),
                advection_field_xy[grid_without_Opoint],
                logical_to_physical_mapping_proxy,
                advection_field_rtheta[grid_without_Opoint]);

        ddc::parallel_for_each(
                ExecSpace(),
                Opoint_grid,
                KOKKOS_LAMBDA(IdxBatched const idx) {
                    ddcHelper::assign_vector_field_element(
                            advection_field_xy,
                            idx,
                            advection_field_xy_centre);
                });

        (*this)(get_field(allfdistribu), get_const_field(advection_field_xy), dt);

        Kokkos::Profiling::popRegion();

        return allfdistribu;
    }

    DFieldFDistribu operator()(
            DFieldFDistribu allfdistribu,
            DVectorConstFieldAdvectionRTheta advection_field_rtheta,
            double dt) const
    {
        using IdxRangeBatchedWithoutR = ddc::remove_dims_of_t<IdxRangeBatched, GridR>;
        Kokkos::Profiling::pushRegion("PolarAdvection");
        IdxRangeBatched grid(get_idx_range(allfdistribu));
        IdxRangeR radial_grid(grid);
        IdxRangeTheta theta_grid(grid);
        IdxRangeBatchedWithoutR no_r_grid(grid);
        IdxRangeBatch no_rtheta_grid(grid);

        LogicalToPhysicalMapping const& logical_to_physical_mapping_proxy
                = m_logical_to_physical_mapping;

        // Test if the first points on R correspond to the O-point.
        bool const first_row_is_o_point = (ddc::coordinate(radial_grid.front()) < 1e-13);

        /*
            If the O-point is not considered as on the grid, we have to be sure that the grid of 
            the advection field on (R, Theta) matches with the grid of the advected function.  
        */
        assert(first_row_is_o_point
               || (IdxRangeR(get_idx_range(advection_field_rtheta)) == radial_grid));

        IdxRangeBatched const grid_without_Opoint(
                radial_grid.remove_first(IdxStep<GridR>(first_row_is_o_point)),
                no_r_grid);

        // Convert advection field on RTheta to advection field on XY
        DVectorFieldMemAdvectionXY advection_field_xy_alloc(grid);
        DVectorFieldAdvectionXY advection_field_xy = get_field(advection_field_xy_alloc);

        // (Ax, Ay) = J (Ar, Atheta)
        copy_to_vector_space<CartesianBasis>(
                ExecSpace(),
                advection_field_xy[grid_without_Opoint],
                logical_to_physical_mapping_proxy,
                advection_field_rtheta[grid_without_Opoint]);

        // Treatment for the O-point.
        if (first_row_is_o_point) {
            IdxRangeBatched const Opoint_grid(radial_grid.take_first(IdxStepR(1)), no_r_grid);
            IdxRangeRTheta const grid_first_ring(
                    radial_grid.take_first(IdxStepR(2)).remove_first(IdxStep<GridR>(1)),
                    theta_grid);


            // Jacobian ill-defined at the O-point, we average the values around the O-point,
            std::size_t ntheta_points = theta_grid.size();
            ddc::parallel_for_each(
                    ExecSpace(),
                    no_rtheta_grid,
                    KOKKOS_LAMBDA(IdxBatch const idx_batch) {
                        DTensor<CartesianBasis> advection_field_xy_average_centre(0, 0);

                        IdxR const idx_r(grid_first_ring.front()); // one ring => one r index.
                        for (IdxTheta const idx_theta : IdxRangeTheta(grid_first_ring)) {
                            ddcHelper::get<DimX>(advection_field_xy_average_centre)
                                    += ddcHelper::get<DimX>(
                                            advection_field_xy)(idx_batch, idx_r, idx_theta);
                            ddcHelper::get<DimY>(advection_field_xy_average_centre)
                                    += ddcHelper::get<DimY>(
                                            advection_field_xy)(idx_batch, idx_r, idx_theta);
                        }

                        ddcHelper::get<DimX>(advection_field_xy_average_centre) /= ntheta_points;
                        ddcHelper::get<DimY>(advection_field_xy_average_centre) /= ntheta_points;

                        // and assign the averaged value to all the points at the O-point.
                        IdxR const idx_r_Opt(radial_grid.front());
                        for (IdxTheta const idx_theta : IdxRangeTheta(grid_first_ring)) {
                            ddcHelper::get<DimX>(
                                    advection_field_xy)(idx_batch, idx_r_Opt, idx_theta)
                                    = ddcHelper::get<DimX>(advection_field_xy_average_centre);
                            ddcHelper::get<DimY>(
                                    advection_field_xy)(idx_batch, idx_r_Opt, idx_theta)
                                    = ddcHelper::get<DimY>(advection_field_xy_average_centre);
                        }
                    });
        }

        (*this)(get_field(allfdistribu), get_const_field(advection_field_xy), dt);

        Kokkos::Profiling::popRegion();

        return allfdistribu;
    }
};
```


