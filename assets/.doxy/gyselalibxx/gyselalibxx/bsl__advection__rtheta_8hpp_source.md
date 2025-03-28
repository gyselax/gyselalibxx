

# File bsl\_advection\_rtheta.hpp

[**File List**](files.md) **>** [**advection**](dir_18bb63d4be19d3ea733e61d8625caf4d.md) **>** [**bsl\_advection\_rtheta.hpp**](bsl__advection__rtheta_8hpp.md)

[Go to the documentation of this file](bsl__advection__rtheta_8hpp.md)


```C++
// SPDX-License-Identifier: MIT
#pragma once

#include "ddc_alias_inline_functions.hpp"
#include "ddc_aliases.hpp"
#include "geometry.hpp"
#include "i_interpolator_2d.hpp"
#include "iadvection_rtheta.hpp"
#include "indexed_tensor.hpp"
#include "metric_tensor_evaluator.hpp"
#include "spline_interpolator_2d.hpp"
#include "spline_polar_foot_finder.hpp"
#include "vector_field.hpp"
#include "vector_field_mem.hpp"



template <class FootFinder, class Mapping>
class BslAdvectionRTheta : public IAdvectionRTheta
{
private:
    using evaluator_type = ddc::SplineEvaluator2D<
            Kokkos::DefaultExecutionSpace,
            Kokkos::DefaultExecutionSpace::memory_space,
            BSplinesR,
            BSplinesTheta,
            GridR,
            GridTheta,
            ddc::NullExtrapolationRule,
            ddc::NullExtrapolationRule,
            ddc::PeriodicExtrapolationRule<Theta>,
            ddc::PeriodicExtrapolationRule<Theta>,
            GridR,
            GridTheta>;
    using PreallocatableSplineInterpolatorType
            = PreallocatableSplineInterpolator2D<SplineRThetaBuilder, evaluator_type>;
    PreallocatableSplineInterpolatorType const& m_interpolator;

    FootFinder const& m_find_feet;

    Mapping const& m_mapping;


public:
    BslAdvectionRTheta(
            PreallocatableSplineInterpolatorType const& function_interpolator,
            FootFinder const& foot_finder,
            Mapping const& mapping)
        : m_interpolator(function_interpolator)
        , m_find_feet(foot_finder)
        , m_mapping(mapping)
    {
    }

    ~BslAdvectionRTheta() override = default;


    DFieldRTheta operator()(
            DFieldRTheta allfdistribu,
            DConstVectorFieldRTheta<X, Y> advection_field_xy,
            double dt) const override
    {
        // Pre-allocate some memory to prevent allocation later in loop
        std::unique_ptr<IInterpolator2D<IdxRangeRTheta, IdxRangeRTheta>> const interpolator_ptr
                = m_interpolator.preallocate();

        // Initialise the feet
        FieldMemRTheta<CoordRTheta> feet_rtheta_alloc(get_idx_range(advection_field_xy));
        FieldRTheta<CoordRTheta> feet_rtheta = get_field(feet_rtheta_alloc);
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


    DFieldRTheta operator()(
            DFieldRTheta allfdistribu,
            DConstVectorFieldRTheta<R, Theta> advection_field_rtheta,
            CoordXY const& advection_field_xy_centre,
            double dt) const override
    {
        Kokkos::Profiling::pushRegion("PolarAdvection");
        IdxRangeRTheta grid(get_idx_range<GridR, GridTheta>(allfdistribu));

        const int npoints_theta = IdxRangeTheta(grid).size();
        IdxRangeRTheta const grid_without_Opoint(grid.remove_first(IdxStepRTheta(1, 0)));
        IdxRangeRTheta const Opoint_grid(grid.take_first(IdxStepRTheta(1, npoints_theta)));


        // Convert advection field on RTheta to advection field on XY
        DVectorFieldMemRTheta<X, Y> advection_field_xy_alloc(grid);
        DVectorFieldRTheta<X, Y> advection_field_xy = get_field(advection_field_xy_alloc);

        Mapping const& mapping_proxy = m_mapping;

        // (Ax, Ay) = J (Ar, Atheta)
        ddc::parallel_for_each(
                Kokkos::DefaultExecutionSpace(),
                grid_without_Opoint,
                KOKKOS_LAMBDA(IdxRTheta const irtheta) {
                    CoordRTheta const coord_rtheta(ddc::coordinate(irtheta));

                    Tensor J = mapping_proxy.jacobian_matrix(coord_rtheta);

                    DVector<X, Y> advec_field_xy = tensor_mul(
                            index<'i', 'j'>(J),
                            index<'j'>(advection_field_rtheta(irtheta)));
                    ddcHelper::get<X>(advection_field_xy)(irtheta)
                            = ddcHelper::get<X>(advec_field_xy);
                    ddcHelper::get<Y>(advection_field_xy)(irtheta)
                            = ddcHelper::get<Y>(advec_field_xy);
                });

        ddc::parallel_for_each(
                Kokkos::DefaultExecutionSpace(),
                Opoint_grid,
                KOKKOS_LAMBDA(IdxRTheta const irtheta) {
                    ddcHelper::get<X>(advection_field_xy)(irtheta)
                            = CoordX(advection_field_xy_centre);
                    ddcHelper::get<Y>(advection_field_xy)(irtheta)
                            = CoordY(advection_field_xy_centre);
                });

        (*this)(get_field(allfdistribu), get_const_field(advection_field_xy), dt);

        Kokkos::Profiling::popRegion();

        return allfdistribu;
    }
};
```


