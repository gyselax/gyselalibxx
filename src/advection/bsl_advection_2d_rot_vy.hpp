// SPDX-License-Identifier: MIT
#pragma once
#include <ddc/ddc.hpp>

#include "ddc_alias_inline_functions.hpp"
#include "ddc_aliases.hpp"
#include "ddc_helper.hpp"
#include "iadvectionrot2d.hpp"
#include "iinterpolator.hpp"
#include "species_info.hpp"
#include <cmath>


/**
 * @brief A class which computes the velocity advection in vy of a 2D shifted rotation. Working for every cartesian geometry.
 */
template <class Geometry, class GridV>
class BslAdvectionVelocityRot2DVy : public IAdvectionVelocityRot2D<Geometry, GridV>
{
    using IdxRangeFdistribu = typename Geometry::IdxRangeFdistribu;
    using IdxRangeSpatial = typename Geometry::IdxRangeSpatial;
    using IdxSpatial = typename IdxRangeSpatial::discrete_element_type;
    using IdxV = Idx<GridV>;
    using DimV = typename GridV::continuous_dimension_type;
    using IdxRangeSpaceVelocity
            = ddc::remove_dims_of_t<typename Geometry::IdxRangeFdistribu, Species>;

    using GridVx = typename Geometry::template the_other_velocity_dim<GridV>;
    using IdxVx = Idx<GridVx>;

private:
    using PreallocatableInterpolatorType = interpolator_on_idx_range_t<
            IPreallocatableInterpolator,
            GridV,
            IdxRangeSpaceVelocity>;
    using InterpolatorType
            = interpolator_on_idx_range_t<IInterpolator, GridV, IdxRangeSpaceVelocity>;
    PreallocatableInterpolatorType const& m_interpolator_v;

public:
    /**
     * @brief Constructor 
     * @param[in] interpolator_v interpolator along the GridV direction which refers to the velocity space.  
     */
    explicit BslAdvectionVelocityRot2DVy(PreallocatableInterpolatorType const& interpolator_v)
        : m_interpolator_v(interpolator_v)
    {
    }

    ~BslAdvectionVelocityRot2DVy() override = default;

    /**
     * @brief Advects fdistribu along GridV for a duration dt.
     * @param[in, out] allfdistribu Reference to the whole distribution function for one species, allocated on the device (ie it lets the choice of the location depend on the build configuration).
     * @param[in] magnetic_field_z Reference to the magnetic field.
     * @param[in] mean_velocity_x Reference to the velocity shift in vx.
     * @param[in] mean_velocity_y Reference to the velocity shift in vy.
     * @param[in] dt Time step
     * @return A reference to the allfdistribu array containing the value of the function at the coordinates.
     */
    Field<double, IdxRangeFdistribu> operator()(
            Field<double, IdxRangeFdistribu> const allfdistribu,
            Field<const double, IdxRangeSpatial> const magnetic_field_z,
            Field<const double, IdxRangeSpatial> const mean_velocity_x,
            Field<const double, IdxRangeSpatial> const mean_velocity_y,
            double const dt) const override
    {
        using IdxRangeBatch = ddc::remove_dims_of_t<IdxRangeFdistribu, Species, GridV>;
        using IdxBatch = typename IdxRangeBatch::discrete_element_type;

        Kokkos::Profiling::pushRegion("BslAdvectionVelocityRot2DVy");
        IdxRangeFdistribu const idx_range = get_idx_range(allfdistribu);
        IdxRange<GridV> const idx_range_v = ddc::select<GridV>(idx_range);
        IdxRange<Species> const idx_range_sp = ddc::select<Species>(idx_range);

        FieldMem<double, typename InterpolatorType::batched_derivs_idx_range_type> derivs_min(
                m_interpolator_v.batched_derivs_idx_range_xmin(
                        ddc::remove_dims_of<Species>(idx_range)));
        FieldMem<double, typename InterpolatorType::batched_derivs_idx_range_type> derivs_max(
                m_interpolator_v.batched_derivs_idx_range_xmax(
                        ddc::remove_dims_of<Species>(idx_range)));
        ddc::parallel_fill(derivs_min, 0.);
        ddc::parallel_fill(derivs_max, 0.);

        // pre-allocate some memory to prevent allocation later in loop
        IdxRangeSpaceVelocity batched_feet_idx_range(idx_range);
        FieldMem<Coord<DimV>, IdxRangeSpaceVelocity> feet_coords_alloc(batched_feet_idx_range);
        Field<Coord<DimV>, IdxRangeSpaceVelocity> feet_coords(get_field(feet_coords_alloc));
        std::unique_ptr<InterpolatorType> const interpolator_v_ptr = m_interpolator_v.preallocate();
        InterpolatorType const& interpolator_v = *interpolator_v_ptr;

        IdxRangeSpatial const idx_range_spatial(get_idx_range(allfdistribu));

        IdxRangeBatch batch_idx_range(idx_range);

        ddc::for_each(idx_range_sp, [&](IdxSp const isp) {
            double const charge_proxy
                    = charge(isp); // TODO: consider proper way to access charge from device
            //double const sqrt_me_on_mspecies = std::sqrt(mass(ielec()) / mass(isp));
            ddc::parallel_for_each(
                    Kokkos::DefaultExecutionSpace(),
                    batch_idx_range,
                    KOKKOS_LAMBDA(IdxBatch const ib) {
                        IdxSpatial const ix(ib);
                        IdxVx const ivx(ib);
                        // compute the displacement
                        //double const dvx
                        //        = -1.0 * charge_proxy * sqrt_me_on_mspecies * (ddc::coordinate(ivx) - mean_velocity_x(ix)) * std::sin(dt*magnetic_field_z(ix)) ;
                        //double const dvx
                               // = charge_proxy / mass(isp) * (ddc::coordinate(ivx) + mean_velocity_x(ix)) * std::sin(-dt*magnetic_field_z(ix)) ;
                        double const dvx
                                = (ddc::coordinate(ivx) + mean_velocity_x(ix)) * std::sin(-dt*charge_proxy/mass(isp)*magnetic_field_z(ix)) ;
                        //std::cout << "check rot 2d vy is: " << std::setprecision(16) <<  mean_velocity_x(ix) << std::endl;
                        // compute the coordinates of the feet
                        for (IdxV const iv : idx_range_v) {
                            feet_coords(iv, ib) = Coord<DimV>(ddc::coordinate(iv) - dvx);
                        }
                    });
            interpolator_v(
                    allfdistribu[isp],
                    get_const_field(feet_coords),
                    get_const_field(derivs_min),
                    get_const_field(derivs_max));
        });

        Kokkos::Profiling::popRegion();
        return allfdistribu;
    }
};

