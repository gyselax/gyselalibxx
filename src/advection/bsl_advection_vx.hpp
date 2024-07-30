// SPDX-License-Identifier: MIT

#pragma once
#include <ddc/ddc.hpp>

#include <ddc_helper.hpp>
#include <iinterpolator.hpp>
#include <species_info.hpp>

#include "ddc_aliases.hpp"
#include "iadvectionvx.hpp"

/**
 * @brief A class which computes the velocity advection along the dimension of interest GridV. Working for every cartesian geometry.
 */
template <class Geometry, class GridV>
class BslAdvectionVelocity : public IAdvectionVelocity<Geometry, GridV>
{
    using FdistribuIdxRange = typename Geometry::FdistribuDDom;
    using SpatialIdxRange = typename Geometry::SpatialDDom;
    using IdxV = Idx<GridV>;
    using DimV = typename GridV::continuous_dimension_type;

private:
    using PreallocatableInterpolatorType = interpolator_on_idx_range_t<
            IPreallocatableInterpolator,
            GridV,
            ddc::cartesian_prod_t<typename Geometry::SpatialDDom, typename Geometry::VelocityDDom>>;
    using InterpolatorType = interpolator_on_idx_range_t<
            IInterpolator,
            GridV,
            ddc::cartesian_prod_t<typename Geometry::SpatialDDom, typename Geometry::VelocityDDom>>;
    PreallocatableInterpolatorType const& m_interpolator_v;

public:
    /**
     * @brief Constructor 
     * @param[in] interpolator_v interpolator along the GridV direction which refers to the velocity space.  
     */
    explicit BslAdvectionVelocity(PreallocatableInterpolatorType const& interpolator_v)
        : m_interpolator_v(interpolator_v)
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
    Field<double, FdistribuIdxRange> operator()(
            Field<double, FdistribuIdxRange> const allfdistribu,
            Field<const double, SpatialIdxRange> const electric_field,
            double const dt) const override
    {
        Kokkos::Profiling::pushRegion("BslAdvectionVelocity");
        FdistribuIdxRange const dom = get_idx_range(allfdistribu);
        IdxRange<GridV> const v_dom = ddc::select<GridV>(dom);
        IdxRange<Species> const sp_dom = ddc::select<Species>(dom);

        FieldMem<double, typename InterpolatorType::batched_derivs_idx_range_type> derivs_min(
                m_interpolator_v.batched_derivs_idx_range_xmin(ddc::remove_dims_of<Species>(dom)));
        FieldMem<double, typename InterpolatorType::batched_derivs_idx_range_type> derivs_max(
                m_interpolator_v.batched_derivs_idx_range_xmax(ddc::remove_dims_of<Species>(dom)));
        ddc::parallel_fill(derivs_min, 0.);
        ddc::parallel_fill(derivs_max, 0.);

        // pre-allocate some memory to prevent allocation later in loop
        ddc::Chunk feet_coords_alloc(
                ddc::remove_dims_of(dom, sp_dom),
                ddc::DeviceAllocator<Coord<DimV>>());
        ddc::ChunkSpan feet_coords = get_field(feet_coords_alloc);
        std::unique_ptr<InterpolatorType> const interpolator_v_ptr = m_interpolator_v.preallocate();
        InterpolatorType const& interpolator_v = *interpolator_v_ptr;

        SpatialIdxRange const spatial_dom(get_idx_range(allfdistribu));

        auto batch_idx_range = ddc::remove_dims_of(ddc::remove_dims_of(dom, sp_dom), v_dom);
        using IdxB = typename decltype(batch_idx_range)::discrete_element_type;
        using IdxSpatial = typename SpatialIdxRange::discrete_element_type;

        ddc::for_each(sp_dom, [&](IdxSp const isp) {
            double const charge_proxy
                    = charge(isp); // TODO: consider proper way to access charge from device
            double const sqrt_me_on_mspecies = std::sqrt(mass(ielec()) / mass(isp));
            ddc::parallel_for_each(
                    Kokkos::DefaultExecutionSpace(),
                    batch_idx_range,
                    KOKKOS_LAMBDA(IdxB const ib) {
                        IdxSpatial const ix(ib);
                        // compute the displacement
                        double const dvx
                                = charge_proxy * sqrt_me_on_mspecies * dt * electric_field(ix);

                        // compute the coordinates of the feet
                        for (IdxV const iv : v_dom) {
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
