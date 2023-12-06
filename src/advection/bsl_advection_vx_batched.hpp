// SPDX-License-Identifier: MIT

#pragma once

#include <ddc/ddc.hpp>

#include <ddc_helper.hpp>
#include <i_interpolator_batched.hpp>
#include <species_info.hpp>

#include "iadvectionvx.hpp"

/**
 * @brief A class which computes the velocity advection along the dimension of interest DDimV. Working for every cartesian geometry.
 */
template <class Geometry, class DDimV>
class BslAdvectionVelocityBatched : public IAdvectionVelocity<Geometry, DDimV>
{
    using DDimSp = typename Geometry::DDimSp;
    using FdistribuDDom = typename Geometry::FdistribuDDom;
    using SpatialDDom = typename Geometry::SpatialDDom;
    using DElemV = ddc::DiscreteElement<DDimV>;
    using DElemSp = ddc::DiscreteElement<DDimSp>;
    using CDimV = typename DDimV::continuous_dimension_type;

private:
    using PreallocatableInterpolatorType = interpolator_on_domain_t<
            IPreallocatableInterpolatorBatched,
            DDimV,
            ddc::cartesian_prod_t<typename Geometry::SpatialDDom, typename Geometry::VelocityDDom>>;
    using InterpolatorType = interpolator_on_domain_t<
            IInterpolatorBatched,
            DDimV,
            ddc::cartesian_prod_t<typename Geometry::SpatialDDom, typename Geometry::VelocityDDom>>;
    PreallocatableInterpolatorType const& m_interpolator_v;

public:
    /**
     * @brief Constructor 
     * @param[in] interpolator_v interpolator along the DDimV direction which refers to the velocity space.  
     */
    explicit BslAdvectionVelocityBatched(PreallocatableInterpolatorType const& interpolator_v)
        : m_interpolator_v(interpolator_v)
    {
    }

    ~BslAdvectionVelocityBatched() override = default;

    /**
     * @brief Advects fdistribu along DDimV for a duration dt.
     * @param[in, out] allfdistribu Reference to the whole distribution function for one species, allocated on the device (ie it lets the choice of the location depend on the build configuration).
     * @param[in] electric_field Reference to the electric field which derives from electrostatic potential, allocated on the device.
     * @param[in] dt Time step
     * @return A reference to the allfdistribu array containing the value of the function at the coordinates.
     */
    device_t<ddc::ChunkSpan<double, FdistribuDDom>> operator()(
            device_t<ddc::ChunkSpan<double, FdistribuDDom>> const allfdistribu,
            device_t<ddc::ChunkSpan<const double, SpatialDDom>> const electric_field,
            double const dt) const override
    {
        Kokkos::Profiling::pushRegion("BslAdvectionVelocity");
        FdistribuDDom const dom = allfdistribu.domain();
        ddc::DiscreteDomain<DDimV> const v_dom = ddc::select<DDimV>(dom);
        ddc::DiscreteDomain<DDimSp> const sp_dom = ddc::select<DDimSp>(dom);

        // pre-allocate some memory to prevent allocation later in loop
        ddc::Chunk feet_coords_alloc(
                ddc::remove_dims_of(dom, sp_dom),
                ddc::DeviceAllocator<ddc::Coordinate<CDimV>>());
        ddc::ChunkSpan feet_coords = feet_coords_alloc.span_view();
        std::unique_ptr<InterpolatorType> const interpolator_v_ptr = m_interpolator_v.preallocate();
        InterpolatorType const& interpolator_v = *interpolator_v_ptr;

        SpatialDDom const spatial_dom(allfdistribu.domain());

        auto c_dom = ddc::remove_dims_of(ddc::remove_dims_of(dom, sp_dom), v_dom);
        using DElemC = typename decltype(c_dom)::discrete_element_type;
        using DElemSpatial = typename SpatialDDom::discrete_element_type;

        ddc::for_each(sp_dom, [&](DElemSp const isp) {
            double const charge_proxy
                    = charge(isp); // TODO: consider proper way to access charge from device
            double const sqrt_me_on_mspecies = std::sqrt(mass(ielec()) / mass(isp));
            ddc::for_each(
                    ddc::policies::parallel_device,
                    c_dom,
                    DDC_LAMBDA(DElemC const ic) {
                        DElemSpatial const ix(ic);
                        // compute the displacement
                        double const dvx
                                = charge_proxy * sqrt_me_on_mspecies * dt * electric_field(ix);

                        // compute the coordinates of the feet
                        for (DElemV const iv : v_dom) {
                            feet_coords(iv, ic) = ddc::Coordinate<CDimV>(ddc::coordinate(iv) - dvx);
                        }
                    });
            interpolator_v(allfdistribu[isp], feet_coords.span_cview());
        });

        Kokkos::Profiling::popRegion();
        return allfdistribu;
    }
};
