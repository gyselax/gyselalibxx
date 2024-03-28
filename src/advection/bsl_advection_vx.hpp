// SPDX-License-Identifier: MIT

#pragma once

#include <ddc/ddc.hpp>

#include <ddc_helper.hpp>
#include <iinterpolator.hpp>
#include <species_info.hpp>

#include "iadvectionvx.hpp"

/**
 * @brief A class which computes the velocity advection along the dimension of interest DDimV.
 */
template <class Geometry, class DDimV>
class BslAdvectionVelocity : public IAdvectionVelocity<Geometry, DDimV>
{
    using DDimSp = typename Geometry::DDimSp;
    using FdistribuDDom = typename Geometry::FdistribuDDom;
    using SpatialDDom = typename Geometry::SpatialDDom;
    using DElemV = ddc::DiscreteElement<DDimV>;
    using DElemSp = ddc::DiscreteElement<DDimSp>;
    using CDimV = typename DDimV::continuous_dimension_type;

private:
    IPreallocatableInterpolator<DDimV> const& m_interpolator_v;

public:
    /**
     * @brief Constructor 
     * @param[in] interpolator_v interpolator along the DDimV direction which refers to the velocity space.  
     */
    explicit BslAdvectionVelocity(IPreallocatableInterpolator<DDimV> const& interpolator_v)
        : m_interpolator_v(interpolator_v)
    {
    }

    ~BslAdvectionVelocity() override = default;

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

        auto allfdistribu_host_alloc = ddc::create_mirror_view_and_copy(allfdistribu);
        ddc::ChunkSpan allfdistribu_host = allfdistribu_host_alloc.span_view();
        auto electric_field_host_alloc = ddc::create_mirror_view_and_copy(electric_field);
        ddc::ChunkSpan electric_field_host = electric_field_host_alloc.span_view();

        FdistribuDDom const dom = allfdistribu.domain();
        ddc::DiscreteDomain<DDimV> const v_dom = ddc::select<DDimV>(dom);
        ddc::DiscreteDomain<DDimSp> const sp_dom = ddc::select<DDimSp>(dom);

        // pre-allocate some memory to prevent allocation later in loop
        ddc::Chunk<ddc::Coordinate<CDimV>, ddc::DiscreteDomain<DDimV>> feet_coords(v_dom);
        ddc::Chunk<double, ddc::DiscreteDomain<DDimV>> contiguous_slice(v_dom);
        std::unique_ptr<IInterpolator<DDimV>> const interpolator_v_ptr
                = m_interpolator_v.preallocate();
        IInterpolator<DDimV> const& interpolator_v = *interpolator_v_ptr;

        SpatialDDom const spatial_dom(allfdistribu.domain());

        auto c_dom = ddc::remove_dims_of(
                ddc::remove_dims_of(ddc::remove_dims_of(dom, sp_dom), spatial_dom),
                v_dom);
        using DElemC = typename decltype(c_dom)::discrete_element_type;
        using DElemSpatial = typename SpatialDDom::discrete_element_type;

        ddc::for_each(c_dom, [&](DElemC const ic) {
            ddc::for_each(sp_dom, [&](DElemSp const isp) {
                double const sqrt_me_on_mspecies = std::sqrt(mass(ielec()) / mass(isp));

                ddc::for_each(spatial_dom, [&](DElemSpatial const ix) {
                    // compute the displacement
                    double const dvx
                            = charge(isp) * sqrt_me_on_mspecies * dt * electric_field_host(ix);

                    // compute the coordinates of the feet
                    ddc::for_each(v_dom, [&](DElemV const iv) {
                        feet_coords(iv) = ddc::Coordinate<CDimV>(ddc::coordinate(iv) - dvx);
                    });

                    // copy the slice in contiguous memory
                    ddc::parallel_deepcopy(contiguous_slice, allfdistribu_host[ic][isp][ix]);

                    // build a spline representation of the data
                    interpolator_v(contiguous_slice, feet_coords.span_cview());

                    // copy back
                    ddc::parallel_deepcopy(allfdistribu_host[ic][isp][ix], contiguous_slice);
                });
            });
        });

        ddc::parallel_deepcopy(allfdistribu, allfdistribu_host);
        Kokkos::Profiling::popRegion();
        return allfdistribu;
    }
};
