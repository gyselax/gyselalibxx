// SPDX-License-Identifier: MIT

#pragma once

#include <ddc/ddc.hpp>

#include <i_interpolator_batched.hpp>
#include <species_info.hpp>

#include "iadvectionx.hpp"

/**
 * @brief A class which computes the spatial advection along the dimension of interest DDimX. Working for every cartesian geometry. 
 */
template <class Geometry, class DDimX>
class BslAdvectionSpatialBatched : public IAdvectionSpatial<Geometry, DDimX>
{
    using DDimSp = typename Geometry::DDimSp;
    using DDimV = typename Geometry::template velocity_dim_for<DDimX>;
    using DDom = typename Geometry::FdistribuDDom;
    using DElemX = ddc::DiscreteElement<DDimX>;
    using DElemV = ddc::DiscreteElement<DDimV>;
    using DElemSp = ddc::DiscreteElement<DDimSp>;
    using DElemSpV = ddc::DiscreteElement<DDimSp, DDimV>;
    using CDimX = typename DDimX::continuous_dimension_type;

private:
    using PreallocatableInterpolatorType = interpolator_on_domain_t<
            IPreallocatableInterpolatorBatched,
            DDimX,
            ddc::cartesian_prod_t<typename Geometry::SpatialDDom, typename Geometry::VelocityDDom>>;
    using InterpolatorType = interpolator_on_domain_t<
            IInterpolatorBatched,
            DDimX,
            ddc::cartesian_prod_t<typename Geometry::SpatialDDom, typename Geometry::VelocityDDom>>;
    PreallocatableInterpolatorType const& m_interpolator_x;

public:
    /**
     * @brief Constructor  
     * @param[in] interpolator_x interpolator along the DDimX direction which refers to the spatial space.  
     */
    explicit BslAdvectionSpatialBatched(PreallocatableInterpolatorType const& interpolator_x)
        : m_interpolator_x(interpolator_x)
    {
    }

    ~BslAdvectionSpatialBatched() override = default;

    /**
     * @brief Advects fdistribu along DDimX for a duration dt.
     * @param[in, out] allfdistribu Reference to the whole distribution function for one species, allocated on the device (ie it lets the choice of the location depend on the build configuration).
     * @param[in] dt Time step
     * @return A reference to the allfdistribu array containing the value of the function at the coordinates.
     */
    device_t<ddc::ChunkSpan<double, DDom>> operator()(
            device_t<ddc::ChunkSpan<double, DDom>> const allfdistribu,
            double const dt) const override
    {
        Kokkos::Profiling::pushRegion("BslAdvectionSpatial");
        DDom const dom = allfdistribu.domain();
        ddc::DiscreteDomain<DDimX> const x_dom = ddc::select<DDimX>(dom);
        ddc::DiscreteDomain<DDimV> const v_dom = ddc::select<DDimV>(dom);
        ddc::DiscreteDomain<DDimSp> const sp_dom = ddc::select<DDimSp>(dom);

        // pre-allocate some memory to prevent allocation later in loop
        ddc::Chunk feet_coords_alloc(
                ddc::remove_dims_of(dom, sp_dom),
                ddc::DeviceAllocator<ddc::Coordinate<CDimX>>());
        ddc::ChunkSpan feet_coords = feet_coords_alloc.span_view();
        std::unique_ptr<InterpolatorType> const interpolator_x_ptr = m_interpolator_x.preallocate();
        InterpolatorType const& interpolator_x = *interpolator_x_ptr;

        auto c_dom = ddc::remove_dims_of(
                dom,
                ddc::DiscreteDomain<DDimSp, DDimX, DDimV>(sp_dom, x_dom, v_dom));
        using DElemC = typename decltype(c_dom)::discrete_element_type;

        for (DElemSp const isp : sp_dom) {
            double const sqrt_me_on_mspecies = std::sqrt(mass(ielec()) / mass(isp));
            ddc::for_each(
                    ddc::policies::parallel_device,
                    c_dom,
                    KOKKOS_LAMBDA(DElemC const ic) {
                        for (DElemV const iv : v_dom) {
                            // compute the displacement
                            double const dx = sqrt_me_on_mspecies * dt * ddc::coordinate(iv);

                            // compute the coordinates of the feet
                            for (DElemX const ix : x_dom) {
                                feet_coords(ix, iv, ic)
                                        = ddc::Coordinate<CDimX>(ddc::coordinate(ix) - dx);
                            }
                        }
                    });
            interpolator_x(allfdistribu[isp], feet_coords.span_cview());
        }

        Kokkos::Profiling::popRegion();
        return allfdistribu;
    }
};
