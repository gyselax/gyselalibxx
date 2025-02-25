// SPDX-License-Identifier: MIT
#pragma once
#include <ddc/ddc.hpp>

#include "ddc_alias_inline_functions.hpp"
#include "ddc_aliases.hpp"
#include "iadvectionx.hpp"
#include "iinterpolator.hpp"
#include "species_info.hpp"

/**
 * @brief A class which computes the spatial advection along the dimension of interest GridX. Working for every Cartesian geometry. 
 */
template <class Geometry, class GridX>
class BslAdvectionSpatial : public IAdvectionSpatial<Geometry, GridX>
{
    using GridV = typename Geometry::template velocity_dim_for<GridX>;
    using IdxRangeFdistrib = typename Geometry::IdxRangeFdistribu;
    using IdxX = Idx<GridX>;
    using IdxV = Idx<GridV>;
    using DimX = typename GridX::continuous_dimension_type;
    using DimV = typename GridV::continuous_dimension_type;
    using IdxRangeSpaceVelocity = ddc::cartesian_prod_t<
            typename Geometry::IdxRangeSpatial,
            typename Geometry::IdxRangeVelocity>;

private:
    using PreallocatableInterpolatorType = interpolator_on_idx_range_t<
            IPreallocatableInterpolator,
            GridX,
            IdxRangeSpaceVelocity>;
    using InterpolatorType
            = interpolator_on_idx_range_t<IInterpolator, GridX, IdxRangeSpaceVelocity>;
    PreallocatableInterpolatorType const& m_interpolator_x;

public:
    /**
     * @brief Constructor  
     * @param[in] interpolator_x interpolator along the GridX direction which refers to the spatial space.  
     */
    explicit BslAdvectionSpatial(PreallocatableInterpolatorType const& interpolator_x)
        : m_interpolator_x(interpolator_x)
    {
    }

    ~BslAdvectionSpatial() override = default;

    /**
     * @brief Advects fdistribu along GridX for a duration dt.
     * @param[in, out] allfdistribu Reference to the whole distribution function for one species, allocated on the device (ie it lets the choice of the location depend on the build configuration).
     * @param[in] dt Time step
     * @return A reference to the allfdistribu array containing the value of the function at the coordinates.
     */
    Field<double, IdxRangeFdistrib> operator()(
            Field<double, IdxRangeFdistrib> const allfdistribu,
            double const dt) const override
    {
        using IdxRangeBatch = ddc::remove_dims_of_t<IdxRangeFdistrib, Species, GridX>;
        using IdxBatch = typename IdxRangeBatch::discrete_element_type;

        Kokkos::Profiling::pushRegion("BslAdvectionSpatial");
        IdxRangeFdistrib const idx_range = get_idx_range(allfdistribu);
        IdxRange<GridX> const x_idx_range = ddc::select<GridX>(idx_range);
        IdxRange<Species> const sp_idx_range = ddc::select<Species>(idx_range);

        // pre-allocate some memory to prevent allocation later in loop
        IdxRangeSpaceVelocity batched_feet_idx_range(idx_range);
        FieldMem<Coord<DimX>, IdxRangeSpaceVelocity> feet_coords_alloc(batched_feet_idx_range);
        Field<Coord<DimX>, IdxRangeSpaceVelocity> feet_coords(get_field(feet_coords_alloc));
        std::unique_ptr<InterpolatorType> const interpolator_x_ptr = m_interpolator_x.preallocate();
        InterpolatorType const& interpolator_x = *interpolator_x_ptr;

        IdxRangeBatch batch_idx_range(idx_range);

        for (IdxSp const isp : sp_idx_range) {
            double const sqrt_me_on_mspecies = std::sqrt(mass(ielec()) / mass(isp));
            ddc::parallel_for_each(
                    Kokkos::DefaultExecutionSpace(),
                    batch_idx_range,
                    KOKKOS_LAMBDA(IdxBatch const ib) {
                        // compute the displacement
                        IdxV const iv(ib);
                        Coord<DimV> const coord_iv = ddc::coordinate(iv);
                        double const dx = sqrt_me_on_mspecies * dt * coord_iv;

                        // compute the coordinates of the feet
                        for (IdxX const ix : x_idx_range) {
                            feet_coords(ix, ib) = Coord<DimX>(ddc::coordinate(ix) - dx);
                        }
                    });
            interpolator_x(allfdistribu[isp], get_const_field(feet_coords));
        }

        Kokkos::Profiling::popRegion();
        return allfdistribu;
    }
};
