#pragma once

#include <vector>

#include <ddc/ddc.hpp>

#include <sll/spline_boundary_conditions.hpp>

/**
 * @brief Helper class for the initialisation of the mesh of interpolation points.
 *
 * A helper class for the initialisation of the mesh of interpolation points. This
 * class should be used when the interpolation points should be located at the
 * knots of the spline. This is possible with any kind of boundary condition except
 * Greville boundary conditions (as there will not be enough interpolation points).
 * In the case of strongly non-uniform splines this choice may result in a less
 * well conditioned problem, however most mathematical stability results are proven
 * with this choice of interpolation points.
 */
template <class BSplines, BoundCond BcXmin, BoundCond BcXmax>
class KnotsAsInterpolationPoints
{
    static_assert(BcXmin != BoundCond::GREVILLE);
    static_assert(BcXmax != BoundCond::GREVILLE);

    using tag_type = typename BSplines::tag_type;

public:
    /**
     * Get the sampling of interpolation points.
     *
     * @return sampling The DDC point sampling of the interpolation points.
     */
    template <typename U = BSplines>
    static auto get_sampling()
    {
        if constexpr (U::is_uniform()) {
            using Sampling = ddc::UniformPointSampling<tag_type>;
            return std::get<0>(
                    Sampling::
                            init(ddc::discrete_space<BSplines>().rmin(),
                                 ddc::discrete_space<BSplines>().rmax(),
                                 ddc::DiscreteVector<ddc::UniformPointSampling<tag_type>>(
                                         ddc::discrete_space<BSplines>().ncells() + 1)));
        } else {
            using Sampling = ddc::NonUniformPointSampling<tag_type>;
            using SamplingImpl = typename Sampling::template Impl<Kokkos::HostSpace>;
            std::vector<double> knots(ddc::discrete_space<BSplines>().npoints());
            for (int i(0); i < ddc::discrete_space<BSplines>().npoints(); ++i) {
                knots[i] = ddc::discrete_space<BSplines>().get_knot(i);
            }
            return SamplingImpl(knots);
        }
    }

    /// The DDC type of the sampling for the interpolation points.
    using interpolation_mesh_type = typename decltype(get_sampling())::discrete_dimension_type;

    /**
     * Get the domain which can be used to access the interpolation points in the sampling.
     *
     * @return domain The DDC discrete domain which maps to the sampling of interpolation points.
     */
    static ddc::DiscreteDomain<interpolation_mesh_type> get_domain()
    {
        int const npoints = ddc::discrete_space<BSplines>().ncells() + !BSplines::is_periodic();
        return ddc::DiscreteDomain<interpolation_mesh_type>(
                ddc::DiscreteElement<interpolation_mesh_type>(0),
                ddc::DiscreteVector<interpolation_mesh_type>(npoints));
    }
};
