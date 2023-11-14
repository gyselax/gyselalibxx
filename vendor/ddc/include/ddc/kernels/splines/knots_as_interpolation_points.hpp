#pragma once

#include <vector>

#include <ddc/ddc.hpp>

#include "spline_boundary_conditions.hpp"

namespace ddc {
template <class BSplines, ddc::BoundCond BcXmin, ddc::BoundCond BcXmax>
class KnotsAsInterpolationPoints
{
    static_assert(BcXmin != ddc::BoundCond::GREVILLE);
    static_assert(BcXmax != ddc::BoundCond::GREVILLE);

    using tag_type = typename BSplines::tag_type;

public:
    template <typename U = BSplines>
    static auto get_sampling()
    {
        if constexpr (U::is_uniform()) {
            using Sampling = ddc::UniformPointSampling<tag_type>;
            using SamplingImpl = typename Sampling::template Impl<Kokkos::HostSpace>;
            return SamplingImpl(
                    ddc::discrete_space<BSplines>().rmin(),
                    ddc::discrete_space<BSplines>().rmax(),
                    ddc::DiscreteVector<ddc::UniformPointSampling<tag_type>>(
                            ddc::discrete_space<BSplines>().ncells() + 1));
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

    using interpolation_mesh_type = typename decltype(get_sampling())::discrete_dimension_type;

    static ddc::DiscreteDomain<interpolation_mesh_type> get_domain()
    {
        int const npoints = ddc::discrete_space<BSplines>().ncells() + !BSplines::is_periodic();
        return ddc::DiscreteDomain<interpolation_mesh_type>(
                ddc::DiscreteElement<interpolation_mesh_type>(0),
                ddc::DiscreteVector<interpolation_mesh_type>(npoints));
    }
};
} // namespace ddc
