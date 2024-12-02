// SPDX-License-Identifier: MIT

#pragma once

#include <cassert>
#include <cmath>

#include <ddc/ddc.hpp>
#include <ddc/kernels/splines.hpp>

#include "ddc_aliases.hpp"

namespace ddcHelper {

template <class BSplines, ddc::BoundCond BcXmin, ddc::BoundCond BcXmax>
class NonUniformInterpolationPoints;

template <class T>
struct is_non_uniform_interpolation_points : std::false_type
{
};

template <class BSplines, ddc::BoundCond BcXmin, ddc::BoundCond BcXmax>
struct is_non_uniform_interpolation_points<NonUniformInterpolationPoints<BSplines, BcXmin, BcXmax>>
    : std::true_type
{
};

template <class T>
inline constexpr bool is_non_uniform_interpolation_points_v
        = is_non_uniform_interpolation_points<T>::value;

/**
 * @brief Helper class for the initialisation of the mesh of interpolation points.
 *
 * A helper class for the initialisation of the mesh of interpolation points. This
 * class should be used when the interpolation points are specified externally. 
 * This is possible with any kind of boundary condition. In the case of Greville 
 * boundary conditions, the second and second to last interpolation points 
 * should be the Greville points. 
 */
template <class BSplines, ddc::BoundCond BcXmin, ddc::BoundCond BcXmax>
class NonUniformInterpolationPoints
{
    using Dim = typename BSplines::continuous_dimension_type;
    using BreakPointsType = std::conditional_t<
            BSplines::is_uniform(),
            ddc::UniformBsplinesKnots<BSplines>,
            ddc::NonUniformBsplinesKnots<BSplines>>;
    using IdxBreakPoints = Idx<BreakPointsType>;

public:
    /// The number of boundary equations at the lower bound in a spline interpolation.
    static constexpr int N_BE_MIN = n_boundary_equations(BcXmin, BSplines::degree());
    /// The number of boundary equations at the upper bound in a spline interpolation.
    static constexpr int N_BE_MAX = n_boundary_equations(BcXmax, BSplines::degree());

private:
    static void check_n_points_in_cell(int n_points_in_cell, IdxBreakPoints current_cell_end_idx)
    {
        if (n_points_in_cell > BSplines::degree() + 1) {
            IdxBreakPoints rmin_idx = ddc::discrete_space<BSplines>().break_point_domain().front();
            int failed_cell = (current_cell_end_idx - rmin_idx).value();
            throw std::runtime_error(
                    "The spline problem is overconstrained. There are "
                    + std::to_string(n_points_in_cell) + " points in the "
                    + std::to_string(failed_cell) + "-th cell.");
        }
    }

public:
    /**
     * Get the sampling of interpolation points.
     *
     * @param[in] interp_points The grid points on which the simulation will evolve.
     * @param[in] TOL The tolerance above which points are considered to be non-uniform.
     *
     * @return sampling The DDC point sampling of the interpolation points.
     */
    template <typename Sampling, typename U = BSplines>
    static auto get_sampling(std::vector<Coord<Dim>>& interp_points, double TOL = 2e-14)
    {
        int const expected_npoints = ddc::discrete_space<BSplines>().nbasis() - N_BE_MIN - N_BE_MAX;
        if (interp_points.size() != expected_npoints) {
            throw std::runtime_error(
                    "Incorrect number of points supplied to NonUniformInterpolationPoints. "
                    "(Received : "
                    + std::to_string(interp_points.size())
                    + ", expected : " + std::to_string(expected_npoints));
        }
        if constexpr (BSplines::is_periodic()) {
            double domain_len = ddc::discrete_space<BSplines>().rmax()
                                - ddc::discrete_space<BSplines>().rmin();
            // Insert final periodic point to ensure all cell sizes are known
            interp_points.push_back(interp_points.back() + domain_len);
        }
        int n_points_in_cell = 0;
        IdxBreakPoints current_cell_end_idx
                = ddc::discrete_space<BSplines>().break_point_domain().front() + 1;
        for (Coord<Dim> point : interp_points) {
            if (point == ddc::coordinate(current_cell_end_idx)) {
                check_n_points_in_cell(n_points_in_cell + 1, current_cell_end_idx);
                n_points_in_cell = 1;
                current_cell_end_idx += 1;
            } else if (point > ddc::coordinate(current_cell_end_idx)) {
                check_n_points_in_cell(n_points_in_cell, current_cell_end_idx);
                n_points_in_cell = 1;
                current_cell_end_idx += 1;
            } else {
                n_points_in_cell += 1;
            }
        }
        check_n_points_in_cell(n_points_in_cell, current_cell_end_idx);

        if constexpr (ddc::is_non_uniform_point_sampling_v<Sampling>) {
            using SamplingImpl = typename Sampling::template Impl<Sampling, Kokkos::HostSpace>;
            return SamplingImpl(interp_points);
        } else {
            using SamplingImpl = typename Sampling::template Impl<Sampling, Kokkos::HostSpace>;
            SamplingImpl result(interp_points[0], interp_points[1] - interp_points[0]);
            IdxRange<Sampling> idx_range(get_domain<Sampling>());
            bool same_points = ddc::transform_reduce(
                    idx_range,
                    true,
                    ddc::reducer::land<bool>(),
                    [&](Idx<Sampling> idx) {
                        return fabs(result.coordinate(idx) - interp_points[idx - idx_range.front()])
                               < TOL;
                    });
            if (!same_points) {
                throw std::runtime_error("Provided points are not uniform");
            }
            return result;
        }
    }


    /**
     * Get the domain which can be used to access the interpolation points in the sampling.
     *
     * @return domain The discrete domain which maps to the sampling of interpolation points.
     */
    template <typename Sampling>
    static IdxRange<Sampling> get_domain()
    {
        int const npoints = ddc::discrete_space<BSplines>().nbasis() - N_BE_MIN - N_BE_MAX;
        return IdxRange<Sampling>(Idx<Sampling>(0), IdxStep<Sampling>(npoints));
    }

    /**
     * The type of the mesh.
     *
     * This is always NonUniformPointSampling.
     */
    using interpolation_discrete_dimension_type = NonUniformGridBase<Dim>;
};
} // namespace ddcHelper
