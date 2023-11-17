#pragma once

#include <cstddef>
#include <type_traits>
#include <vector>

#include <ddc/ddc.hpp>

#include <sll/spline_boundary_conditions.hpp>

/**
 * A class which provides helper functions to initialise the Greville points from a B-Spline definition.
 *
 * @tparam BSplines The bspline class relative to which the Greville points will be calculated.
 * @tparam BcXmin The (left) boundary condition that will be used to build the splines.
 * @tparam BcXmax The (right) boundary condition that will be used to build the splines.
 */
template <class BSplines, BoundCond BcXmin, BoundCond BcXmax>
class GrevilleInterpolationPoints
{
    using tag_type = typename BSplines::tag_type;

    template <typename U = BSplines, class = std::enable_if_t<U::is_uniform()>>
    static auto uniform_greville_points()
    {
        using Sampling = ddc::UniformPointSampling<tag_type>;
        using SamplingImpl = typename Sampling::template Impl<Kokkos::HostSpace>;

        double constexpr shift = (BSplines::degree() % 2 == 0) ? 0.5 : 0.0;
        double dx
                = (ddc::discrete_space<BSplines>().rmax() - ddc::discrete_space<BSplines>().rmin())
                  / ddc::discrete_space<BSplines>().ncells();
        return SamplingImpl(
                ddc::Coordinate<tag_type>(ddc::discrete_space<BSplines>().rmin() + shift * dx),
                ddc::Coordinate<tag_type>(dx));
    }

    template <typename U = BSplines, class = std::enable_if_t<!U::is_uniform()>>
    static auto non_uniform_greville_points()
    {
        using Sampling = ddc::NonUniformPointSampling<tag_type>;
        using SamplingImpl = typename Sampling::template Impl<Kokkos::HostSpace>;

        int n_greville_points = ddc::discrete_space<BSplines>().nbasis();
        if constexpr (U::is_periodic()) {
            n_greville_points += 1;
        }

        std::vector<double> greville_points(n_greville_points);
        ddc::DiscreteDomain<BSplines> bspline_domain
                = ddc::discrete_space<BSplines>().full_domain().take_first(
                        ddc::DiscreteVector<BSplines>(ddc::discrete_space<BSplines>().nbasis()));

        ddc::for_each(bspline_domain, [&](ddc::DiscreteElement<BSplines> ib) {
            // Define the Greville points from the bspline knots
            greville_points[ib.uid()] = 0.0;
            for (std::size_t i(0); i < BSplines::degree(); ++i) {
                greville_points[ib.uid()]
                        += ddc::discrete_space<BSplines>().get_support_knot_n(ib, i + 1);
            }
            greville_points[ib.uid()] /= BSplines::degree();
        });

        std::vector<double> temp_knots(BSplines::degree());
        // Use periodicity to ensure all points are in the domain
        if constexpr (U::is_periodic()) {
            int npoints(0);
            // Count the number of interpolation points that need shifting to preserve the ordering
            while (greville_points[npoints] < ddc::discrete_space<BSplines>().rmin()) {
                temp_knots[npoints]
                        = greville_points[npoints] + ddc::discrete_space<BSplines>().length();
                npoints++;
            }
            // Shift the points
            for (std::size_t i = 0; i < ddc::discrete_space<BSplines>().nbasis() - npoints; ++i) {
                greville_points[i] = greville_points[i + npoints];
            }
            for (int i = 0; i < npoints; ++i) {
                greville_points[ddc::discrete_space<BSplines>().nbasis() - npoints + i]
                        = temp_knots[i];
            }

            // Save a periodic point to initialise the domain size
            greville_points[n_greville_points - 1]
                    = greville_points[0] + ddc::discrete_space<BSplines>().length();
        }

        return SamplingImpl(greville_points);
    }

    static constexpr int N_BE_MIN = n_boundary_equations(BcXmin, BSplines::degree());
    static constexpr int N_BE_MAX = n_boundary_equations(BcXmax, BSplines::degree());
    template <class U>
    static constexpr bool is_uniform_mesh_v
            = U::is_uniform() && ((N_BE_MIN != 0 && N_BE_MAX != 0) || U::is_periodic());

public:
    /**
     * Get the UniformPointSampling defining the Greville points.
     *
     * This function is called when the result is a UniformPointSampling. This is the case
     * when uniform splines are used with an odd degree and with boundary conditions which
     * do not introduce additional interpolation points.
     *
     * @returns The mesh of uniform Greville points.
     */
    template <
            typename U = BSplines,
            std::enable_if_t<
                    is_uniform_mesh_v<U>,
                    bool> = true> // U must be in condition for SFINAE
    static auto get_sampling()
    {
        return uniform_greville_points();
    }

    /**
     * Get the NonUniformPointSampling defining the Greville points.
     *
     * @returns The mesh of non-uniform Greville points.
     */
    template <
            typename U = BSplines,
            std::enable_if_t<
                    !is_uniform_mesh_v<U>,
                    bool> = true> // U must be in condition for SFINAE
    static auto get_sampling()
    {
        using Sampling = ddc::NonUniformPointSampling<tag_type>;
        using SamplingImpl = typename Sampling::template Impl<Kokkos::HostSpace>;
        if constexpr (U::is_uniform()) {
            auto points_wo_bcs = uniform_greville_points();
            int const n_break_points = ddc::discrete_space<BSplines>().ncells() + 1;
            int const npoints = ddc::discrete_space<BSplines>().nbasis() - N_BE_MIN - N_BE_MAX;
            std::vector<double> points_with_bcs(npoints);

            // Construct Greville-like points at the edge
            if constexpr (BcXmin == BoundCond::GREVILLE) {
                for (std::size_t i(0); i < BSplines::degree() / 2 + 1; ++i) {
                    points_with_bcs[i]
                            = (BSplines::degree() - i) * ddc::discrete_space<BSplines>().rmin();
                    for (std::size_t j(0); j < i; ++j) {
                        points_with_bcs[i] += ddc::discrete_space<BSplines>().get_support_knot_n(
                                ddc::DiscreteElement<BSplines>(i),
                                BSplines::degree() - j);
                    }
                    points_with_bcs[i] /= BSplines::degree();
                }
            } else {
                points_with_bcs[0] = points_wo_bcs.coordinate(
                        ddc::DiscreteElement<ddc::UniformPointSampling<tag_type>>(0));
            }

            int const n_start = (BcXmin == BoundCond::GREVILLE) ? BSplines::degree() / 2 + 1 : 1;
            int const domain_size = n_break_points - 2;
            ddc::DiscreteDomain<ddc::UniformPointSampling<tag_type>> const
                    domain(ddc::DiscreteElement<ddc::UniformPointSampling<tag_type>>(1),
                           ddc::DiscreteVector<ddc::UniformPointSampling<tag_type>>(domain_size));

            // Copy central points
            ddc::for_each(domain, [&](auto ip) {
                points_with_bcs[ip.uid() + n_start - 1] = points_wo_bcs.coordinate(ip);
            });

            // Construct Greville-like points at the edge
            if constexpr (BcXmax == BoundCond::GREVILLE) {
                for (std::size_t i(0); i < BSplines::degree() / 2 + 1; ++i) {
                    points_with_bcs[npoints - 1 - i]
                            = (BSplines::degree() - i) * ddc::discrete_space<BSplines>().rmax();
                    for (std::size_t j(0); j < i; ++j) {
                        points_with_bcs[npoints - 1 - i]
                                += ddc::discrete_space<BSplines>().get_support_knot_n(
                                        ddc::DiscreteElement<BSplines>(
                                                ddc::discrete_space<BSplines>().nbasis() - 1 - i),
                                        j + 1);
                    }
                    points_with_bcs[npoints - 1 - i] /= BSplines::degree();
                }
            } else {
                points_with_bcs[npoints - 1] = points_wo_bcs.coordinate(
                        ddc::DiscreteElement<ddc::UniformPointSampling<tag_type>>(
                                ddc::discrete_space<BSplines>().ncells() - 1
                                + BSplines::degree() % 2));
            }
            return SamplingImpl(points_with_bcs);
        } else {
            auto points_wo_bcs = non_uniform_greville_points();
            if constexpr (N_BE_MIN == 0 && N_BE_MAX == 0) {
                return points_wo_bcs;
            } else {
                // All points are Greville points. Extract unnecessary points near the boundary
                std::vector<double> points_with_bcs(points_wo_bcs.size() - N_BE_MIN - N_BE_MAX);
                int constexpr n_start = N_BE_MIN;

                using length = ddc::DiscreteVector<ddc::NonUniformPointSampling<tag_type>>;

                ddc::DiscreteDomain<ddc::NonUniformPointSampling<tag_type>> const
                        domain(ddc::DiscreteElement<ddc::NonUniformPointSampling<tag_type>>(
                                       n_start),
                               length(points_with_bcs.size()));

                points_with_bcs[0] = points_wo_bcs.coordinate(domain.front());
                ddc::for_each(domain.remove(length(1), length(1)), [&](auto ip) {
                    points_with_bcs[ip.uid() - n_start] = points_wo_bcs.coordinate(ip);
                });
                points_with_bcs[points_with_bcs.size() - 1]
                        = points_wo_bcs.coordinate(domain.back());

                return SamplingImpl(points_with_bcs);
            }
        }
    }

    /**
     * The type of the mesh.
     *
     * This is either NonUniformPointSampling or UniformPointSampling.
     */
    using interpolation_mesh_type = typename decltype(get_sampling())::discrete_dimension_type;

    /**
     * Get the domain which gives us access to all of the Greville points.
     *
     * @returns The domain of the Greville points.
     */
    static ddc::DiscreteDomain<interpolation_mesh_type> get_domain()
    {
        int const npoints = ddc::discrete_space<BSplines>().nbasis() - N_BE_MIN - N_BE_MAX;
        return ddc::DiscreteDomain<interpolation_mesh_type>(
                ddc::DiscreteElement<interpolation_mesh_type>(0),
                ddc::DiscreteVector<interpolation_mesh_type>(npoints));
    }
};
