// SPDX-License-Identifier: MIT

#pragma once

#include <ddc/ddc.hpp>

#include "coord_transformation_tools.hpp"
#include "ddc_aliases.hpp"
#include "indexed_tensor.hpp"
#include "inverse_jacobian_matrix.hpp"
#include "l_norm_tools.hpp"
#include "multipatch_type.hpp"
#include "tensor.hpp"
#include "types.hpp"
#include "vector_index_tools.hpp"


/**
 * @brief Patch locator using a Newton algorithm to invert the coordinate transformation
 * of each patch.
 *
 * This operator locates on which patch a given physical coordinate is. No assumption is made
 * on the shape of the patches (convex or not). The coordinate transformations are only assumed
 * to be continuous and with a non-singular Jacobian inside the patches.
 *
 * @tparam MultipatchMapping A MultipatchType type containing the coordinate
 *          transformations (logical -> physical) on each patch.
 * @tparam ExecSpace The space (CPU/GPU) where the coordinate transformations can be evaluated.
 */
template <class MultipatchMapping, class ExecSpace = Kokkos::DefaultExecutionSpace>
class PatchLocatorNewton;


/**
 * @brief Patch locator using a Newton algorithm to invert the coordinate transformation
 * of each patch.
 *
 * @tparam Patches The patches.
 * @tparam Mappings The coordinate transformations (logical -> physical) on each patch.
 * @tparam ExecSpace The space (CPU/GPU) where the coordinate transformations can be evaluated.
 */
template <class... Patches, class... Mappings, class ExecSpace>
class PatchLocatorNewton<
        detail::MultipatchType<ddc::TypeSeq<Patches...>, ddc::TypeSeq<Mappings...>>,
        ExecSpace>
{
public:
    /// @brief Sequence ddc::TypeSeq of patch tags.
    using patch_seq = ddc::TypeSeq<Patches...>;

    /// @brief The space (CPU/GPU) where the coordinate transformations can be evaluated.
    using exec_space = ExecSpace;

    /// @brief Default value to define outside domain (not a patch).
    static constexpr int outside_domain = -1;

private:
    using mapping_on_first_patch = ddc::type_seq_element_t<0, ddc::TypeSeq<Mappings...>>;

    using ArgBasis = ddc::to_type_seq_t<typename mapping_on_first_patch::CoordArg>;
    using ResBasis = ddc::to_type_seq_t<typename mapping_on_first_patch::CoordResult>;

    static_assert(
            (std::is_same_v<ArgBasis, ddc::to_type_seq_t<typename Mappings::CoordArg>> && ...),
            "The argument continuous dimensions of all the coordinate transformations should be "
            "the same type on all the patches.");

    static_assert(
            (std::is_same_v<ResBasis, ddc::to_type_seq_t<typename Mappings::CoordResult>> && ...),
            "The result continuous dimensions of all the coordinate transformations should be "
            "the same type on all the patches.");

    using DimArg1 = ddc::type_seq_element_t<0, ArgBasis>;
    using DimArg2 = ddc::type_seq_element_t<1, ArgBasis>;

    using DimRes1 = ddc::type_seq_element_t<0, ResBasis>;
    using DimRes2 = ddc::type_seq_element_t<1, ResBasis>;

    using CoordArg = Coord<DimArg1, DimArg2>;
    using CoordRes = Coord<DimRes1, DimRes2>;

    using MultipatchMapping
            = detail::MultipatchType<ddc::TypeSeq<Patches...>, ddc::TypeSeq<Mappings...>>;

    /// @brief Type of the coordinate transformation on the patch P.
    template <class P>
    using MappingOnPatch = typename MultipatchMapping::template TypeOnPatch<P>;

private:
    static constexpr std::size_t n_patches = ddc::type_seq_size_v<patch_seq>;

    /// @brief Number of starting points tried on a patch before rejecting it.
    static constexpr int n_starts = 3;

    /// @brief Number of samples per logical dimension used to choose the starting points.
    static constexpr int n_coarse_samples = 5;

    MultipatchMapping m_coord_transformations;

    std::size_t m_counter_max;
    double m_TOL;

    /// Bounding box {xmin, xmax, ymin, ymax} of each patch in the physical domain.
    Kokkos::Array<Kokkos::Array<double, 4>, n_patches> m_boxes;
    /// Logical domain {dim1_min, dim1_max, dim2_min, dim2_max} of each patch.
    Kokkos::Array<Kokkos::Array<double, 4>, n_patches> m_logical_bounds;

public:
    /**
     * @brief Instantiate the patch locator solving with a Newton algorithm.
     *
     * @param coord_transformations A MultipatchType of coordinate transformations from the logical
     *          domain to the physical domain of each patch.
     * @param counter_max Maximum number of loops for the Newton algorithm.
     * @param TOL Tolerance on the error for the Newton algorithm.
     * @param n_samples_per_edge Number of samples on each edge of the patches used to compute
     *          their bounding boxes.
     */
    explicit PatchLocatorNewton(
            MultipatchMapping const& coord_transformations,
            std::size_t const counter_max = 20,
            double const TOL = 1e-12,
            int const n_samples_per_edge = 64)
        : m_coord_transformations(coord_transformations)
        , m_counter_max(counter_max)
        , m_TOL(TOL)
    {
        set_boxes(n_samples_per_edge);
    }

    /**
     * @brief Get the patch where the given physical coordinate is.
     *
     * @param[in] coord The given physical coordinate.
     * @return [int] The patch index where the physical coordinate. If the coordinate
     *              is outside of the domain, it returns a negative value.
     */
    KOKKOS_INLINE_FUNCTION int operator()(CoordRes const coord) const
    {
        CoordArg coord_log = no_guess();
        return find_patch<patch_seq>(coord_log, coord);
    }

    /**
     * @brief Get the patch where the given physical coordinate is.
     *
     * @param[in, out] coord_log The logical coordinate corresponding to the given
     *              physical coordinate. The given value is used as initial condition
     *              in the Newton algorithm (NaN if no initial guess is known).
     * @param[in] coord The given physical coordinate.
     * @return [int] The patch index where the physical coordinate. If the coordinate
     *              is outside of the domain, it returns a negative value.
     */
    KOKKOS_INLINE_FUNCTION int operator()(CoordArg& coord_log, CoordRes const coord) const
    {
        return find_patch<patch_seq>(coord_log, coord);
    }

    /**
     * @brief Get the patch where the given physical coordinate is, testing first the patch
     * PatchInit.
     *
     * @tparam PatchInit The first patch tested.
     * @param[in] coord The given physical coordinate.
     * @return [int] The patch index where the physical coordinate. If the coordinate
     *              is outside of the domain, it returns a negative value.
     */
    template <class PatchInit>
    KOKKOS_INLINE_FUNCTION int find(CoordRes const coord) const
    {
        CoordArg coord_log = no_guess();
        using patch_seq_without_init = ddc::type_seq_remove_t<patch_seq, ddc::TypeSeq<PatchInit>>;
        using patch_seq_starting_with_init
                = ddc::type_seq_merge_t<ddc::TypeSeq<PatchInit>, patch_seq_without_init>;
        return find<patch_seq_starting_with_init>(coord_log, coord);
    }

    /**
     * @brief Get the patch where the given physical coordinate is, testing first the patch
     * PatchInit.
     *
     * @tparam PatchInit The first patch tested.
     * @param[in, out] coord_log The logical coordinate corresponding to the given
     *              physical coordinate. The given value is used as initial condition
     *              in the Newton algorithm (NaN if no initial guess is known).
     * @param[in] coord The given physical coordinate.
     * @return [int] The patch index where the physical coordinate. If the coordinate
     *              is outside of the domain, it returns a negative value.
     */
    template <class PatchInit>
    KOKKOS_INLINE_FUNCTION int find(CoordArg& coord_log, CoordRes const coord) const
    {
        // type_seq_remove_t<A, B> removes from A the elements of B.
        using patch_seq_without_init = ddc::type_seq_remove_t<patch_seq, ddc::TypeSeq<PatchInit>>;
        using patch_seq_starting_with_init
                = ddc::type_seq_merge_t<ddc::TypeSeq<PatchInit>, patch_seq_without_init>;
        return find_patch<patch_seq_starting_with_init>(coord_log, coord);
    }

public:
    /**
     * @brief Compute the bounding box of each patch by sampling the image of its boundary.
     *
     * The extrema of a bounded closed set of the plane are reached on its boundary, so this
     * works for any patch shape (convex or not). A margin equal to the largest distance
     * between two consecutive samples is added to account for the curvature of the boundary
     * between two samples.
     *
     * The coordinate transformations are evaluated on ExecSpace (they may only be accessible
     * on device, e.g. spline mappings).
     * This function is public only because CUDA forbids lambdas in private member functions.
     *
     * @param n_samples Number of samples on each edge of the patches.
     */
    template <std::size_t I = 0>
    void set_boxes(int const n_samples)
    {
        if constexpr (I < n_patches) {
            using P = ddc::type_seq_element_t<I, patch_seq>;
            double const dim1_min = ddc::discrete_space<typename P::BSplines1>().rmin();
            double const dim1_max = ddc::discrete_space<typename P::BSplines1>().rmax();
            double const dim2_min = ddc::discrete_space<typename P::BSplines2>().rmin();
            double const dim2_max = ddc::discrete_space<typename P::BSplines2>().rmax();
            m_logical_bounds[I] = {dim1_min, dim1_max, dim2_min, dim2_max};

            // Local copy so that the lambda does not capture `this` (host pointer).
            MappingOnPatch<P> const mapping = m_coord_transformations.template get<P>();
            Kokkos::Array<double, 4> const bounds = m_logical_bounds[I];

            Kokkos::MinMaxScalar<double> mm_x;
            Kokkos::MinMaxScalar<double> mm_y;
            double max_step = 0.;

            // Segment k of the boundary: 4 edges with n_samples segments each.
            Kokkos::parallel_reduce(
                    "PatchLocatorNewton::set_boxes",
                    Kokkos::RangePolicy<ExecSpace>(ExecSpace(), 0, 4 * n_samples),
                    KOKKOS_LAMBDA(
                            int const k,
                            Kokkos::MinMaxScalar<double>& lmm_x,
                            Kokkos::MinMaxScalar<double>& lmm_y,
                            double& lmax_step) {
                        int const edge = k / n_samples;
                        int const sample = k % n_samples;
                        Kokkos::Array<CoordRes, 2> const coords
                                = {mapping(boundary_point(bounds, edge, sample, n_samples)),
                                   mapping(boundary_point(bounds, edge, sample + 1, n_samples))};
                        for (int c = 0; c < 2; ++c) {
                            double const x = ddc::get<DimRes1>(coords[c]);
                            double const y = ddc::get<DimRes2>(coords[c]);
                            lmm_x.min_val = Kokkos::min(lmm_x.min_val, x);
                            lmm_x.max_val = Kokkos::max(lmm_x.max_val, x);
                            lmm_y.min_val = Kokkos::min(lmm_y.min_val, y);
                            lmm_y.max_val = Kokkos::max(lmm_y.max_val, y);
                        }
                        lmax_step = Kokkos::max(lmax_step, norm_inf(coords[1] - coords[0]));
                    },
                    Kokkos::MinMax<double>(mm_x),
                    Kokkos::MinMax<double>(mm_y),
                    Kokkos::Max<double>(max_step));

            // Margin: the image of the boundary can bulge between two samples.
            m_boxes[I]
                    = {mm_x.min_val - max_step,
                       mm_x.max_val + max_step,
                       mm_y.min_val - max_step,
                       mm_y.max_val + max_step};

            set_boxes<I + 1>(n_samples);
        }
    }

private:
    /**
     * @brief Logical coordinate of the sample t_idx (in [0, n_samples]) on an edge of the
     * logical domain {dim1_min, dim1_max, dim2_min, dim2_max}.
     * Edges:
     *      0 = {dim2 = dim2_min}, 
     *      1 = {dim2 = dim2_max}, 
     *      2 = {dim1 = dim1_min},
     *      3 = {dim1 = dim1_max}.
     */
    KOKKOS_INLINE_FUNCTION static CoordArg boundary_point(
            Kokkos::Array<double, 4> const& bounds,
            int const edge,
            int const t_idx,
            int const n_samples)
    {
        double const t = double(t_idx) / n_samples;
        double const dim1 = (edge == 2)   ? bounds[0]
                            : (edge == 3) ? bounds[1]
                                          : bounds[0] + t * (bounds[1] - bounds[0]);
        double const dim2 = (edge == 0)   ? bounds[2]
                            : (edge == 1) ? bounds[3]
                                          : bounds[2] + t * (bounds[3] - bounds[2]);
        return CoordArg(dim1, dim2);
    }

private:
    /// @brief Logical coordinate meaning "no initial guess".
    KOKKOS_INLINE_FUNCTION static CoordArg no_guess()
    {
        double const nan = Kokkos::Experimental::quiet_NaN_v<double>;
        return CoordArg(nan, nan);
    }

    /// @brief Index of the patch P in patch_seq.
    template <class P>
    static constexpr std::size_t patch_idx = ddc::type_seq_rank_v<P, patch_seq>;

    /**
     * @brief Find recursively on which patch the physical coordinate is located (patches
     * tested in the order of PatchSeq) and update the logical coordinate.
     */
    template <class PatchSeq, std::size_t I = 0>
    KOKKOS_INLINE_FUNCTION int find_patch(CoordArg& coord_log, CoordRes const& coord) const
    {
        if constexpr (I == n_patches) {
            return outside_domain;
        } else {
            using P = ddc::type_seq_element_t<I, PatchSeq>;
            // Each patch starts from the same initial guess (not from the end point of a
            // failed Newton on the previous patch).
            CoordArg guess = coord_log;
            if (is_in_patch<P>(guess, coord)) {
                coord_log = guess;
                // Index in patch_seq (not in PatchSeq, which may be reordered).
                return static_cast<int>(patch_idx<P>);
            }
            return find_patch<PatchSeq, I + 1>(coord_log, coord);
        }
    }

    /**
     * @brief Check if the physical coordinate is located inside the patch P.
     *
     * 1. Fast rejection with the bounding box of the patch.
     * 2. Damped Newton from the given initial guess (if not NaN).
     * 3. If it fails, damped Newton from the n_starts closest points of a coarse sampling
     *    of the patch (robustness for non-convex patches).
     *
     * @param[in, out] logical_coord Initial guess, and logical coordinate of the physical
     *          coordinate if it is found on the patch.
     */
    template <class P>
    KOKKOS_INLINE_FUNCTION bool is_in_patch(CoordArg& logical_coord, CoordRes const& physical_coord)
            const
    {
        // Quick check to elimate patch if it is not inside the smallest box containing the patch.
        if (!is_in_box<P>(physical_coord)) {
            return false;
        }

        // Fix initial condition if is NaN.
        if (!Kokkos::isnan(ddc::get<DimArg1>(logical_coord))
            && !Kokkos::isnan(ddc::get<DimArg2>(logical_coord))) {
            CoordArg guess = logical_coord;
            if (damped_newton<P>(guess, physical_coord)) {
                logical_coord = guess;
                return true;
            }
        }

        // Apply a damped Newton algorithm.
        Kokkos::Array<CoordArg, n_starts> const starts = starting_points<P>(physical_coord);
        for (int s = 0; s < n_starts; ++s) {
            CoordArg guess = starts[s];
            if (damped_newton<P>(guess, physical_coord)) {
                logical_coord = guess;
                return true;
            }
        }
        return false;
    }

    /// @brief Check if the physical coordinate is inside the bounding box of the patch P.
    template <class P>
    KOKKOS_INLINE_FUNCTION bool is_in_box(CoordRes const& coord) const
    {
        Kokkos::Array<double, 4> const& box = m_boxes[patch_idx<P>];
        double const x = ddc::get<DimRes1>(coord);
        double const y = ddc::get<DimRes2>(coord);
        return (box[0] <= x) && (x <= box[1]) && (box[2] <= y) && (y <= box[3]);
    }

    /**
     * @brief Damped Newton algorithm to invert the coordinate transformation of the patch P.
     *
     * If a Newton step does not decrease the error, the step is halved (up to a few times).
     * If the error cannot be decreased anymore, the algorithm stops: the point is either
     * outside the patch or the Newton is stuck.
     *
     * @return True if the error is lower than the tolerance.
     */
    template <class P>
    KOKKOS_INLINE_FUNCTION bool damped_newton(
            CoordArg& logical_coord,
            CoordRes const& physical_coord) const
    {
        constexpr int max_halving = 6;
        MappingOnPatch<P> const mapping = m_coord_transformations.template get<P>();
        InverseJacobianMatrix const get_inverse_jacobian(mapping);

        replace_inside<P>(logical_coord);
        CoordRes computed_physical_coord = mapping(logical_coord);
        double error = norm_inf(computed_physical_coord - physical_coord);

        for (std::size_t counter = 0; (counter < m_counter_max) && !(error <= m_TOL); ++counter) {
            DVector<DimRes1, DimRes2> const residual(computed_physical_coord - physical_coord);
            DVector<DimArg1, DimArg2> const delta = tensor_mul(
                    index<'i', 'j'>(get_inverse_jacobian(logical_coord)),
                    index<'j'>(residual));
            double const delta1 = ddcHelper::get<DimArg1>(delta);
            double const delta2 = ddcHelper::get<DimArg2>(delta);

            double step = 1.;
            CoordArg trial = logical_coord;
            CoordRes trial_physical_coord = computed_physical_coord;
            double trial_error = error;
            for (int h = 0; h <= max_halving; ++h) {
                trial = logical_coord - CoordArg(step * delta1, step * delta2);
                replace_inside<P>(trial);
                trial_physical_coord = mapping(trial);
                trial_error = norm_inf(trial_physical_coord - physical_coord);
                if (trial_error < error) {
                    break;
                }
                step *= 0.5;
            }
            // No decrease (also if NaN): stuck against a bound or outside the patch.
            if (!(trial_error < error)) {
                break;
            }
            logical_coord = trial;
            computed_physical_coord = trial_physical_coord;
            error = trial_error;
        }
        return error <= m_TOL;
    }

    /**
     * @brief Get the n_starts logical coordinates of a coarse uniform sampling of the patch P
     * whose images are the closest to the physical coordinate, each one refined by a
     * coarse-to-fine search.
     */
    template <class P>
    KOKKOS_INLINE_FUNCTION Kokkos::Array<CoordArg, n_starts> starting_points(
            CoordRes const& physical_coord) const
    {
        MappingOnPatch<P> const mapping = m_coord_transformations.template get<P>();
        Kokkos::Array<double, 4> const& bounds = m_logical_bounds[patch_idx<P>];
        double const dim1_min = bounds[0];
        double const dim2_min = bounds[2];
        double const step1 = (bounds[1] - dim1_min) / (n_coarse_samples - 1);
        double const step2 = (bounds[3] - dim2_min) / (n_coarse_samples - 1);

        // Keep the n_starts closest samples, sorted by increasing distance.
        Kokkos::Array<CoordArg, n_starts> best;
        Kokkos::Array<double, n_starts> best_dist;
        for (int s = 0; s < n_starts; ++s) {
            best[s] = CoordArg(dim1_min, dim2_min);
            best_dist[s] = Kokkos::Experimental::infinity_v<double>;
        }
        for (int i = 0; i < n_coarse_samples; ++i) {
            for (int j = 0; j < n_coarse_samples; ++j) {
                CoordArg const sample(dim1_min + i * step1, dim2_min + j * step2);
                double const dist = norm_inf(mapping(sample) - physical_coord);
                // Insertion in the sorted list.
                int pos = n_starts;
                while ((pos > 0) && (dist < best_dist[pos - 1])) {
                    if (pos < n_starts) {
                        best[pos] = best[pos - 1];
                        best_dist[pos] = best_dist[pos - 1];
                    }
                    --pos;
                }
                if (pos < n_starts) {
                    best[pos] = sample;
                    best_dist[pos] = dist;
                }
            }
        }

        for (int s = 0; s < n_starts; ++s) {
            best[s] = refine<P>(best[s], best_dist[s], 0.5 * step1, 0.5 * step2, physical_coord);
        }
        return best;
    }

    /**
     * @brief Coarse-to-fine search of the closest point around a starting point.
     *
     * At each level, the mapping is evaluated on a 3x3 stencil around the current point
     * and the stencil is recentred on the closest point and halved.
     */
    template <class P>
    KOKKOS_INLINE_FUNCTION CoordArg
    refine(CoordArg const& start,
           double best_dist,
           double delta_dim1,
           double delta_dim2,
           CoordRes const& physical_coord) const
    {
        constexpr int n_levels = 4;
        MappingOnPatch<P> const mapping = m_coord_transformations.template get<P>();
        Kokkos::Array<double, 4> const& bounds = m_logical_bounds[patch_idx<P>];

        double dim1_mid = ddc::get<DimArg1>(start);
        double dim2_mid = ddc::get<DimArg2>(start);
        for (int level = 0; level < n_levels; ++level) {
            double best1 = dim1_mid;
            double best2 = dim2_mid;
            for (int i = -1; i <= 1; ++i) {
                for (int j = -1; j <= 1; ++j) {
                    if ((i == 0) && (j == 0)) {
                        continue; // Already evaluated.
                    }
                    double const dim1
                            = Kokkos::clamp(dim1_mid + i * delta_dim1, bounds[0], bounds[1]);
                    double const dim2
                            = Kokkos::clamp(dim2_mid + j * delta_dim2, bounds[2], bounds[3]);
                    double const dist = norm_inf(mapping(CoordArg(dim1, dim2)) - physical_coord);
                    if (dist < best_dist) {
                        best_dist = dist;
                        best1 = dim1;
                        best2 = dim2;
                    }
                }
            }
            dim1_mid = best1;
            dim2_mid = best2;
            delta_dim1 *= 0.5;
            delta_dim2 *= 0.5;
        }
        return CoordArg(dim1_mid, dim2_mid);
    }

    /// @brief Replace the logical coordinate inside the logical domain of the patch P.
    template <class P>
    KOKKOS_INLINE_FUNCTION void replace_inside(CoordArg& coord_log) const
    {
        Kokkos::Array<double, 4> const& bounds = m_logical_bounds[patch_idx<P>];
        double const dim1 = put_inside<
                P::BSplines1::is_periodic()>(ddc::get<DimArg1>(coord_log), bounds[0], bounds[1]);
        double const dim2 = put_inside<
                P::BSplines2::is_periodic()>(ddc::get<DimArg2>(coord_log), bounds[2], bounds[3]);
        coord_log = CoordArg(dim1, dim2);
    }

    /**
     * @brief Put a 1D logical coordinate inside [xmin, xmax].
     * Periodic dimension: the coordinate is wrapped (any number of periods).
     * Non-periodic dimension: the coordinate is clamped to the closest bound.
     */
    template <bool Periodic>
    KOKKOS_INLINE_FUNCTION static double put_inside(
            double const x,
            double const xmin,
            double const xmax)
    {
        if constexpr (Periodic) {
            double const L = xmax - xmin;
            return x - L * Kokkos::floor((x - xmin) / L);
        } else {
            return Kokkos::clamp(x, xmin, xmax);
        }
    }
};


// To help the template deduction.
template <class MultipatchMapping>
PatchLocatorNewton(MultipatchMapping const& coord_transformations)
        -> PatchLocatorNewton<MultipatchMapping>;

template <class MultipatchMapping>
PatchLocatorNewton(
        MultipatchMapping const& coord_transformations,
        std::size_t const counter_max,
        double const TOL) -> PatchLocatorNewton<MultipatchMapping>;
