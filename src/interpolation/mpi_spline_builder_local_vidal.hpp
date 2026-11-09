// SPDX-License-Identifier: MIT
#pragma once

#include <mpi.h>

#include <array>
#include <cassert>
#include <cmath>
#include <optional>

#include <ddc/ddc.hpp>
#include <ddc/kernels/splines.hpp>

#include "ddc_alias_inline_functions.hpp"
#include "ddc_aliases.hpp"
#include "deriv_reconstruction_coeffs.hpp"
#include "geometry_descriptors.hpp"
#include "i_interpolation_builder.hpp"
#include "mpitools.hpp"

template <
        class GlobalBSplines,
        class LocalBSplinesLower,
        class LocalBSplines,
        class LocalBSplinesUpper>
void construct_local_bsplines_and_neighbours(MPI_Comm comm)
{
    using Knots = ddc::knot_discrete_dimension_t<GlobalBSplines>;

    int rank, comm_size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &comm_size);

    IdxRange<Knots> idx_range_bp = ddc::discrete_space<GlobalBSplines>().break_point_domain();

    std::size_t const global_ncells = ddc::discrete_space<GlobalBSplines>().ncells();

    if (global_ncells % comm_size != 0) {
        std::ostringstream error_msg;
        error_msg << "The provided BSpline cannot be split equally over the specified number "
                     "of MPI ranks ("
                  << global_ncells << " is not divisible by " << comm_size << ")";
        throw std::runtime_error(error_msg.str());
    }

    IdxStep<Knots> const local_ncells(global_ncells / comm_size);

    bool constexpr is_periodic = GlobalBSplines::is_periodic();

    if constexpr (ddc::is_uniform_bsplines_v<GlobalBSplines>) {
        if (is_periodic or rank != 0) {
            int offset = rank == 0 ? comm_size - 1 : rank - 1;
            Idx<Knots> idx_bp_front = idx_range_bp.front() + offset * local_ncells;
            ddc::init_discrete_space<LocalBSplinesLower>(
                    ddc::coordinate(idx_bp_front),
                    ddc::coordinate(idx_bp_front + local_ncells),
                    local_ncells);
        }
        ddc::init_discrete_space<LocalBSplines>(
                ddc::coordinate(idx_range_bp.front() + rank * local_ncells),
                ddc::coordinate(idx_range_bp.front() + (rank + 1) * local_ncells),
                local_ncells);
        if (is_periodic or rank != comm_size - 1) {
            Idx<Knots> idx_bp_front
                    = idx_range_bp.front()
                      + static_cast<int>(rank != comm_size - 1) * (rank + 1) * local_ncells;
            ddc::init_discrete_space<LocalBSplinesUpper>(
                    ddc::coordinate(idx_bp_front),
                    ddc::coordinate(idx_bp_front + local_ncells),
                    local_ncells);
        }
    } else {
        if (is_periodic or rank != 0) {
            int offset = rank == 0 ? comm_size - 1 : rank - 1;
            Idx<Knots> idx_bp_front = idx_range_bp.front() + offset * local_ncells;
            ddc::init_discrete_space<LocalBSplinesLower>(idx_bp_front, idx_bp_front + local_ncells);
        }
        ddc::init_discrete_space<LocalBSplines>(
                idx_range_bp.begin() + rank * local_ncells,
                idx_range_bp.begin() + (rank + 1) * local_ncells);
        if (is_periodic or rank != comm_size - 1) {
            Idx<Knots> idx_bp_front = idx_range_bp.front()
                                      + static_cast<int>(!is_periodic) * (rank + 1) * local_ncells;
            ddc::init_discrete_space<LocalBSplinesUpper>(idx_bp_front, idx_bp_front + local_ncells);
        }
    }
}

/**
 * @brief A builder which computes, on the local subdomain owned by one MPI rank,
 * a cubic spline approximation of a function whose interpolation mesh is split
 * across MPI ranks.
 *
 * A plain ddc::SplineBuilder cannot be applied independently on each rank's
 * subdomain because a cubic spline needs a Hermite (derivative) boundary
 * condition at each end, and the correct value of that derivative at a rank
 * boundary depends on data owned by the neighbouring rank.
 *
 * TODO:
 * This class approximates that boundary derivative locally, using a fixed-width
 * linear combination of the 10 nearest function values on each side of the
 * boundary (see ):
 * @f[
 *      s'(x_i) = \sum_{j=-10}^{10} w_{j}\, f_{i+j}, \qquad w_{-j} = -w_j,\ w_0 = 0
 * @f]
 * where @f$x_i@f$ is the point shared by the two ranks. Only the immediate left/right
 * MPI neighbour is involved (rank-1 / rank+1 in the communicator, or the wrap-around
 * rank at the two ends when the boundary condition is PERIODIC).
 *
 * @tparam ExecSpace The Kokkos execution space on which the spline approximation is performed.
 * @tparam MemorySpace The Kokkos memory space on which the data is stored.
 * @tparam BSplines The B-spline basis. Must be cubic (degree 3) and uniform.
 * @tparam InterpolationGrid1D The local interpolation mesh.
 * @tparam BCLower The boundary condition applied at the lower bound of the *global*
 * domain (only relevant to the rank which owns that bound). Must be
 * @c ddc::SplineBuilderClosure::PERIODIC, @c HERMITE or @c HOMOGENEOUS_HERMITE.
 * @tparam BCUpper The boundary condition applied at the upper bound of the *global*
 * domain (only relevant to the rank which owns that bound). Same restriction as @c BCLower.
 * @tparam Solver The linear solver used internally by the wrapped @c ddc::SplineBuilder.
 */
template <
        class ExecSpace,
        class MemorySpace,
        class LocalBSplinesLower,
        class LocalBSplines,
        class LocalBSplinesUpper,
        class BSplines,
        class InterpolationGrid1D,
        ddc::SplineBuilderClosure BCLower,
        ddc::SplineBuilderClosure BCUpper,
        ddc::SplineSolver Solver = ddc::SplineSolver::LAPACK>
class MPISplineBuilderLocalVidal
{
    static_assert(
            BSplines::degree() == 3,
            "MPISplineBuilderLocalVidal only supports cubic (degree 3) B-splines.");
    static_assert(
            BSplines::is_uniform(),
            "MPISplineBuilderLocalVidal only supports uniform B-splines.");

    /// @brief Check that a boundary condition is one of the values this class supports.
    static constexpr bool is_valid_spline_bc(ddc::SplineBuilderClosure bc)
    {
        return (bc == ddc::SplineBuilderClosure::PERIODIC)
               || (bc == ddc::SplineBuilderClosure::HERMITE)
               || (bc == ddc::SplineBuilderClosure::HOMOGENEOUS_HERMITE);
    }

    static_assert(
            is_valid_spline_bc(BCLower),
            "BCLower must be PERIODIC, HERMITE or HOMOGENEOUS_HERMITE.");
    static_assert(
            is_valid_spline_bc(BCUpper),
            "BCUpper must be PERIODIC, HERMITE or HOMOGENEOUS_HERMITE.");
    static_assert(
            (BCLower == ddc::SplineBuilderClosure::PERIODIC)
                    == (BCUpper == ddc::SplineBuilderClosure::PERIODIC),
            "PERIODIC must be specified on both bounds or neither.");

public:
    /// @brief The type of the Kokkos execution space used by this class.
    using exec_space = ExecSpace;

    /// @brief The type of the Kokkos memory space used by this class.
    using memory_space = MemorySpace;

    /// @brief The data type that the data is saved on.
    using data_type = double;

    /// @brief The discrete dimension on which interpolation points are defined.
    using interpolation_grid_type = InterpolationGrid1D;

    /// @brief The continuous dimension of interest.
    using continuous_dimension_type = typename interpolation_grid_type::continuous_dimension_type;

    /// @brief The type of the 1D index range of the local interpolation mesh.
    using interpolation_idx_range_type = IdxRange<InterpolationGrid1D>;

    /// @brief The type of the Deriv dimension at the rank/domain boundaries.
    using deriv_type = ddc::Deriv<continuous_dimension_type>;

    /// @brief The index range for the interpolation coefficients (1D representative).
    using coeff_idx_range_type = IdxRange<BSplines>;

    /// @brief The batched index range type with interpolation_grid_type replaced by the B-spline basis.
    template <class IdxRangeBatchedInterpolation>
    using batched_basis_idx_range_type = ddc::
            replace_dim_of_t<IdxRangeBatchedInterpolation, interpolation_grid_type, LocalBSplines>;

    /// @brief The batched index range type with interpolation_grid_type replaced by deriv_type.
    template <class IdxRangeBatchedInterpolation>
    using batched_derivs_idx_range_type = ddc::
            replace_dim_of_t<IdxRangeBatchedInterpolation, interpolation_grid_type, deriv_type>;

private:
    /// @brief The internal, per-rank spline builder. Always uses HERMITE on both ends:
    /// whichever boundary (true global edge or interior rank boundary) is being handled,
    /// this class always supplies an explicit derivative value for it.
    using internal_builder_type = ddc::SplineBuilder<
            ExecSpace,
            MemorySpace,
            LocalBSplines,
            interpolation_grid_type,
            ddc::SplineBuilderClosure::HERMITE,
            ddc::SplineBuilderClosure::HERMITE,
            Solver>;

    // Helper types
    using IdxInterpGrid = Idx<interpolation_grid_type>;
    using IdxStepInterpGrid = IdxStep<interpolation_grid_type>;

private:
    std::optional<DerivReconstructionCoeffs<
            InterpolationGrid1D,
            InterpolationGrid1D,
            LocalBSplinesLower,
            LocalBSplines,
            Extremity::BACK,
            Extremity::FRONT>>
            m_reconstruction_coeffs_lower;
    std::optional<DerivReconstructionCoeffs<
            InterpolationGrid1D,
            InterpolationGrid1D,
            LocalBSplines,
            LocalBSplinesUpper,
            Extremity::BACK,
            Extremity::FRONT>>
            m_reconstruction_coeffs_upper;

    // MPI tags identifying the two logical exchange channels (a rank's lower-boundary
    // local sum, and a rank's upper-boundary local sum).
    static constexpr int s_tag_lower = 0;
    static constexpr int s_tag_upper = 1;

    // MPI information
    MPI_Comm m_comm;

    // Ranks of the adjacent cells
    int m_rank_below;
    int m_rank_above;
    int m_rank;

    // spline builder information
    interpolation_idx_range_type m_local_idx_range;
    internal_builder_type m_spline_builder;

    // The spacing between two neighbouring interpolation points (uniform mesh).
    double m_dx;

public:
    /**
     * @brief Instantiate MPISplineBuilderLocalVidal.
     *
     * @param local_idx_range The 1D index range of the interpolation mesh owned by this rank.
     * Must contain at least 11 points so that the local linear combination never reaches
     * outside the local subdomain.
     * @param comm The MPI communicator over which the global domain is split. Rank i is
     * assumed to own the subdomain immediately to the left of rank i+1.
     */
    MPISplineBuilderLocalVidal(
            interpolation_idx_range_type const& local_idx_range,
            interpolation_idx_range_type const& global_idx_range,
            MPI_Comm comm,
            int n_neighbours)
        : m_comm(comm)
        , m_local_idx_range(local_idx_range)
        , m_spline_builder(local_idx_range)
        , m_dx(ddc::discrete_space<interpolation_grid_type>().step())
    {
        int rank, comm_size;
        MPI_Comm_rank(m_comm, &rank);
        MPI_Comm_size(m_comm, &comm_size);

        m_rank = rank;

        assert(local_idx_range.size() >= n_neighbours + 1);

        bool constexpr is_periodic = (BCLower == ddc::SplineBuilderClosure::PERIODIC);

        // Save the ranks of the adjacent cells
        if (rank > 0) {
            m_rank_below = rank - 1;
            m_reconstruction_coeffs_lower.emplace(
                    interpolation_idx_range_type {
                            local_idx_range.front() - n_neighbours + 1,
                            IdxStepInterpGrid {n_neighbours}},
                    local_idx_range.take_first(IdxStepInterpGrid {n_neighbours}));
        } else if (is_periodic) {
            // Handle boundary conditions
            m_rank_below = comm_size - 1;
            interpolation_idx_range_type local_idx_range_lower(
                    global_idx_range.back() - n_neighbours + 1,
                    IdxStepInterpGrid {n_neighbours});
            m_reconstruction_coeffs_lower.emplace(
                    local_idx_range_lower,
                    local_idx_range.take_first(IdxStepInterpGrid {n_neighbours}));
        } else {
            m_rank_below = -1;
        }

        if (rank < comm_size - 1) {
            m_rank_above = rank + 1;
            m_reconstruction_coeffs_upper.emplace(
                    local_idx_range.take_last(IdxStepInterpGrid {n_neighbours}),
                    interpolation_idx_range_type {
                            local_idx_range.back(),
                            IdxStepInterpGrid {n_neighbours}});
        } else if (is_periodic) {
            // Handle boundary conditions
            m_rank_above = 0;
            m_reconstruction_coeffs_upper.emplace(
                    local_idx_range.take_last(IdxStepInterpGrid {n_neighbours}),
                    global_idx_range.take_first(IdxStepInterpGrid {n_neighbours}));
        } else {
            m_rank_above = -1;
        }
    }

    /**
     * @brief Compute the spline coefficients for a function.
     *
     * Approximates the derivative at the lower and upper boundaries of this rank's local
     * domain (via the local linear combination described in the class
     * documentation, exchanging one scalar per batch index with each real neighbour rank),
     * then calls the internal ddc::SplineBuilder with these as Hermite boundary conditions.
     *
     * @param[out] coeffs The coefficients of the spline computed by this builder.
     * @param[in] vals The values of the function on this rank's local interpolation mesh.
     * @param[in] derivs_xmin The true physical derivative at the lower bound of the *global*
     * domain. Required (and used) only on the rank owning that bound, and only when
     * BCLower == HERMITE.
     * @param[in] derivs_xmax The true physical derivative at the upper bound of the *global*
     * domain. Required (and used) only on the rank owning that bound, and only when
     * BCUpper == HERMITE.
     */
    template <
            class IdxRangeBatchedInterpolation,
            class LayoutCoeffs,
            class LayoutVals,
            class LayoutDerivs = Kokkos::layout_right>
    void operator()(
            DField<batched_basis_idx_range_type<IdxRangeBatchedInterpolation>,
                   memory_space,
                   LayoutCoeffs> coeffs,
            DConstField<IdxRangeBatchedInterpolation, memory_space, LayoutVals> local_vals,
            std::optional<DConstField<
                    batched_derivs_idx_range_type<IdxRangeBatchedInterpolation>,
                    memory_space,
                    LayoutDerivs>> global_derivs_xmin
            = std::nullopt,
            std::optional<DConstField<
                    batched_derivs_idx_range_type<IdxRangeBatchedInterpolation>,
                    memory_space,
                    LayoutDerivs>> global_derivs_xmax
            = std::nullopt) const
    {
        using IdxRangeBatch
                = ddc::detail::convert_type_seq_to_discrete_domain_t<ddc::type_seq_remove_t<
                        ddc::to_type_seq_t<IdxRangeBatchedInterpolation>,
                        ddc::to_type_seq_t<interpolation_idx_range_type>>>;
        using IdxBatch = typename IdxRangeBatch::discrete_element_type;
        Idx<deriv_type> const first_deriv(1);
        IdxRangeBatch idx_range_batch(get_idx_range(local_vals));

        using DerivsIdxRange = batched_derivs_idx_range_type<IdxRangeBatchedInterpolation>;

        IdxRange<deriv_type> const deriv_idx_range(Idx<deriv_type>(1), IdxStep<deriv_type>(1));
        DerivsIdxRange const local_derivs_domain(deriv_idx_range, get_idx_range(local_vals));

        // Allocate memory for derivatives at each edge of the function on the rank.
        // This memory is equivalent to the size of the batch index
        DFieldMem<DerivsIdxRange, memory_space> local_derivs_xmin_alloc(local_derivs_domain);
        DFieldMem<DerivsIdxRange, memory_space> local_derivs_xmax_alloc(local_derivs_domain);
        DFieldMem<DerivsIdxRange, memory_space> neighbour_derivs_xmin_alloc(local_derivs_domain);
        DFieldMem<DerivsIdxRange, memory_space> neighbour_derivs_xmax_alloc(local_derivs_domain);

        DField<DerivsIdxRange, memory_space> local_derivs_xmin(local_derivs_xmin_alloc);
        DField<DerivsIdxRange, memory_space> local_derivs_xmax(local_derivs_xmax_alloc);
        DField<DerivsIdxRange, memory_space> neighbour_derivs_xmin(neighbour_derivs_xmin_alloc);
        DField<DerivsIdxRange, memory_space> neighbour_derivs_xmax(neighbour_derivs_xmax_alloc);

        int const count = static_cast<int>(local_derivs_domain.size());
        MPI_Request lower_request;
        MPI_Request upper_request;

        // Compute the part of the linear combination which uses local function values and send
        // the result of the calculation to the neighbouring ranks.
        compute_local_deriv_component(
                local_derivs_xmin,
                lower_request,
                s_tag_lower,
                local_vals,
                m_local_idx_range.front(),
                m_rank_below,
                BCLower,
                global_derivs_xmin);
        compute_local_deriv_component(
                local_derivs_xmax,
                upper_request,
                s_tag_upper,
                local_vals,
                m_local_idx_range.back(),
                m_rank_above,
                BCUpper,
                global_derivs_xmax);

        // Collect calculation from neighbouring ranks and complete the linear combination to obtain the
        // approximations of the derivatives.
        if (m_rank_below >= 0) {
            MPI_Recv(
                    neighbour_derivs_xmin.data_handle(),
                    count,
                    MPI_DOUBLE,
                    m_rank_below,
                    s_tag_upper,
                    m_comm,
                    MPI_STATUS_IGNORE);
            MPI_Wait(&lower_request, MPI_STATUS_IGNORE);
            ddc::host_for_each(
                    //exec_space(),
                    idx_range_batch,
                    KOKKOS_LAMBDA(IdxBatch const idx_b) {
                        double const total = local_derivs_xmin(first_deriv, idx_b)
                                             + neighbour_derivs_xmin(first_deriv, idx_b);
                        local_derivs_xmin(first_deriv, idx_b) = total;
                    });
            Kokkos::fence();
        }
        if (m_rank_above >= 0) {
            MPI_Recv(
                    neighbour_derivs_xmax.data_handle(),
                    count,
                    MPI_DOUBLE,
                    m_rank_above,
                    s_tag_lower,
                    m_comm,
                    MPI_STATUS_IGNORE);
            MPI_Wait(&upper_request, MPI_STATUS_IGNORE);
            ddc::host_for_each(
                    //exec_space(),
                    idx_range_batch,
                    KOKKOS_LAMBDA(IdxBatch const idx_b) {
                        double const total = local_derivs_xmax(first_deriv, idx_b)
                                             + neighbour_derivs_xmax(first_deriv, idx_b);
                        local_derivs_xmax(first_deriv, idx_b) = total;
                    });
            Kokkos::fence();
        }

        // Use the approximations to calculate the spline interpolation
        m_spline_builder(
                coeffs,
                local_vals,
                std::optional(get_const_field(local_derivs_xmin)),
                std::optional(get_const_field(local_derivs_xmax)));
    }

    /**
     * @brief Get the whole index range on which derivatives on the lower boundary are defined.
     * @param batched_interpolation_domain The whole index range on which the interpolation points are defined.
     * @return The index range for the Derivs values.
     */
    template <class IdxRangeBatchedInterpolation>
    batched_derivs_idx_range_type<IdxRangeBatchedInterpolation> batched_derivs_xmin_domain(
            IdxRangeBatchedInterpolation const& batched_interpolation_domain) const noexcept
    {
        IdxRange<deriv_type> const deriv_idx_range(Idx<deriv_type>(1), IdxStep<deriv_type>(1));
        return batched_derivs_idx_range_type<
                IdxRangeBatchedInterpolation>(deriv_idx_range, batched_interpolation_domain);
    }

    /**
     * @brief Get the whole index range on which derivatives on the upper boundary are defined.
     * @param batched_interpolation_domain The whole index range on which the interpolation points are defined.
     * @return The index range for the Derivs values.
     */
    template <class IdxRangeBatchedInterpolation>
    batched_derivs_idx_range_type<IdxRangeBatchedInterpolation> batched_derivs_xmax_domain(
            IdxRangeBatchedInterpolation const& batched_interpolation_domain) const noexcept
    {
        IdxRange<deriv_type> const deriv_idx_range(Idx<deriv_type>(1), IdxStep<deriv_type>(1));
        return batched_derivs_idx_range_type<
                IdxRangeBatchedInterpolation>(deriv_idx_range, batched_interpolation_domain);
    }

    /**
     * @brief Compute the part of the linear combination which uses local function values and send
     * the result of the calculation to the neighbouring ranks.
     * This method should be private but cannot be due to cuda restrictions.
     *
     * @param[out] local_derivs The field to fill.
     * @param[out] request The MPI request resulting from the MPI_Isend call.
     * @param[in] send_tag The MPI tag used to send this boundary's local sum to @p neighbour_rank.
     * @param[in] local_vals The local function values.
     * @param[in] boundary_idx The local index of the boundary point (front() or back() of
     * the local interpolation index range).
     * @param[in] neighbour_rank The rank across this boundary, or -1 if there is none (true
     * global domain edge under a non-periodic boundary condition).
     * @param[in] bc The boundary condition selected for this end of the *global* domain
     * (only actually used when @p neighbour_rank == -1).
     * @param[in] global_derivs The physical derivative supplied by the caller for this end of
     * the *global* domain (only used when @p neighbour_rank == -1 and @p bc == HERMITE).
     */
    template <
            class IdxRangeBatchedDerivInterpolation,
            class IdxRangeBatchedInterpolation,
            class LayoutVals,
            class LayoutDerivs>
    void compute_local_deriv_component(
            DField<IdxRangeBatchedDerivInterpolation, memory_space> local_derivs,
            MPI_Request& request,
            int send_tag,
            DConstField<IdxRangeBatchedInterpolation, memory_space, LayoutVals> local_vals,
            Idx<interpolation_grid_type> const& boundary_idx,
            int neighbour_rank,
            ddc::SplineBuilderClosure bc,
            std::optional<DConstField<
                    IdxRangeBatchedDerivInterpolation,
                    memory_space,
                    LayoutDerivs>> const& global_derivs) const
    {
        using IdxDerivBatch = typename IdxRangeBatchedDerivInterpolation::discrete_element_type;
        using IdxRangeBatch = ddc::remove_dims_of_t<IdxRangeBatchedDerivInterpolation, deriv_type>;
        using IdxBatch = typename IdxRangeBatch::discrete_element_type;

        if (neighbour_rank < 0) {
            // True edge of a non-periodic global domain: no MPI communication.
            if (bc == ddc::SplineBuilderClosure::HOMOGENEOUS_HERMITE) {
                ddc::parallel_fill(local_derivs, 0.0);
            } else {
                assert(global_derivs.has_value());
                ddc::parallel_deepcopy(local_derivs, global_derivs.value());
            }
            Kokkos::fence();
            return;
        }

        if (send_tag == s_tag_lower) {
            assert(m_reconstruction_coeffs_lower.has_value());
            //const std::source_location location = std::source_location::current();
            ddc::host_for_each(
                    //location.function_name(),
                    //exec_space(),
                    get_idx_range(local_derivs),
                    KOKKOS_LAMBDA(IdxDerivBatch const idx_db) {
                        IdxBatch idx_b(idx_db);
                        // TODO: parallel vs host
                        local_derivs(idx_db)
                                = m_reconstruction_coeffs_lower->get_approx_deriv_patch2_component(
                                        local_vals[idx_b]);
                    });
            Kokkos::fence();
        } else {
            assert(m_reconstruction_coeffs_upper.has_value());
            ddc::host_for_each(
                    //location.function_name(),
                    //exec_space(),
                    get_idx_range(local_derivs),
                    KOKKOS_LAMBDA(IdxDerivBatch const idx_db) {
                        IdxBatch idx_b(idx_db);
                        // TODO: parallel vs host
                        local_derivs(idx_db)
                                = m_reconstruction_coeffs_upper->get_approx_deriv_patch1_component(
                                        local_vals[idx_b]);
                    });
            Kokkos::fence();
        }

        int const count = static_cast<int>(local_derivs.size());

        Kokkos::fence("Ensure derivative is calculated before sending data");
        MPI_Isend(
                local_derivs.data_handle(),
                count,
                MPI_DOUBLE,
                neighbour_rank,
                send_tag,
                m_comm,
                &request);
    }
};
