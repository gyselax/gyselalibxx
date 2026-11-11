// SPDX-License-Identifier: MIT
#pragma once
#include <mpi.h>

#include <cassert>
#include <limits>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <vector>

#include <ddc/ddc.hpp>
#include <ddc/kernels/splines.hpp>

#include <Kokkos_Profiling_ScopedRegion.hpp>

#include "ddc_alias_inline_functions.hpp"
#include "ddc_aliases.hpp"
#include "ddc_helper.hpp"
#include "i_interpolation_builder.hpp"
#include "mpilayout.hpp"
#include "mpitools.hpp"
#include "mpitransposealltoall.hpp"
#include "transpose.hpp"

/**
 * @brief A 1D spline builder for data distributed across MPI processes along the
 * interpolation dimension.
 *
 * The input values are distributed (evenly) across the MPI processes of the communicator
 * along the interpolation dimension. In order to build the spline the data is transposed
 * with an MPI Alltoall such that it is instead distributed across the batch dimensions.
 * The wrapped (serial) spline builder is then used to compute the spline coefficients,
 * and these coefficients are transposed back with an MPI Alltoallv such that they are
 * distributed along the B-spline dimension.
 *
 * Only periodic B-splines are currently handled. The number of cells must be divisible by
 * the number of MPI processes. If nloc = ncells / comm_size then the rank r holds:
 * - the values at the interpolation points [r*nloc, (r+1)*nloc) (see local_interpolation_domain());
 * - the spline coefficients [r*nloc, r*nloc + nloc + degree). These are all the coefficients
 *   of the B-splines which are non-zero on the cells handled by the rank. Coefficients are
 *   therefore duplicated across neighbouring ranks (degree coefficients overlap).
 *
 * As with ddc::SplineBuilder, the batch dimensions are only known when operator() is called.
 * There must be at least one batch dimension, and the batch dimensions must be distributable
 * evenly over the MPI processes (see MPILayout).
 *
 * The overlapping blocks of coefficients are sent directly from the builder's output using
 * overlapping regions of the MPI_Alltoallv send buffer, so no duplication copy is needed.
 * The coefficients are only reordered before sending if the B-spline dimension is not the
 * outermost dimension.
 *
 * This class satisfies the concepts::InterpolationBuilder1D concept.
 *
 * @tparam Builder The serial 1D spline builder (e.g. ddc::SplineBuilder) used to compute the
 *                  coefficients once the data is local along the interpolation dimension.
 *                  Its coefficient domain must contain the periodic duplicates
 *                  (ncells + degree coefficients).
 */
template <concepts::InterpolationBuilder1D Builder>
class MPIBuilderAlltoall
{
public:
    /// @brief The type of the Kokkos execution space used by this class.
    using exec_space = typename Builder::exec_space;

    /// @brief The type of the Kokkos memory space used by this class.
    using memory_space = typename Builder::memory_space;

    /// @brief The data type that the data is saved on.
    using data_type = typename InterpolationBuilderTraits<Builder>::data_type;

    /// @brief The type of the B-splines on which the coefficients are defined.
    using bsplines_type = typename Builder::bsplines_type;

    /// @brief The type of the interpolation discrete dimension.
    using interpolation_grid_type = interpolation_grid_t<Builder>;

    /// @brief The type of the interpolation discrete dimension (ddc::SplineBuilder name).
    using interpolation_discrete_dimension_type = interpolation_grid_type;

    /// @brief The type of the index range for the 1D interpolation mesh.
    using interpolation_idx_range_type = IdxRange<interpolation_grid_type>;

    /// @brief The type of the index range for the 1D interpolation mesh (ddc::SplineBuilder name).
    using interpolation_domain_type = interpolation_idx_range_type;

    /// @brief The index range for the interpolation coefficients.
    using coeff_idx_range_type = IdxRange<bsplines_type>;

    /// @brief The type of the continuous dimension of interest.
    using continuous_dimension_type = typename interpolation_grid_type::continuous_dimension_type;

    /// @brief The type of the Deriv dimension at the boundaries.
    using deriv_type = ddc::Deriv<continuous_dimension_type>;

    /// @brief The batched index range type with the interpolation grid replaced by the B-splines.
    template <class IdxRangeBatchedInterp>
    using batched_basis_idx_range_type
            = ddc::replace_dim_of_t<IdxRangeBatchedInterp, interpolation_grid_type, bsplines_type>;

    /// @brief The batched index range type with the interpolation grid replaced by the B-splines (ddc::SplineBuilder name).
    template <class IdxRangeBatchedInterp>
    using batched_spline_domain_type = batched_basis_idx_range_type<IdxRangeBatchedInterp>;

    /// @brief The batched index range type with the interpolation grid replaced by the deriv dimension.
    template <class IdxRangeBatchedInterp>
    using batched_derivs_idx_range_type
            = ddc::replace_dim_of_t<IdxRangeBatchedInterp, interpolation_grid_type, deriv_type>;

    /// @brief The batched index range type with the interpolation grid replaced by the deriv dimension (ddc::SplineBuilder name).
    template <class IdxRangeBatchedInterp>
    using batched_derivs_domain_type = batched_derivs_idx_range_type<IdxRangeBatchedInterp>;

private:
    static_assert(
            bsplines_type::is_periodic(),
            "MPIBuilderAlltoall currently only handles periodic B-splines");

private:
    Builder const& m_builder;
    MPI_Comm m_comm;
    int m_rank;
    int m_comm_size;
    // Number of cells (and interpolation points) handled by each rank.
    int m_nloc;
    // Number of coefficients held by each rank (m_nloc + degree).
    int m_ncoeffs_loc;
    // The B-splines on which the coefficients provided by the builder are defined.
    IdxRange<bsplines_type> m_builder_spline_range;

public:
    /**
     * @brief Construct the distributed spline builder.
     *
     * @param builder The serial spline builder. It must be defined on the full (global)
     *              interpolation domain. A reference to this builder is stored.
     * @param comm The MPI communicator across which the data is distributed.
     */
    MPIBuilderAlltoall(Builder const& builder, MPI_Comm comm)
        : m_builder(builder)
        , m_comm(comm)
        , m_rank(get_comm_rank(comm))
        , m_comm_size(get_comm_size(comm))
        , m_nloc(get_n_local_cells(m_comm_size))
        , m_ncoeffs_loc(m_nloc + bsplines_type::degree())
        , m_builder_spline_range(ddc::select<bsplines_type>(
                  m_builder.batched_spline_domain(m_builder.interpolation_domain())))
    {
        if (m_builder_spline_range.size() < m_comm_size * m_nloc + bsplines_type::degree()) {
            throw std::runtime_error(
                    "The builder's coefficient domain must contain the periodic duplicates "
                    "(ncells + degree coefficients)");
        }
    }

    /**
     * @brief Get the section of the interpolation domain held by the current MPI rank.
     *
     * The values passed to operator() must be defined on this index range (combined with
     * any batch dimensions).
     *
     * @return The local interpolation index range.
     */
    interpolation_idx_range_type local_interpolation_domain() const noexcept
    {
        interpolation_idx_range_type const global = m_builder.interpolation_domain();
        return interpolation_idx_range_type(
                global.front() + m_rank * m_nloc,
                IdxStep<interpolation_grid_type>(m_nloc));
    }

    /**
     * @brief Get the MPI communicator.
     *
     * @return The MPI communicator across which the data is distributed.
     */
    MPI_Comm communicator() const
    {
        return m_comm;
    }

    /**
     * @brief Get the whole (global) interpolation domain.
     *
     * @return The interpolation domain of the wrapped builder.
     */
    interpolation_idx_range_type interpolation_domain() const noexcept
    {
        return m_builder.interpolation_domain();
    }

    /**
     * @brief Get the index range on which the coefficients held by the current rank are defined.
     *
     * @param batched_interpolation_domain The local index range on which the values are defined.
     *
     * @return The local index range for the spline coefficients.
     */
    template <class IdxRangeBatchedInterp>
    batched_basis_idx_range_type<IdxRangeBatchedInterp> batched_spline_domain(
            IdxRangeBatchedInterp const& batched_interpolation_domain) const noexcept
    {
        IdxRange<bsplines_type> local_spline_range(
                m_builder_spline_range.front() + m_rank * m_nloc,
                IdxStep<bsplines_type>(m_ncoeffs_loc));
        return batched_basis_idx_range_type<
                IdxRangeBatchedInterp>(local_spline_range, batched_interpolation_domain);
    }

    /**
     * @brief Get the whole domain on which derivatives on lower boundary are defined.
     *
     * Derivatives are not used (only periodic boundary conditions are handled) so this
     * domain is empty along the deriv dimension.
     *
     * @param batched_interpolation_domain The whole domain on which the interpolation points are defined.
     *
     * @return The domain for the Derivs values.
     */
    template <class IdxRangeBatchedInterp>
    batched_derivs_idx_range_type<IdxRangeBatchedInterp> batched_derivs_xmin_domain(
            IdxRangeBatchedInterp const& batched_interpolation_domain) const noexcept
    {
        IdxRange<deriv_type> empty_deriv_range(Idx<deriv_type>(0), IdxStep<deriv_type>(0));
        return batched_derivs_idx_range_type<
                IdxRangeBatchedInterp>(empty_deriv_range, batched_interpolation_domain);
    }

    /**
     * @brief Get the whole domain on which derivatives on upper boundary are defined.
     *
     * Derivatives are not used (only periodic boundary conditions are handled) so this
     * domain is empty along the deriv dimension.
     *
     * @param batched_interpolation_domain The whole domain on which the interpolation points are defined.
     *
     * @return The domain for the Derivs values.
     */
    template <class IdxRangeBatchedInterp>
    batched_derivs_idx_range_type<IdxRangeBatchedInterp> batched_derivs_xmax_domain(
            IdxRangeBatchedInterp const& batched_interpolation_domain) const noexcept
    {
        IdxRange<deriv_type> empty_deriv_range(Idx<deriv_type>(0), IdxStep<deriv_type>(0));
        return batched_derivs_idx_range_type<
                IdxRangeBatchedInterp>(empty_deriv_range, batched_interpolation_domain);
    }

    /**
     * @brief Compute the spline coefficients of a function distributed across MPI ranks.
     *
     * This is a collective operation over the communicator. All ranks must call it with the
     * same batch index ranges.
     *
     * @param[out] coeffs The local spline coefficients. They must be defined on
     *              batched_spline_domain(get_idx_range(vals)).
     * @param[in] vals The local values of the function at the interpolation points. They must
     *              be defined on local_interpolation_domain() and the (full) batch index ranges.
     * @param[in] derivs_xmin Unused (periodic boundary conditions). Must be std::nullopt.
     * @param[in] derivs_xmax Unused (periodic boundary conditions). Must be std::nullopt.
     */
    template <class IdxRangeBatchedInterp>
    void operator()(
            Field<data_type, batched_basis_idx_range_type<IdxRangeBatchedInterp>, memory_space>
                    coeffs,
            ConstField<data_type, IdxRangeBatchedInterp, memory_space> vals,
            std::optional<ConstField<
                    data_type,
                    batched_derivs_idx_range_type<IdxRangeBatchedInterp>,
                    memory_space>> derivs_xmin
            = std::nullopt,
            std::optional<ConstField<
                    data_type,
                    batched_derivs_idx_range_type<IdxRangeBatchedInterp>,
                    memory_space>> derivs_xmax
            = std::nullopt) const
    {
        static_assert(
                ddc::in_tags_v<interpolation_grid_type, ddc::to_type_seq_t<IdxRangeBatchedInterp>>);

        using IdxRangeBatch = ddc::remove_dims_of_t<IdxRangeBatchedInterp, interpolation_grid_type>;
        static_assert(
                IdxRangeBatch::rank() > 0,
                "MPIBuilderAlltoall requires at least one batch dimension");

        using IdxRangeBatchedSpline = batched_basis_idx_range_type<IdxRangeBatchedInterp>;

        using LayoutValsSplit = MPILayout<IdxRangeBatchedInterp, interpolation_grid_type>;
        using LayoutBatchSplit = mpi_layout_t<IdxRangeBatchedInterp, IdxRangeBatch>;
        using ValsTranspose = MPITransposeAllToAll<LayoutValsSplit, LayoutBatchSplit>;

        // MPI tags (e.g. MPIDim<GridY>) describing the rank where a batch block is found.
        using batch_mpi_dims = ddcHelper::
                apply_template_to_type_seq_t<MPIDim, ddc::to_type_seq_t<IdxRangeBatch>>;
        // Coefficients with the B-splines as the outermost dimension so that the block sent to
        // each rank is contiguous.
        using IdxRangeSend = ddc::combine_t<coeff_idx_range_type, IdxRangeBatch>;
        // Received coefficients: one block per source rank.
        using IdxRangeRecv = ddc::combine_t<
                ddc::detail::convert_type_seq_to_discrete_domain_t<batch_mpi_dims>,
                IdxRangeSend>;
        // Output coefficients with the batch dimensions split by the MPI tags.
        using IdxRangeOutMPI
                = ddc::detail::convert_type_seq_to_discrete_domain_t<insert_mpi_tags_into_seq_t<
                        batch_mpi_dims,
                        ddc::to_type_seq_t<IdxRangeBatchedSpline>>>;

        assert(!derivs_xmin.has_value() && !derivs_xmax.has_value());
        assert(ddc::select<interpolation_grid_type>(get_idx_range(vals))
               == local_interpolation_domain());
        assert(get_idx_range(coeffs) == batched_spline_domain(get_idx_range(vals)));

        Kokkos::Profiling::ScopedRegion region("(GSLX) MPIBuilderAlltoall");

        // The batch dimensions are not distributed in the input layout so the global index
        // range can be reconstructed locally.
        IdxRangeBatch const batch_idx_range(get_idx_range(vals));
        IdxRangeBatchedInterp const
                global_idx_range(m_builder.interpolation_domain(), batch_idx_range);

        // 1. Transpose the values so they are distributed along the batch dimensions.
        ValsTranspose const vals_transpose(global_idx_range, m_comm);
        assert(vals_transpose.template get_local_idx_range<LayoutValsSplit>()
               == get_idx_range(vals));
        IdxRangeBatchedInterp const vals_batch_distrib_idx_range
                = vals_transpose.template get_local_idx_range<LayoutBatchSplit>();
        FieldMem<data_type, IdxRangeBatchedInterp, memory_space> vals_split_batch_alloc(
                vals_batch_distrib_idx_range);
        vals_transpose.template transpose_to<
                LayoutBatchSplit>(exec_space(), get_field(vals_split_batch_alloc), vals);

        // 2. Build the spline on the data which is now local along the interpolation dimension.
        FieldMem<data_type, IdxRangeBatchedSpline, memory_space> coeffs_split_batch_alloc(
                m_builder.batched_spline_domain(vals_batch_distrib_idx_range));
        m_builder(get_field(coeffs_split_batch_alloc), get_const_field(vals_split_batch_alloc));

        // 3. Ensure the B-splines are the outermost dimension so the (overlapping) block sent
        //    to each rank is contiguous. This is a view (no copy) if this is already the case.
        // This is not done earlier as the position of the batch dims of the coeffs must match
        // the position of the batch dims of the values.
        auto send_buffer = ddcHelper::create_transpose_mirror_view_and_copy<
                IdxRangeSend>(exec_space(), get_field(coeffs_split_batch_alloc));

        // 4. Describe how the batch dimensions are distributed (as in MPITransposeAllToAll the
        //    batch blocks are labelled using the local batch index range of rank 0).
        LayoutBatchSplit batch_layout;
        IdxRangeBatch const batch_ref_idx_range(
                batch_layout.distribute_idx_range(global_idx_range, m_comm_size, 0));
        auto const batch_mpi_idx_range = get_mpi_idx_range(batch_ref_idx_range, batch_idx_range);
        assert(batch_mpi_idx_range.size() == m_comm_size);

        // 5. Create a view on the output with the batch dimensions split by the MPI tags and
        //    a buffer to receive the blocks from each rank (MPI tags outermost).
        IdxRangeOutMPI const coeffs_mpi_idx_range(
                batch_mpi_idx_range,
                ddc::select<bsplines_type>(get_idx_range(coeffs)),
                batch_ref_idx_range);
        Field<data_type, IdxRangeOutMPI, memory_space>
                coeffs_mpi(coeffs.data_handle(), coeffs_mpi_idx_range);
        auto recv_buffer
                = ddcHelper::create_transpose_mirror<IdxRangeRecv>(exec_space(), coeffs_mpi);

        // 6. Send each rank the coefficients of the B-splines supported on its cells.
        //    The blocks sent to neighbouring ranks overlap by degree coefficients.
        std::size_t const batch_loc_size = batch_ref_idx_range.size();
        std::size_t const block_size = m_ncoeffs_loc * batch_loc_size;
        if (m_builder_spline_range.size() * batch_loc_size > std::numeric_limits<int>::max()
            || m_comm_size * block_size > std::numeric_limits<int>::max()) {
            throw std::runtime_error("The coefficient buffers are too large for MPI_Alltoallv");
        }
        std::vector<int> const counts(m_comm_size, static_cast<int>(block_size));
        std::vector<int> send_displs(m_comm_size);
        std::vector<int> recv_displs(m_comm_size);
        for (int p = 0; p < m_comm_size; ++p) {
            send_displs[p] = static_cast<int>(p * m_nloc * batch_loc_size);
            recv_displs[p] = static_cast<int>(p * block_size);
        }

        exec_space().fence("MPIBuilderAlltoall: fence before MPI_Alltoallv");
        MPI_Alltoallv(
                get_const_field(send_buffer).data_handle(),
                counts.data(),
                send_displs.data(),
                MPI_type_descriptor_t<data_type>,
                get_field(recv_buffer).data_handle(),
                counts.data(),
                recv_displs.data(),
                MPI_type_descriptor_t<data_type>,
                m_comm);

        // 7. Reorder the received blocks into the output layout.
        if constexpr (!ddc::is_borrowed_chunk_v<decltype(recv_buffer)>) {
            transpose_layout(exec_space(), coeffs_mpi, get_const_field(recv_buffer));
        }
    }

private:
    static int get_comm_rank(MPI_Comm comm)
    {
        int rank;
        MPI_Comm_rank(comm, &rank);
        return rank;
    }

    static int get_comm_size(MPI_Comm comm)
    {
        int size;
        MPI_Comm_size(comm, &size);
        return size;
    }

    static int get_n_local_cells(int comm_size)
    {
        int const ncells = ddc::discrete_space<bsplines_type>().ncells();
        if (ncells % comm_size != 0) {
            std::ostringstream error_msg;
            error_msg << "The number of cells (" << ncells
                      << ") must be divisible by the number of MPI ranks (" << comm_size << ")";
            throw std::runtime_error(error_msg.str());
        }
        return ncells / comm_size;
    }

    template <class... BatchDims>
    static IdxRange<MPIDim<BatchDims>...> get_mpi_idx_range(
            IdxRange<BatchDims...> local_idx_range,
            IdxRange<BatchDims...> global_idx_range)
    {
        Idx<MPIDim<BatchDims>...> start(Idx<MPIDim<BatchDims>>(0)...);
        IdxStep<MPIDim<BatchDims>...> size(IdxStep<MPIDim<BatchDims>>(
                ddc::select<BatchDims>(global_idx_range).size()
                / ddc::select<BatchDims>(local_idx_range).size())...);
        return IdxRange<MPIDim<BatchDims>...>(start, size);
    }
};
