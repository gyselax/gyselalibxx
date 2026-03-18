// SPDX-License-Identifier: MIT

#pragma once

#include "ddc_aliases.hpp"
#include "i_interpolation_builder.hpp"

/**
 * @brief A builder class for copying data.
 *
 * A class which contains an operator () which can be used to build an interpolation
 * of a function. This class handles the case where no calculations are necessary and
 * the data simply needs to be copied.
 *
 * @tparam ExecSpace The Kokkos execution space on which the spline approximation is performed.
 * @tparam MemorySpace The Kokkos memory space on which the data (interpolation function and splines coefficients) is stored.
 * @tparam DataType The data type of the field values and coefficients.
 * @tparam InterpolationGrid The discrete dimension on which interpolation points are defined.
 * @tparam Basis The basis on which the interpolation is constructed.
 */
template <class ExecSpace, class MemorySpace, class DataType, class InterpolationGrid, class Basis>
class IdentityInterpolationBuilder
{
public:
    /// @brief The type of the Kokkos execution space used by this class.
    using exec_space = ExecSpace;

    /// @brief The type of the Kokkos memory space used by this class.
    using memory_space = MemorySpace;

    /// @brief The data type that the data is saved on.
    using data_type = DataType;

    /// @brief The type of the interpolation continuous dimension (continuous dimension of interest) used by this class.
    using continuous_dimension_type = typename InterpolationGrid::continuous_dimension_type;

    /// @brief The type of the interpolation discrete dimension (discrete dimension of interest) used by this class.
    using interpolation_grid_type = InterpolationGrid;

    /// @brief The type of the domain for the 1D interpolation mesh used by this class.
    using interpolation_idx_range_type = IdxRange<interpolation_grid_type>;

    /// @brief The grid on which the interpolation coefficients should be provided.
    using basis_domain_type = typename Basis::template Impl<Basis, MemorySpace>::knot_grid;

    /// @brief The type of the Deriv dimension at the boundaries.
    using deriv_type = ddc::Deriv<continuous_dimension_type>;

    /// @brief The batched domain type with interpolation_grid_type replaced by basis_domain_type.
    template <class BatchedInterpolationIdxRange>
    using batched_basis_idx_range_type = ddc::replace_dim_of_t<
            BatchedInterpolationIdxRange,
            interpolation_grid_type,
            basis_domain_type>;

    /// @brief The batched domain type with interpolation_grid_type replaced by deriv_type.
    template <class BatchedInterpolationIdxRange>
    using batched_derivs_idx_range_type = ddc::
            replace_dim_of_t<BatchedInterpolationIdxRange, interpolation_grid_type, deriv_type>;

public:
    /// @brief The number of equations defining the boundary condition at the lower bound.
    static constexpr int s_nbe_xmin = 0;

    /// @brief The number of equations defining the boundary condition at the upper bound.
    static constexpr int s_nbe_xmax = 0;

public:
    IdentityInterpolationBuilder() = default;

    /**
     * @brief Compute the interpolation coefficients for a function.
     *
     * No calculations are necessary and the data simply needs to be copied
     * from vals to coeffs.
     *
     * @param[out] coeffs The coefficients of the interpolation computed by this builder.
     * @param[in] vals The values of the function on the interpolation mesh.
     * @param[in] derivs_xmin The values of the derivatives at the lower boundary
     *                  (unused in this class).
     * @param[in] derivs_xmax The values of the derivatives at the upper boundary
     *                  (unused in this class).
     */
    template <class BatchedInterpolationIdxRange>
    void operator()(
            Field<DataType,
                  batched_basis_idx_range_type<BatchedInterpolationIdxRange>,
                  memory_space> coeffs,
            ConstField<DataType, BatchedInterpolationIdxRange, memory_space> vals,
            std::optional<ConstField<
                    DataType,
                    batched_derivs_idx_range_type<BatchedInterpolationIdxRange>,
                    memory_space>> derivs_xmin
            = std::nullopt,
            std::optional<ConstField<
                    DataType,
                    batched_derivs_idx_range_type<BatchedInterpolationIdxRange>,
                    memory_space>> derivs_xmax
            = std::nullopt) const
    {
        IdxRange<basis_domain_type> bp_idx_range
                = ddc::discrete_space<Basis>().break_point_domain().remove_last(
                        IdxStep<basis_domain_type>(static_cast<int>(Basis::is_periodic())));
        Kokkos::deep_copy(
                coeffs[bp_idx_range].allocation_kokkos_view(),
                vals.allocation_kokkos_view());
        if constexpr (Basis::is_periodic()) {
            IdxRange<basis_domain_type> extended_domain(
                    ddc::discrete_space<Basis>().full_domain().remove_first(
                            bp_idx_range.extents()));
            typename BatchedInterpolationIdxRange::discrete_vector_type nrepeat(
                    extended_domain.size());
            BatchedInterpolationIdxRange repeat_domain(get_idx_range(vals).take_first(nrepeat));
            Kokkos::deep_copy(
                    coeffs[extended_domain].allocation_kokkos_view(),
                    vals[repeat_domain].allocation_kokkos_view());
        }
    }

    /**
     * @brief Get the whole domain on which derivatives on lower boundary are defined.
     *
     * This is only used with BoundCond::HERMITE boundary conditions.
     *
     * @param batched_interpolation_domain The whole domain on which the interpolation points are defined.
     *
     * @return The domain for the Derivs values.
     */
    template <class BatchedInterpolationIdxRange>
    batched_derivs_idx_range_type<BatchedInterpolationIdxRange> batched_derivs_xmin_domain(
            BatchedInterpolationIdxRange const& batched_interpolation_domain) const noexcept
    {
        IdxRange<deriv_type> empty_deriv_range(Idx<deriv_type>(0), IdxStep<deriv_type>(0));
        return batched_derivs_idx_range_type<
                BatchedInterpolationIdxRange>(empty_deriv_range, batched_interpolation_domain);
    }

    /**
     * @brief Get the whole domain on which derivatives on upper boundary are defined.
     *
     * This is only used with BoundCond::HERMITE boundary conditions.
     *
     * @param batched_interpolation_domain The whole domain on which the interpolation points are defined.
     *
     * @return The domain for the Derivs values.
     */
    template <class BatchedInterpolationIdxRange>
    batched_derivs_idx_range_type<BatchedInterpolationIdxRange> batched_derivs_xmax_domain(
            BatchedInterpolationIdxRange const& batched_interpolation_domain) const noexcept
    {
        IdxRange<deriv_type> empty_deriv_range(Idx<deriv_type>(0), IdxStep<deriv_type>(0));
        return batched_derivs_idx_range_type<
                BatchedInterpolationIdxRange>(empty_deriv_range, batched_interpolation_domain);
    }

    /**
     * @brief Get the whole domain on which the interpolation coefficients are defined.
     *
     * @param batched_interpolation_domain The whole domain on which the interpolation points are defined.
     *
     * @return The domain for the interpolation coefficients.
     */
    template <class BatchedInterpolationIdxRange>
    batched_basis_idx_range_type<BatchedInterpolationIdxRange> batched_basis_idx_range(
            BatchedInterpolationIdxRange const& batched_interpolation_domain) const noexcept
    {
        return batched_basis_idx_range_type<BatchedInterpolationIdxRange>(
                ddc::discrete_space<Basis>().full_domain(),
                batched_interpolation_domain);
    }
};
