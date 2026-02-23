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
 * @tparam InterpolationGrid The discrete dimension on which interpolation points are defined.
 * @tparam Basis The basis on which the interpolation is constructed.
 */
template <
        class ExecSpace,
        class MemorySpace,
        class DataType,
        class InterpolationGrid,
        class Basis,
        class BatchedInterpolationIdxRange = IdxRange<InterpolationGrid>>
class IdentityInterpolationBuilder
    : public IInterpolationBuilder<
              ExecSpace,
              MemorySpace,
              DataType,
              InterpolationGrid,
              Basis,
              BatchedInterpolationIdxRange>
{
    using base_type = IInterpolationBuilder<
            ExecSpace,
            MemorySpace,
            DataType,
            InterpolationGrid,
            Basis,
            BatchedInterpolationIdxRange>;

public:
    /// @brief The type of the Kokkos execution space used by this class.
    using typename base_type::exec_space;

    /// @brief The type of the Kokkos memory space used by this class.
    using typename base_type::memory_space;

    /// @brief The type of the interpolation continuous dimension (continuous dimension of interest) used by this class.
    using typename base_type::continuous_dimension_type;

    /// @brief The type of the interpolation discrete dimension (discrete dimension of interest) used by this class.
    using typename base_type::interpolation_grid_type;
    /// @brief The type of the domain for the 1D interpolation mesh used by this class.
    using typename base_type::interpolation_domain_type;

    /// @brief The grid on which the interpolation coefficients should be provided.
    using typename base_type::basis_domain_type;

public:
    /// @brief The number of equations defining the boundary condition at the lower bound.
    static constexpr int s_nbc_xmin = 0;

    /// @brief The number of equations defining the boundary condition at the upper bound.
    static constexpr int s_nbc_xmax = 0;

public:
    IdentityInterpolationBuilder() = default;

    /**
     * @brief Compute the interpolation coefficients for a function.
     *
     * No calculations are necessary and the data simply needs to be copied
     * from vals to coeffs.
     *
     * @param[out] coeffs The coefficients of the spline computed by this SplineBuilder.
     * @param[in] vals The values of the function on the interpolation mesh.
     * @param[in] derivs_xmin The values of the derivatives at the lower boundary
     *                  (unused in this class).
     * @param[in] derivs_xmax The values of the derivatives at the upper boundary
     *                  (unused in this class).
     */
    void operator()(
            Field<DataType,
                  typename base_type::batched_basis_domain_type<BatchedInterpolationIdxRange>,
                  memory_space> coeffs,
            ConstField<DataType, BatchedInterpolationIdxRange, memory_space> vals,
            std::optional<ConstField<
                    DataType,
                    typename base_type::batched_derivs_domain_type<BatchedInterpolationIdxRange>,
                    memory_space>> derivs_xmin
            = std::nullopt,
            std::optional<ConstField<
                    DataType,
                    typename base_type::batched_derivs_domain_type<BatchedInterpolationIdxRange>,
                    memory_space>> derivs_xmax
            = std::nullopt) const final
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
            typename BatchedInterpolationIdxRange::discrete_vector_type nrepeat(extended_domain.size());
            BatchedInterpolationIdxRange repeat_domain(get_idx_range(vals).take_first(nrepeat));
            Kokkos::deep_copy(
                    coeffs[extended_domain].allocation_kokkos_view(),
                    vals[repeat_domain].allocation_kokkos_view());
        }
    }
};
