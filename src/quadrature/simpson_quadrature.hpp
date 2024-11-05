// SPDX-License-Identifier: MIT

#pragma once
#include <ddc/ddc.hpp>

#include "ddc_aliases.hpp"
#include "geometry_descriptors.hpp"
#include "quadrature_coeffs_nd.hpp"
#include "trapezoid_quadrature.hpp"

/**
 * @brief Get the Simpson coefficients in 1D.
 *
 * Calculate the quadrature coefficients for the Simpson method defined on the provided index range.
 * The non-uniform form of the Simpson quadrature is:
 *
 * (x_3-x_1)(2-(x_3-x_2)/(x_2-x_1)) / 6 f(x_1) + (x_3-x_1)^3 / 6 / (x_3-x_2) / (x_2-x_1) f(x_2) + (x_3-x_1)(2-(x_2-x_1)/(x_3-x_2)) / 6 f(x_3)
 *
 * @param[in] coefficients
 * 	The field where the quadrature coefficients should be saved.
 */
template <class ExecSpace, class Grid1D>
void fill_simpson_quadrature_coefficients_1d(
        DField<IdxRange<Grid1D>, typename ExecSpace::memory_space> coefficients)
{
    IdxRange<Grid1D> idx_range = get_idx_range(coefficients);
    if constexpr (Grid1D::continuous_dimension_type::PERIODIC) {
        if (idx_range.size() % 2 != 0) {
            throw std::runtime_error(
                    "Simpson quadrature requires 2n quadrature points for a periodic direction.");
        }
    } else {
        if (idx_range.size() % 2 == 0) {
            throw std::runtime_error("Simpson quadrature requires 2n+1 quadrature points for a "
                                     "non-periodic direction.");
        }
    }

    Kokkos::parallel_for(
            "centre_left",
            Kokkos::RangePolicy<ExecSpace>(
                    0,
                    idx_range.size() / 2 - int(Grid1D::continuous_dimension_type::PERIODIC)),
            KOKKOS_LAMBDA(const int i) {
                Idx<Grid1D> idx_c = idx_range.front() + i * 2 + 1;
                double const dx_l = distance_at_left(idx_c);
                double const dx_r = distance_at_right(idx_c);
                double const dx_sum = dx_l + dx_r;
                coefficients(idx_c - 1) = dx_sum * (2. * dx_l - dx_r) / (6. * dx_l);
                coefficients(idx_c) = dx_sum * dx_sum * dx_sum / (6. * dx_l * dx_r);
            });

    // The incrementation is done in a separate loop to avoid race conditions
    Kokkos::parallel_for(
            "overlap_increment",
            Kokkos::RangePolicy<ExecSpace>(
                    0,
                    idx_range.size() / 2 - int(Grid1D::continuous_dimension_type::PERIODIC) - 1),
            KOKKOS_LAMBDA(const int i) {
                Idx<Grid1D> idx_c = idx_range.front() + i * 2 + 1;
                double const dx_l = distance_at_left(idx_c);
                double const dx_r = distance_at_right(idx_c);
                double const dx_sum = dx_l + dx_r;
                coefficients(idx_c + 1) += dx_sum * (2. * dx_r - dx_l) / (6. * dx_r);
            });

    Kokkos::parallel_for(
            "bounds",
            Kokkos::RangePolicy<ExecSpace>(0, 1),
            KOKKOS_LAMBDA(const int i) {
                Idx<Grid1D> idx_c
                        = idx_range.back() - 1 - int(Grid1D::continuous_dimension_type::PERIODIC);
                double const dx_l = distance_at_left(idx_c);
                double const dx_r = distance_at_right(idx_c);
                double const dx_sum = dx_l + dx_r;
                coefficients(idx_c + 1) = dx_sum * (2. * dx_r - dx_l) / (6. * dx_r);
            });

    if constexpr (Grid1D::continuous_dimension_type::PERIODIC) {
        Kokkos::parallel_for(
                "bounds",
                Kokkos::RangePolicy<ExecSpace>(0, 1),
                KOKKOS_LAMBDA(const int i) {
                    Idx<Grid1D> idx_c = idx_range.back();
                    double const dx_l = distance_at_left(idx_c);
                    double const dx_r = distance_at_right(idx_c);
                    double const dx_sum = dx_l + dx_r;
                    coefficients(idx_c - 1) += dx_sum * (2. * dx_l - dx_r) / (6. * dx_l);
                    coefficients(idx_c) = dx_sum * dx_sum * dx_sum / (6. * dx_l * dx_r);
                    coefficients(idx_range.front()) += dx_sum * (2. * dx_r - dx_l) / (6. * dx_r);
                });
    }
}

/**
 * @brief Get the Simpson coefficients in 1D.
 *
 * Calculate the quadrature coefficients for the Simpson method defined on the provided index range.
 * The non-uniform form of the Simpson quadrature is:
 *
 * (x_3-x_1)(2-(x_3-x_2)/(x_2-x_1)) / 6 f(x_1) + (x_3-x_1)^3 / 6 / (x_3-x_2) / (x_2-x_1) f(x_2) + (x_3-x_1)(2-(x_2-x_1)/(x_3-x_2)) / 6 f(x_3)
 *
 * @param[in] idx_range
 * 	The index range on which the quadrature will be carried out.
 *
 * @return The quadrature coefficients for the Simpson method defined on the provided index range.
 */
template <class ExecSpace, class Grid1D>
DFieldMem<IdxRange<Grid1D>, typename ExecSpace::memory_space> simpson_quadrature_coefficients_1d(
        IdxRange<Grid1D> const& idx_range)
{
    DFieldMem<IdxRange<Grid1D>, typename ExecSpace::memory_space> coefficients_alloc(idx_range);
    fill_simpson_quadrature_coefficients_1d<ExecSpace>(get_field(coefficients_alloc));
    return coefficients_alloc;
}

/**
 * @brief Get the Simpson coefficients in 1D.
 *
 * If the number of grid points is not compatible with the Simpson quadrature scheme then
 * use a trapezoid formula over one cell at the specified extremity.
 *
 * @param[in] idx_range
 * 	The index range on which the quadrature will be carried out.
 * @param[in] trapezoid_extremity
 *  The extremity where the trapezoid quadrature may be used.
 *
 * @return The quadrature coefficients for the Simpson method defined on the provided index range.
 */
template <class ExecSpace, class Grid1D>
DFieldMem<IdxRange<Grid1D>, typename ExecSpace::memory_space>
simpson_trapezoid_quadrature_coefficients_1d(
        IdxRange<Grid1D> const& idx_range,
        Extremity trapezoid_extremity)
{
    static_assert(
            !Grid1D::continuous_dimension_type::PERIODIC,
            "The extremity is non-sensical in a Periodic dimension");
    try {
        return simpson_quadrature_coefficients_1d<ExecSpace>(idx_range);
    } catch (const std::runtime_error& error) {
        DFieldMem<IdxRange<Grid1D>, typename ExecSpace::memory_space> coefficients_alloc(idx_range);
        DField<IdxRange<Grid1D>, typename ExecSpace::memory_space> coefficients(
                get_field(coefficients_alloc));
        IdxStep<Grid1D> npts_to_remove(1);
        if (trapezoid_extremity == Extremity::FRONT) {
            fill_simpson_quadrature_coefficients_1d<ExecSpace>(
                    coefficients_alloc[idx_range.remove_first(npts_to_remove)]);
            Kokkos::parallel_for(
                    "trapezoid_region",
                    Kokkos::RangePolicy<ExecSpace>(0, 1),
                    KOKKOS_LAMBDA(const int i) {
                        Idx<Grid1D> idx_r = idx_range.back();
                        double dx_cell = distance_at_left(idx_r);
                        coefficients(idx_r) = dx_cell;
                        coefficients(idx_r - 1) += dx_cell;
                    });
        } else {
            fill_simpson_quadrature_coefficients_1d<ExecSpace>(
                    coefficients_alloc[idx_range.remove_last(npts_to_remove)]);
            Kokkos::parallel_for(
                    "trapezoid_region",
                    Kokkos::RangePolicy<ExecSpace>(0, 1),
                    KOKKOS_LAMBDA(const int i) {
                        Idx<Grid1D> idx_l = idx_range.front();
                        double dx_cell = distance_at_right(idx_l);
                        coefficients(idx_l) = dx_cell;
                        coefficients(idx_l + 1) += dx_cell;
                    });
        }
        return coefficients_alloc;
    }
}

/**
 * @brief Get the simpson coefficients in ND.
 *
 * Calculate the quadrature coefficients for the simpson method defined on the provided index range.
 *
 * @tparam ExecSpace Execution space, depends on Kokkos.
 *
 * @param[in] idx_range
 * 	The idx_range on which the quadrature will be carried out.
 *
 * @return The quadrature coefficients for the trapezoid method defined on the provided idx_range.
 *         The allocation place (host or device ) will depend on the ExecSpace.
 */
template <class ExecSpace, class... ODims>
DFieldMem<IdxRange<ODims...>, typename ExecSpace::memory_space> simpson_quadrature_coefficients(
        IdxRange<ODims...> const& idx_range)
{
    return quadrature_coeffs_nd<ExecSpace, ODims...>(
            idx_range,
            (std::function<DFieldMem<IdxRange<ODims>, typename ExecSpace::memory_space>(
                     IdxRange<ODims>)>(simpson_quadrature_coefficients_1d<ExecSpace, ODims>))...);
}
