// SPDX-License-Identifier: MIT

#pragma once

#include <ddc/ddc.hpp>

#include "ddc_aliases.hpp"
#include "edge.hpp"
#include "edge_transformation.hpp"
#include "types.hpp"

/**
 * @brief Compute the derivative of an equivalent global spline 
 * at the interface between two patches. 
 * 
 * For a given Interface, this operator compute the coefficients 
 * a, b and c of the following relation: 
 * @f$ s'(X_I) = c + a s'(X_{I+1}) + b s'(X_{I+1})@f$, 
 * 
 * with 
 * * s'(X_I) the derivative at the Interface of the two local splines on the patches. 
 * * s'(X_{I+1}) the derivative at the next (or right) Interface of the local splines 
 *      on patch 1. 
 * * s'(X_{I-1}) the derivative at the previous (or left) Interface of the local splines 
 *      on patch 2. 
 * * a,b scalars depending on the meshes of the two patches. 
 * * c a linear combination of the function values on the two patches. 
 *      The weights in the combination linear depends on the  meshes of the two patches. 
 *      It can be written as a sum: @f$ c = \sum_{k = -N2}^{N1} \omega_k f_k@f$, 
 *      with @f$ \omega_k @f$ the weights, 
 *          @f$ f_k @f$ the function values, 
 *          @f$ N1 @f$ the number of cells in patch 1, and
 *          @f$ N2 @f$ the number of cells in patch 2. 
 * 
 * Scheme of the two patches: 
 *   X_{I-1}  X_I   X_{I+1}
 *      _______ _______
 *     |       |       |
 *     |   2   |   1   |
 *     |_______|_______|
 *          + ← → +
 *          - → ← -
 * Orientation of an equivalent global mesh : →
 * 
 * All the formulae and more details are given in the README.md. 
 * 
 * @tparam Interface The interface between two patches where we want 
 * to compute the derivatives. 
 * @tparam Bound1 The interpolation condition for the spline on the patch 1
 * of the interface. By default the value is set to ddc::BoundCond::HERMITE. 
 * If you want to apply a special treatment for the additional interpolation 
 * point in case of use of interpolation points as closure conditions, please 
 * precise "ddc::BoundCond::GREVILLE" in this template parameter. 
 * @tparam Bound2 The interpolation condition for the spline on the patch 2
 * of the interface. Same comment as for Bound1.  
 * 
 * @warning The applied method only works for interpolation points located on 
 * the break points. An exception can be made for an addition interpolation point
 * on the boundary cells. 
 * Please be sure to initialise the discrete space of your Grid on the break points
 * (especially in the non-uniform case). 
 */
template <
        class Interface,
        ddc::BoundCond Bound1 = ddc::BoundCond::HERMITE,
        ddc::BoundCond Bound2 = ddc::BoundCond::HERMITE>
class SingleInterfaceDerivativesCalculator
{
    static_assert(
            (!std::is_same_v<typename Interface::Edge1, OutsideEdge>)&&(
                    !std::is_same_v<typename Interface::Edge2, OutsideEdge>),
            "The interface cannot be an interface with the outside domain.");

    using EdgePerpGrid1 = typename Interface::Edge1::perpendicular_grid;
    using EdgePerpGrid2 = typename Interface::Edge2::perpendicular_grid;

    using EdgeParGrid1 = typename Interface::Edge1::parallel_grid;
    using EdgeParGrid2 = typename Interface::Edge2::parallel_grid;

    using Patch1 = typename Interface::Edge1::associated_patch;
    using Patch2 = typename Interface::Edge2::associated_patch;

    using IdxRange2D_1 = typename Patch1::IdxRange12;
    using IdxRange2D_2 = typename Patch2::IdxRange12;

    using IdxRange1DPerp_1 = IdxRange<EdgePerpGrid1>;
    using IdxRange1DPerp_2 = IdxRange<EdgePerpGrid2>;
    using IdxRange1DPar_1 = IdxRange<EdgeParGrid1>;
    using IdxRange1DPar_2 = IdxRange<EdgeParGrid2>;

    using Idx1D_1 = Idx<EdgePerpGrid1>;
    using Idx1D_2 = Idx<EdgePerpGrid2>;

    static constexpr bool is_cell_bound_1_with_extra_interpol_pt
            = (Bound1 == ddc::BoundCond::GREVILLE);
    static constexpr bool is_cell_bound_2_with_extra_interpol_pt
            = (Bound2 == ddc::BoundCond::GREVILLE);
    // TODO: Maybe add a condition to be sure that we are taking all the cells of the
    // domain before applying special treatment.


    IdxRange1DPerp_1 const m_idx_range_perp_1;
    IdxRange1DPerp_2 const m_idx_range_perp_2;

    Extremity const m_extremity_1;
    Extremity const m_extremity_2;

    double m_coeff_deriv_patch_2;
    double m_coeff_deriv_patch_1;

    DFieldMem<IdxRange1DPerp_1, Kokkos::HostSpace> m_coeff_values_patch_1_alloc;
    DFieldMem<IdxRange1DPerp_2, Kokkos::HostSpace> m_coeff_values_patch_2_alloc;

    DField<IdxRange1DPerp_1, Kokkos::HostSpace> m_coeff_values_patch_1;
    DField<IdxRange1DPerp_2, Kokkos::HostSpace> m_coeff_values_patch_2;

public:
    /**
     * @brief Instantiate SingleInterfaceDerivativesCalculator. 
     * 
     * @anchor SingleInterfaceDerivativesCalculatorInstantiator
     * 
     * It computes the coefficients a and b, and the weights @f$\omega_k@f$
     * and stores the values in the class. 
     * If you want to compute the derivatives with an exact formula, 
     * please provide the index ranges on the full domain of the patches. 
     * If you want to use an approximation (e.g. use only 5 cells), 
     * please provide the index ranges on a reduced number of cells (e.g. 5 cells). 
     *  
     * @param idx_range_1d_1 1D index range perpendicular to the Interface, 
     * on the patch 1. 
     * @param idx_range_1d_2 1D index range perpendicular to the Interface, 
     * on the patch 2. 
     */
    SingleInterfaceDerivativesCalculator(
            IdxRange1DPerp_1 const& idx_range_1d_1,
            IdxRange1DPerp_2 const& idx_range_1d_2)
        : m_idx_range_perp_1(idx_range_1d_1)
        , m_idx_range_perp_2(idx_range_1d_2)
        , m_extremity_1(Interface::Edge1::extremity)
        , m_extremity_2(Interface::Edge2::extremity)
        , m_coeff_values_patch_1_alloc(m_idx_range_perp_1)
        , m_coeff_values_patch_2_alloc(m_idx_range_perp_2)
        , m_coeff_values_patch_1(m_coeff_values_patch_1_alloc)
        , m_coeff_values_patch_2(m_coeff_values_patch_2_alloc)
    {
        if constexpr (
                ddc::is_uniform_point_sampling_v<
                        EdgePerpGrid1> && ddc::is_uniform_point_sampling_v<EdgePerpGrid2>) {
            set_coefficients_uniform_per_patch_case();
        } else {
            set_coefficients_non_uniform_case();
        }
    }

    /**
     * @brief Instantiate SingleInterfaceDerivativesCalculator. 
     * See @ref SingleInterfaceDerivativesCalculatorInstantiator.
     * @param idx_range_1d_1 1D index range perpendicular to the Interface, 
     * on the patch 1. 
     * @param idx_range_1d_2 1D index range perpendicular to the Interface, 
     * on the patch 2.
     */
    SingleInterfaceDerivativesCalculator(
            IdxRange1DPerp_2 const& idx_range_1d_2,
            IdxRange1DPerp_1 const& idx_range_1d_1)
        : SingleInterfaceDerivativesCalculator(idx_range_1d_1, idx_range_1d_2)
    {
    }

    /**
     * @brief Instantiate SingleInterfaceDerivativesCalculator. 
     * See @ref SingleInterfaceDerivativesCalculatorInstantiator.
     * @param idx_range_1 2D index range on the patch 1. 
     * @param idx_range_2 2D index range on the patch 2. 
     */
    SingleInterfaceDerivativesCalculator(
            IdxRange2D_1 const& idx_range_1,
            IdxRange2D_2 const& idx_range_2)
        : SingleInterfaceDerivativesCalculator(
                IdxRange1DPerp_1(idx_range_1),
                IdxRange1DPerp_2(idx_range_2))
    {
    }

    /**
     * @brief Instantiate SingleInterfaceDerivativesCalculator. 
     * See @ref SingleInterfaceDerivativesCalculatorInstantiator.
     * @param idx_range_1 2D index range on the patch 1. 
     * @param idx_range_2 2D index range on the patch 2. 
     */
    SingleInterfaceDerivativesCalculator(
            IdxRange2D_2 const& idx_range_2,
            IdxRange2D_1 const& idx_range_1)
        : SingleInterfaceDerivativesCalculator(
                IdxRange1DPerp_1(idx_range_1),
                IdxRange1DPerp_2(idx_range_2))
    {
    }


    /**
     * @brief Get the coefficient (a) in front of the derivative on the patch 1 of the given Interface.
     * @return The value of the coefficient a. 
     */
    double get_coeff_deriv_patch_1() const
    {
        return m_coeff_deriv_patch_1;
    }

    /**
     * @brief Get the coefficient (b) in front of the derivative on the patch 2 of the given Interface.
     * @return The value of the coefficient b. 
     */
    double get_coeff_deriv_patch_2() const
    {
        return m_coeff_deriv_patch_2;
    }

    /**
     * @brief Get the linear combination of the function values (c).
     * 
     * @anchor get_function_coefficients
     * 
     * @param function_1 Function values at the interpolation points on patch 1. 
     * @param function_2 Function values at the interpolation points on patch 2. 
     * @return the linear combination of the function values (c).
     * 
     * Remark: for the approximation case, the cells were selected by the given index 
     * ranges in the instantiation. Here, function_1 and function_2 can be given on the 
     * whole domain and the operator will select the correct values. But, it could 
     * be more optimised to directly give a slice on the selected cells. 
     * E.g.
     * get_function_coefficients(function_1[idx_interface_1][selected_cells_idx_range_1], 
     *                           function_2[idx_interface_2][selected_cells_idx_range_2]); 
     * and
     * get_function_coefficients(function_1[idx_interface_1], 
     *                           function_2[idx_interface_2]); 
     * will return the same value. 
     */
    template <class Layout1, class Layout2>
    double get_function_coefficients(
            DConstField<IdxRange1DPerp_1, Kokkos::HostSpace, Layout1> const& function_1,
            DConstField<IdxRange1DPerp_2, Kokkos::HostSpace, Layout2> const& function_2) const
    {
        // The function needs to be continuous at the interface.
        Idx1D_1 interface_idx_1
                = (m_extremity_1 == FRONT) ? m_idx_range_perp_1.front() : m_idx_range_perp_1.back();
        Idx1D_2 interface_idx_2
                = (m_extremity_2 == FRONT) ? m_idx_range_perp_2.front() : m_idx_range_perp_2.back();
        assert(abs(function_1(interface_idx_1) - function_2(interface_idx_2)) < 1e-13);

        double coeff_values = 0;
        ddc::for_each(m_idx_range_perp_1, [&](auto const& idx) {
            coeff_values += function_1(idx) * m_coeff_values_patch_1(idx);
        });

        // To avoid counting twice the value at the interface.
        IdxRange1DPerp_2 idx_range_perp_2_without_interface
                = (m_extremity_2 == FRONT)
                          ? m_idx_range_perp_2.remove_first(IdxStep<EdgePerpGrid2>(1))
                          : m_idx_range_perp_2.remove_last(IdxStep<EdgePerpGrid2>(1));
        ddc::for_each(idx_range_perp_2_without_interface, [&](auto const& idx) {
            coeff_values += function_2(idx) * m_coeff_values_patch_2(idx);
        });
        return coeff_values;
    }

    /**
     * @brief Get the linear combination of the function values (c).
     * See @ref get_function_coefficients.
     * @param function_1 Function values at the interpolation points on patch 1. 
     * @param function_2 Function values at the interpolation points on patch 2. 
     * @return the linear combination of the function values (c).
     */
    template <class Layout1, class Layout2>
    double get_function_coefficients(
            DConstField<IdxRange1DPerp_2, Kokkos::HostSpace, Layout2> const& function_2,
            DConstField<IdxRange1DPerp_1, Kokkos::HostSpace, Layout1> const& function_1) const
    {
        return this->get_function_coefficients(function_1, function_2);
    }

private:
    /**
     * @brief Compute the coefficients a, b and the weights omega applying the recursive
     * formula. 
     */
    void set_coefficients_non_uniform_case()
    {
        // Memory allocation ---------------------------------------------------------------------
        DFieldMem<IdxRange1DPerp_1, Kokkos::HostSpace> coeff_c0_1_alloc(m_idx_range_perp_1);
        DFieldMem<IdxRange1DPerp_2, Kokkos::HostSpace> coeff_c0_2_alloc(m_idx_range_perp_2);

        DFieldMem<IdxRange1DPerp_1, Kokkos::HostSpace> coeff_c1_1_alloc(m_idx_range_perp_1);
        DFieldMem<IdxRange1DPerp_2, Kokkos::HostSpace> coeff_c1_2_alloc(m_idx_range_perp_2);

        DFieldMem<IdxRange1DPerp_1, Kokkos::HostSpace> coeff_c2_1_alloc(m_idx_range_perp_1);
        DFieldMem<IdxRange1DPerp_2, Kokkos::HostSpace> coeff_c2_2_alloc(m_idx_range_perp_2);


        DField<IdxRange1DPerp_1, Kokkos::HostSpace> coeff_c0_1(coeff_c0_1_alloc);
        DField<IdxRange1DPerp_2, Kokkos::HostSpace> coeff_c0_2(coeff_c0_2_alloc);

        DField<IdxRange1DPerp_1, Kokkos::HostSpace> coeff_c1_1(coeff_c1_1_alloc);
        DField<IdxRange1DPerp_2, Kokkos::HostSpace> coeff_c1_2(coeff_c1_2_alloc);

        DField<IdxRange1DPerp_1, Kokkos::HostSpace> coeff_c2_1(coeff_c2_1_alloc);
        DField<IdxRange1DPerp_2, Kokkos::HostSpace> coeff_c2_2(coeff_c2_2_alloc);


        int const n_points_1 = m_idx_range_perp_1.size();
        int const n_points_2 = m_idx_range_perp_2.size();

        // Please provide at least 3 cells.
        assert(n_points_1 > 2);
        assert(n_points_2 > 2);

        // Initialise fields
        ddc::parallel_fill(coeff_c0_1, 0.0);
        ddc::parallel_fill(coeff_c0_2, 0.0);
        ddc::parallel_fill(coeff_c1_1, 0.0);
        ddc::parallel_fill(coeff_c1_2, 0.0);
        ddc::parallel_fill(coeff_c2_1, 0.0);
        ddc::parallel_fill(coeff_c2_2, 0.0);

        // Values at n-1 in the recursion
        double a0;
        double b0;

        // Values at n in the recursion
        double a1;
        double b1;

        // Values at n+1 in the recursion
        double a2;
        double b2;

        double alpha;
        double beta;
        std::array<double, 3> gammas;

        // Initialisation recursivity: set (a0,b0,c0) and (a1,b1,c1) -----------------------------
        Idx1D_1 interface_idx_1
                = (m_extremity_1 == FRONT) ? m_idx_range_perp_1.front() : m_idx_range_perp_1.back();
        Idx1D_2 interface_idx_2
                = (m_extremity_2 == FRONT) ? m_idx_range_perp_2.front() : m_idx_range_perp_2.back();

        // Set (a0,b0,c0) ---
        Idx1D_1 idx_1 = interface_idx_1;
        Idx1D_1 idx_1_minus;
        Idx1D_1 idx_1_plus = get_increment_idx(interface_idx_1, 1);

        Idx1D_2 idx_2 = interface_idx_2;
        Idx1D_2 idx_2_minus;
        Idx1D_2 idx_2_plus = get_increment_idx(interface_idx_2, 1);

        std::array<double, 2> cell_lengths = get_cell_lengths(idx_1);

        alpha = get_alpha(cell_lengths);
        beta = get_beta(cell_lengths);
        gammas = get_gammas(cell_lengths);

        a0 = alpha;
        b0 = beta;

        coeff_c0_1(idx_1_plus) = gammas[2];
        coeff_c0_1(idx_1) = gammas[1];
        coeff_c0_2(idx_2) = gammas[1];
        coeff_c0_2(idx_2_plus) = gammas[0];

        // Set (a1,b1,c1) ---
        idx_1 = get_increment_idx(interface_idx_1, 1);
        idx_1_minus = get_increment_idx(interface_idx_1, 0);
        idx_1_plus = get_increment_idx(interface_idx_1, 2);

        cell_lengths = get_cell_lengths(idx_1);

        alpha = get_alpha(cell_lengths);
        beta = get_beta(cell_lengths);
        gammas = get_gammas(cell_lengths);

        a1 = alpha * a0 / (1 - a0 * beta);
        b1 = b0 / (1 - a0 * beta);

        // c1 = (c0 + a0 * gamma) / (1 - a0 * beta);
        ddc::for_each(m_idx_range_perp_1, [&](auto const& idx) {
            coeff_c1_1(idx) = 1. / (1 - a0 * beta) * coeff_c0_1(idx);
        });
        ddc::for_each(m_idx_range_perp_2, [&](auto const& idx) {
            coeff_c1_2(idx) = 1. / (1 - a0 * beta) * coeff_c0_2(idx);
        });
        coeff_c1_1(idx_1_plus) += 1. / (1 - a0 * beta) * a0 * gammas[2];
        coeff_c1_1(idx_1) += 1. / (1 - a0 * beta) * a0 * gammas[1];
        coeff_c1_1(idx_1_minus) += 1. / (1 - a0 * beta) * a0 * gammas[0];
        coeff_c1_2(idx_2) += 1. / (1 - a0 * beta) * a0 * gammas[0];

        // Set (a2,b2,c2) ---
        std::tie(a2, b2) = std::tie(a1, b1);
        ddc::parallel_deepcopy(coeff_c2_1, get_const_field(coeff_c1_1));
        ddc::parallel_deepcopy(coeff_c2_2, get_const_field(coeff_c1_2));

        // Computing coefficients (c, a, b) on patch 1 -------------------------------------------
        // If is_cell_bound_1_with_extra_interpol_pt is true, the last cell is treated differently.
        int n1_cells = is_cell_bound_1_with_extra_interpol_pt ? n_points_1 - 3 : n_points_1 - 1;

        // for (int n(2); n < n1_cells; n++) {
        //     idx_1 = get_increment_idx(interface_idx_1, n);
        //     idx_1_minus = get_increment_idx(interface_idx_1, n - 1);
        //     idx_1_plus = get_increment_idx(interface_idx_1, n + 1);

        //     cell_lengths = get_cell_lengths(idx_1);

        //     alpha = get_alpha(cell_lengths);
        //     beta = get_beta(cell_lengths);
        //     gammas = get_gammas(cell_lengths);

        //     double denom = 1 - beta * a1 / a0;

        //     // Set (a2,b2,c2) ---
        //     a2 = alpha * a1 / denom;
        //     b2 = (b1 - beta * a1 / a0 * b0) / denom;

        //     ddc::for_each(m_idx_range_perp_1, [&](auto const& idx) {
        //         coeff_c2_1(idx) = 1. / denom * (coeff_c1_1(idx) - beta * a1 / a0 * coeff_c0_1(idx));
        //     });
        //     ddc::for_each(m_idx_range_perp_2, [&](auto const& idx) {
        //         coeff_c2_2(idx) = 1. / denom * (coeff_c1_2(idx) - beta * a1 / a0 * coeff_c0_2(idx));
        //     });
        //     coeff_c2_1(idx_1_plus) += 1. / denom * a1 * gammas[2];
        //     coeff_c2_1(idx_1) += 1. / denom * a1 * gammas[1];
        //     coeff_c2_1(idx_1_minus) += 1. / denom * a1 * gammas[0];


        //     // Update (a0,b0,c0) and (a1,b1,c1).
        //     std::tie(a0, b0) = std::tie(a1, b1);
        //     std::tie(a1, b1) = std::tie(a2, b2);

        //     ddc::parallel_deepcopy(coeff_c0_1, get_const_field(coeff_c1_1));
        //     ddc::parallel_deepcopy(coeff_c0_2, get_const_field(coeff_c1_2));

        //     ddc::parallel_deepcopy(coeff_c1_1, get_const_field(coeff_c2_1));
        //     ddc::parallel_deepcopy(coeff_c1_2, get_const_field(coeff_c2_2));
        // }
        std::array<double, 3> coeff_patch1 = {a0, a1, a2};
        std::array<double, 3> coeff_patch2 = {b0, b1, b2};
        std::array<DField<IdxRange1DPerp_1, Kokkos::HostSpace>, 3> coeff_c_patch1
                = {coeff_c0_1, coeff_c1_1, coeff_c2_1};
        std::array<DField<IdxRange1DPerp_2, Kokkos::HostSpace>, 3> coeff_c_patch2
                = {coeff_c0_2, coeff_c1_2, coeff_c2_2};
        recursion(
                n1_cells,
                interface_idx_1,
                coeff_patch1,
                coeff_patch2,
                coeff_c_patch1,
                coeff_c_patch2);
        a0 = coeff_patch1[0];
        a1 = coeff_patch1[1];
        a2 = coeff_patch1[2];
        b0 = coeff_patch2[0];
        b1 = coeff_patch2[1];
        b2 = coeff_patch2[2];
        coeff_c0_1 = coeff_c_patch1[0];
        coeff_c1_1 = coeff_c_patch1[1];
        coeff_c2_1 = coeff_c_patch1[2];
        coeff_c0_2 = coeff_c_patch2[0];
        coeff_c1_2 = coeff_c_patch2[1];
        coeff_c2_2 = coeff_c_patch2[2];

        // Correction if we use an additional interpolation point as closure in the last cell.
        if (is_cell_bound_1_with_extra_interpol_pt) {
            // int n = n_points_1 - 3;
            // idx_1_minus = get_increment_idx(interface_idx_1, n - 1);
            // idx_1 = get_increment_idx(interface_idx_1, n);
            // idx_1_plus = get_increment_idx(interface_idx_1, n + 2);
            // // Additional interpolation point.
            // Idx1D_1 idx_1_star = get_increment_idx(interface_idx_1, n + 1);

            // cell_lengths = get_cell_lengths(idx_1);
            // cell_lengths[1] += get_cell_lengths(idx_1_star)[1];

            // alpha = get_alpha(cell_lengths);
            // beta = get_beta(cell_lengths);
            // gammas = get_gammas(cell_lengths);

            // // Correction of gammas, alpha and beta
            // double const x = get_cell_lengths(idx_1_star)[1] / fabs(cell_lengths[1]);
            // double const H1 = (1 - x) * (1 - x) * (1 + 2 * x);
            // double const H0 = x * x * (3 - 2 * x);
            // double const K1 = (1 - x) * (1 - x) * x;
            // double const K0 = x * x * (x - 1);

            // beta = beta / (1 + alpha * K0 / K1);

            // double denom_gamma = 1 / (1 + alpha * K0 / K1);

            // double const gamma_plus = denom_gamma * (gammas[2] - alpha * H1 / K1 / cell_lengths[1]);
            // double const gamma_star = denom_gamma * alpha / K1 / cell_lengths[1];
            // double const gamma_idx = denom_gamma * (gammas[1] - alpha * H0 / K1 / cell_lengths[1]);
            // double const gamma_minus = denom_gamma * gammas[0];

            // alpha = 0;

            // double denom = 1 - beta * a1 / a0;

            // // Set (a2,b2,c2) ---
            // a2 = alpha * a1 / denom;
            // b2 = (b1 - beta * a1 / a0 * b0) / denom;

            // // c2 = (c1 + a1 * gamma - beta * a1 / a0 * c0) / denom;
            // ddc::for_each(m_idx_range_perp_1, [&](auto const& idx) {
            //     coeff_c2_1(idx) = 1. / denom * (coeff_c1_1(idx) - beta * a1 / a0 * coeff_c0_1(idx));
            // });
            // ddc::for_each(m_idx_range_perp_2, [&](auto const& idx) {
            //     coeff_c2_2(idx) = 1. / denom * (coeff_c1_2(idx) - beta * a1 / a0 * coeff_c0_2(idx));
            // });
            // coeff_c2_1(idx_1_plus) += 1. / denom * a1 * gamma_plus;
            // coeff_c2_1(idx_1) += 1. / denom * a1 * gamma_idx;
            // coeff_c2_1(idx_1_star) += 1. / denom * a1 * gamma_star;
            // coeff_c2_1(idx_1_minus) += 1. / denom * a1 * gamma_minus;

            std::array<double, 3> coeff_patch1 = {a0, a1, a2};
            std::array<double, 3> coeff_patch2 = {b0, b1, b2};
            std::array<DField<IdxRange1DPerp_1, Kokkos::HostSpace>, 3> coeff_c_patch1
                    = {coeff_c0_1, coeff_c1_1, coeff_c2_1};
            std::array<DField<IdxRange1DPerp_2, Kokkos::HostSpace>, 3> coeff_c_patch2
                    = {coeff_c0_2, coeff_c1_2, coeff_c2_2};
            correction_boundary(
                    n_points_1,
                    interface_idx_1,
                    coeff_patch1,
                    coeff_patch2,
                    coeff_c_patch1,
                    coeff_c_patch2);
            a0 = coeff_patch1[0];
            a1 = coeff_patch1[1];
            a2 = coeff_patch1[2];
            b0 = coeff_patch2[0];
            b1 = coeff_patch2[1];
            b2 = coeff_patch2[2];
            coeff_c0_1 = coeff_c_patch1[0];
            coeff_c1_1 = coeff_c_patch1[1];
            coeff_c2_1 = coeff_c_patch1[2];
            coeff_c0_2 = coeff_c_patch2[0];
            coeff_c1_2 = coeff_c_patch2[1];
            coeff_c2_2 = coeff_c_patch2[2];
        }


        // Initialisation second recursivity: set (a0,b0,c0) and (a1,b1,c1) ----------------------
        // Set (a0,b0,c0) ---
        idx_2_minus = get_increment_idx(interface_idx_2, 0);
        idx_2 = get_increment_idx(interface_idx_2, 1);
        idx_2_plus = get_increment_idx(interface_idx_2, 2);

        cell_lengths = get_cell_lengths(idx_2);

        alpha = get_alpha(cell_lengths);
        beta = get_beta(cell_lengths);
        gammas = get_gammas(cell_lengths);

        std::tie(a0, b0) = std::tie(a2, b2);
        ddc::parallel_deepcopy(coeff_c0_1, get_const_field(coeff_c2_1));
        ddc::parallel_deepcopy(coeff_c0_2, get_const_field(coeff_c2_2));

        // Set (a1,b1,c1) ---
        a1 = a0 / (1 - b0 * alpha);
        b1 = b0 * beta / (1 - b0 * alpha);

        // c1 = (c0 + b0 * gamma) / (1 - b0 * alpha);
        ddc::for_each(m_idx_range_perp_1, [&](auto const& idx) {
            coeff_c1_1(idx) = 1. / (1 - b0 * alpha) * coeff_c0_1(idx);
        });
        ddc::for_each(m_idx_range_perp_2, [&](auto const& idx) {
            coeff_c1_2(idx) = 1. / (1 - b0 * alpha) * coeff_c0_2(idx);
        });
        coeff_c1_2(idx_2_plus) += 1. / (1 - b0 * alpha) * b0 * gammas[0];
        coeff_c1_2(idx_2) += 1. / (1 - b0 * alpha) * b0 * gammas[1];
        coeff_c1_2(idx_2_minus) += 1. / (1 - b0 * alpha) * b0 * gammas[2];
        coeff_c1_1(interface_idx_1) += 1. / (1 - b0 * alpha) * b0 * gammas[2];

        // Computing coefficients (c, a, b) on patch 2 -------------------------------------------
        // If is_cell_bound_2_with_extra_interpol_pt is true, the last cell is treated differently.
        int n2_cells = is_cell_bound_2_with_extra_interpol_pt ? n_points_2 - 3 : n_points_2 - 1;

        // for (int n(2); n < n2_cells; n++) {
        //     idx_2_minus = get_increment_idx(interface_idx_2, n - 1);
        //     idx_2 = get_increment_idx(interface_idx_2, n);
        //     idx_2_plus = get_increment_idx(interface_idx_2, n + 1);

        //     cell_lengths = get_cell_lengths(idx_2);

        //     alpha = get_alpha(cell_lengths);
        //     beta = get_beta(cell_lengths);
        //     gammas = get_gammas(cell_lengths);

        //     double denom = 1 - alpha * b1 / b0;

        //     // Set (a2,b2,c2) ---
        //     a2 = (a1 - alpha * b1 / b0 * a0) / denom;
        //     b2 = b1 * beta / denom;

        //     // c2 = (c1 + b1 * gamma - alpha * b1 / b0 * c0) / denom;
        //     ddc::for_each(m_idx_range_perp_1, [&](auto const& idx) {
        //         coeff_c2_1(idx)
        //                 = 1. / denom * (coeff_c1_1(idx) - alpha * b1 / b0 * coeff_c0_1(idx));
        //     });
        //     ddc::for_each(m_idx_range_perp_2, [&](auto const& idx) {
        //         coeff_c2_2(idx)
        //                 = 1. / denom * (coeff_c1_2(idx) - alpha * b1 / b0 * coeff_c0_2(idx));
        //     });
        //     coeff_c2_2(idx_2_plus) += 1. / denom * b1 * gammas[0];
        //     coeff_c2_2(idx_2) += 1. / denom * b1 * gammas[1];
        //     coeff_c2_2(idx_2_minus) += 1. / denom * b1 * gammas[2];


        //     // Update (a0,b0,c0) and (a1,b1,c1).
        //     std::tie(a0, b0) = std::tie(a1, b1);
        //     std::tie(a1, b1) = std::tie(a2, b2);

        //     ddc::parallel_deepcopy(coeff_c0_1, get_const_field(coeff_c1_1));
        //     ddc::parallel_deepcopy(coeff_c0_2, get_const_field(coeff_c1_2));

        //     ddc::parallel_deepcopy(coeff_c1_1, get_const_field(coeff_c2_1));
        //     ddc::parallel_deepcopy(coeff_c1_2, get_const_field(coeff_c2_2));
        // }
        coeff_patch1 = {a0, a1, a2};
        coeff_patch2 = {b0, b1, b2};
        coeff_c_patch1 = {coeff_c0_1, coeff_c1_1, coeff_c2_1};
        coeff_c_patch2 = {coeff_c0_2, coeff_c1_2, coeff_c2_2};
        recursion(
                n2_cells,
                interface_idx_2,
                coeff_patch2,
                coeff_patch1,
                coeff_c_patch1,
                coeff_c_patch2);
        a0 = coeff_patch1[0];
        a1 = coeff_patch1[1];
        a2 = coeff_patch1[2];
        b0 = coeff_patch2[0];
        b1 = coeff_patch2[1];
        b2 = coeff_patch2[2];
        coeff_c0_1 = coeff_c_patch1[0];
        coeff_c1_1 = coeff_c_patch1[1];
        coeff_c2_1 = coeff_c_patch1[2];
        coeff_c0_2 = coeff_c_patch2[0];
        coeff_c1_2 = coeff_c_patch2[1];
        coeff_c2_2 = coeff_c_patch2[2];

        // Correction if we use an additional interpolation point as closure in the last cell.
        if (is_cell_bound_2_with_extra_interpol_pt) {
            // int n = n_points_2 - 3;
            // idx_2_minus = get_increment_idx(interface_idx_2, n - 1);
            // idx_2 = get_increment_idx(interface_idx_2, n);
            // idx_2_plus = get_increment_idx(interface_idx_2, n + 2);
            // // Additional interpolation point.
            // Idx1D_2 idx_2_star = get_increment_idx(interface_idx_2, n + 1);

            // cell_lengths = get_cell_lengths(idx_2);
            // cell_lengths[0] += get_cell_lengths(idx_2_star)[0];

            // alpha = get_alpha(cell_lengths);
            // beta = get_beta(cell_lengths);
            // gammas = get_gammas(cell_lengths);

            // // Correction of gammas, alpha and beta
            // double const x = get_cell_lengths(idx_2_star)[1] / fabs(cell_lengths[0]);
            // double const H1 = (1 - x) * (1 - x) * (1 + 2 * x);
            // double const H0 = x * x * (3 - 2 * x);
            // double const K1 = (1 - x) * (1 - x) * x;
            // double const K0 = x * x * (x - 1);

            // alpha = alpha / (1 + beta * K1 / K0);

            // double denom_gamma = 1 / (1 + beta * K1 / K0);

            // double const gamma_plus = denom_gamma * (gammas[0] - beta * H0 / K0 / cell_lengths[0]);
            // double const gamma_star = denom_gamma * beta / K0 / cell_lengths[0];
            // double const gamma_idx = denom_gamma * (gammas[1] - beta * H1 / K0 / cell_lengths[0]);
            // double const gamma_minus = denom_gamma * gammas[2];

            // beta = 0;

            // double denom = 1 - alpha * b1 / b0;

            // // Set (a2,b2,c2) ---
            // a2 = (a1 - alpha * b1 / b0 * a0) / denom;
            // b2 = b1 * beta / denom;

            // // c2 = (c1 + b1 * gamma - alpha * b1 / b0 * c0) / denom;
            // ddc::for_each(m_idx_range_perp_1, [&](auto const& idx) {
            //     coeff_c2_1(idx)
            //             = 1. / denom * (coeff_c1_1(idx) - alpha * b1 / b0 * coeff_c0_1(idx));
            // });
            // ddc::for_each(m_idx_range_perp_2, [&](auto const& idx) {
            //     coeff_c2_2(idx)
            //             = 1. / denom * (coeff_c1_2(idx) - alpha * b1 / b0 * coeff_c0_2(idx));
            // });
            // coeff_c2_2(idx_2_plus) += 1. / denom * b1 * gamma_plus;
            // coeff_c2_2(idx_2_star) += 1. / denom * b1 * gamma_star;
            // coeff_c2_2(idx_2) += 1. / denom * b1 * gamma_idx;
            // coeff_c2_2(idx_2_minus) += 1. / denom * b1 * gamma_minus;

            std::array<double, 3> coeff_patch1 = {a0, a1, a2};
            std::array<double, 3> coeff_patch2 = {b0, b1, b2};
            std::array<DField<IdxRange1DPerp_1, Kokkos::HostSpace>, 3> coeff_c_patch1
                    = {coeff_c0_1, coeff_c1_1, coeff_c2_1};
            std::array<DField<IdxRange1DPerp_2, Kokkos::HostSpace>, 3> coeff_c_patch2
                    = {coeff_c0_2, coeff_c1_2, coeff_c2_2};
            correction_boundary(
                    n_points_2,
                    interface_idx_2,
                    coeff_patch2,
                    coeff_patch1,
                    coeff_c_patch1,
                    coeff_c_patch2);
            a0 = coeff_patch1[0];
            a1 = coeff_patch1[1];
            a2 = coeff_patch1[2];
            b0 = coeff_patch2[0];
            b1 = coeff_patch2[1];
            b2 = coeff_patch2[2];
            coeff_c0_1 = coeff_c_patch1[0];
            coeff_c1_1 = coeff_c_patch1[1];
            coeff_c2_1 = coeff_c_patch1[2];
            coeff_c0_2 = coeff_c_patch2[0];
            coeff_c1_2 = coeff_c_patch2[1];
            coeff_c2_2 = coeff_c_patch2[2];
        }

        // Store the final values.
        m_coeff_deriv_patch_1 = a2;
        m_coeff_deriv_patch_2 = b2;
        ddc::parallel_deepcopy(m_coeff_values_patch_1, get_const_field(coeff_c2_1));
        ddc::parallel_deepcopy(m_coeff_values_patch_2, get_const_field(coeff_c2_2));
    }

    template <class Idx>
    void recursion(
            int const n_cells,
            Idx const interface_idx,
            std::array<double, 3>& coeff_patch,
            std::array<double, 3>& coeff_other_patch,
            std::array<DField<IdxRange1DPerp_1, Kokkos::HostSpace>, 3>& sum_patch1,
            std::array<DField<IdxRange1DPerp_2, Kokkos::HostSpace>, 3>& sum_patch2)
    {
        for (int n(2); n < n_cells; n++) {
            Idx const idx_minus = get_increment_idx(interface_idx, n - 1);
            Idx const idx = get_increment_idx(interface_idx, n);
            Idx const idx_plus = get_increment_idx(interface_idx, n + 1);

            std::array<double, 2> const cell_lengths = get_cell_lengths(idx);

            // Local coefficient further from the interface
            double const coeff_further = std::is_same_v<Idx, Idx1D_1> ? get_alpha(cell_lengths)
                                                                      : get_beta(cell_lengths);
            // Local coefficient closer from the interface
            double const coeff_closer = std::is_same_v<Idx, Idx1D_1> ? get_beta(cell_lengths)
                                                                     : get_alpha(cell_lengths);
            std::array<double, 3> const gammas = get_gammas(cell_lengths);

            double const denom = 1 - coeff_closer * coeff_patch[1] / coeff_patch[0];

            // Set (a2,b2,c2) ---
            coeff_patch[2] = coeff_further * coeff_patch[1] / denom;
            coeff_other_patch[2]
                    = (coeff_other_patch[1]
                       - coeff_closer * coeff_patch[1] / coeff_patch[0] * coeff_other_patch[0])
                      / denom;

            ddc::for_each(m_idx_range_perp_1, [&](auto const& idx_perp) {
                sum_patch1[2](idx_perp) = 1. / denom
                                          * (sum_patch1[1](idx_perp)
                                             - coeff_closer * coeff_patch[1] / coeff_patch[0]
                                                       * sum_patch1[0](idx_perp));
            });
            ddc::for_each(m_idx_range_perp_2, [&](auto const& idx_perp) {
                sum_patch2[2](idx_perp) = 1. / denom
                                          * (sum_patch2[1](idx_perp)
                                             - coeff_closer * coeff_patch[1] / coeff_patch[0]
                                                       * sum_patch2[0](idx_perp));
            });
            if constexpr (std::is_same_v<Idx, Idx1D_1>) {
                sum_patch1[2](idx_plus) += 1. / denom * coeff_patch[1] * gammas[2];
                sum_patch1[2](idx) += 1. / denom * coeff_patch[1] * gammas[1];
                sum_patch1[2](idx_minus) += 1. / denom * coeff_patch[1] * gammas[0];
            } else {
                sum_patch2[2](idx_plus) += 1. / denom * coeff_patch[1] * gammas[0];
                sum_patch2[2](idx) += 1. / denom * coeff_patch[1] * gammas[1];
                sum_patch2[2](idx_minus) += 1. / denom * coeff_patch[1] * gammas[2];
            }



            // Update coefficients.
            std::tie(coeff_patch[0], coeff_other_patch[0])
                    = std::tie(coeff_patch[1], coeff_other_patch[1]);
            std::tie(coeff_patch[1], coeff_other_patch[1])
                    = std::tie(coeff_patch[2], coeff_other_patch[2]);

            ddc::parallel_deepcopy(sum_patch1[0], get_const_field(sum_patch1[1]));
            ddc::parallel_deepcopy(sum_patch2[0], get_const_field(sum_patch2[1]));

            ddc::parallel_deepcopy(sum_patch1[1], get_const_field(sum_patch1[2]));
            ddc::parallel_deepcopy(sum_patch2[1], get_const_field(sum_patch2[2]));
        }
    }

    template <class Idx>
    void correction_boundary(
            int const n_points,
            Idx const interface_idx,
            std::array<double, 3>& coeff_patch,
            std::array<double, 3>& coeff_other_patch,
            std::array<DField<IdxRange1DPerp_1, Kokkos::HostSpace>, 3>& sum_patch1,
            std::array<DField<IdxRange1DPerp_2, Kokkos::HostSpace>, 3>& sum_patch2)
    {
        int const n = n_points - 3;
        Idx const idx_minus = get_increment_idx(interface_idx, n - 1);
        Idx const idx = get_increment_idx(interface_idx, n);
        Idx const idx_plus = get_increment_idx(interface_idx, n + 2);
        // Additional interpolation point.
        Idx const idx_star = get_increment_idx(interface_idx, n + 1);

        std::array<double, 2> cell_lengths = get_cell_lengths(idx);
        if constexpr (std::is_same_v<Idx, Idx1D_1>) {
            cell_lengths[1] += get_cell_lengths(idx_star)[1];
        } else {
            cell_lengths[0] += get_cell_lengths(idx_star)[0];
        }

        // Local coefficient further from the interface
        double coeff_further
                = std::is_same_v<Idx, Idx1D_1> ? get_alpha(cell_lengths) : get_beta(cell_lengths);
        // Local coefficient closer from the interface
        double coeff_closer
                = std::is_same_v<Idx, Idx1D_1> ? get_beta(cell_lengths) : get_alpha(cell_lengths);
        std::array<double, 3> const gammas = get_gammas(cell_lengths);

        // Correction of gammas, alpha and beta
        double x;
        if constexpr (std::is_same_v<Idx, Idx1D_1>) {
            x = get_cell_lengths(idx_star)[1] / fabs(cell_lengths[1]);
        } else {
            x = get_cell_lengths(idx_star)[1] / fabs(cell_lengths[0]);
        }
        double const H1 = (1 - x) * (1 - x) * (1 + 2 * x);
        double const H0 = x * x * (3 - 2 * x);
        double const K1 = (1 - x) * (1 - x) * x;
        double const K0 = x * x * (x - 1);

        double gamma_plus;
        double gamma_star;
        double gamma_idx;
        double gamma_minus;
        if constexpr (std::is_same_v<Idx, Idx1D_1>) {
            coeff_closer = coeff_closer / (1 + coeff_further * K0 / K1);

            double denom_gamma = 1 / (1 + coeff_further * K0 / K1);

            gamma_plus = denom_gamma * (gammas[2] - coeff_further * H1 / K1 / cell_lengths[1]);
            gamma_star = denom_gamma * coeff_further / K1 / cell_lengths[1];
            gamma_idx = denom_gamma * (gammas[1] - coeff_further * H0 / K1 / cell_lengths[1]);
            gamma_minus = denom_gamma * gammas[0];
        } else {
            coeff_closer = coeff_closer / (1 + coeff_further * K1 / K0);

            double denom_gamma = 1 / (1 + coeff_further * K1 / K0);

            gamma_plus = denom_gamma * (gammas[0] - coeff_further * H0 / K0 / cell_lengths[0]);
            gamma_star = denom_gamma * coeff_further / K0 / cell_lengths[0];
            gamma_idx = denom_gamma * (gammas[1] - coeff_further * H1 / K0 / cell_lengths[0]);
            gamma_minus = denom_gamma * gammas[2];
        }

        coeff_further = 0;

        double denom = 1 - coeff_closer * coeff_patch[1] / coeff_patch[0];

        // Set (coeff_patch[2],coeff_other_patch[2],c2) ---
        coeff_patch[2] = coeff_further * coeff_patch[1] / denom;
        coeff_other_patch[2]
                = (coeff_other_patch[1]
                   - coeff_closer * coeff_patch[1] / coeff_patch[0] * coeff_other_patch[0])
                  / denom;

        // c2 = (c1 + coeff_patch[1] * gamma - coeff_closer * coeff_patch[1] / coeff_patch[0] * c0) / denom;
        ddc::for_each(m_idx_range_perp_1, [&](auto const& idx_perp) {
            sum_patch1[2](idx_perp) = 1. / denom
                                      * (sum_patch1[1](idx_perp)
                                         - coeff_closer * coeff_patch[1] / coeff_patch[0]
                                                   * sum_patch1[0](idx_perp));
        });
        ddc::for_each(m_idx_range_perp_2, [&](auto const& idx_perp) {
            sum_patch2[2](idx_perp) = 1. / denom
                                      * (sum_patch2[1](idx_perp)
                                         - coeff_closer * coeff_patch[1] / coeff_patch[0]
                                                   * sum_patch2[0](idx_perp));
        });
        if constexpr (std::is_same_v<Idx, Idx1D_1>) {
            sum_patch1[2](idx_plus) += 1. / denom * coeff_patch[1] * gamma_plus;
            sum_patch1[2](idx) += 1. / denom * coeff_patch[1] * gamma_idx;
            sum_patch1[2](idx_star) += 1. / denom * coeff_patch[1] * gamma_star;
            sum_patch1[2](idx_minus) += 1. / denom * coeff_patch[1] * gamma_minus;
        } else {
            sum_patch2[2](idx_plus) += 1. / denom * coeff_patch[1] * gamma_plus;
            sum_patch2[2](idx_star) += 1. / denom * coeff_patch[1] * gamma_star;
            sum_patch2[2](idx) += 1. / denom * coeff_patch[1] * gamma_idx;
            sum_patch2[2](idx_minus) += 1. / denom * coeff_patch[1] * gamma_minus;
        }
    }

    /**
     * @brief Compute the coefficients a, b and the weights omega applying the explicit
     * formula for the uniform case. 
     */
    void set_coefficients_uniform_per_patch_case()
    {
        int const n_cells_1 = m_idx_range_perp_1.size() - 1;
        int const n_cells_2 = m_idx_range_perp_2.size() - 1;

        Idx1D_1 interface_idx_1
                = (m_extremity_1 == FRONT) ? m_idx_range_perp_1.front() : m_idx_range_perp_1.back();
        Idx1D_2 interface_idx_2
                = (m_extremity_2 == FRONT) ? m_idx_range_perp_2.front() : m_idx_range_perp_2.back();

        Idx1D_1 idx_1 = interface_idx_1;
        Idx1D_2 idx_2 = interface_idx_2;

        std::array<double, 2> cell_lengths = get_cell_lengths(idx_1);

        double const a_11 = -0.5 * cell_lengths[0] / (cell_lengths[0] + cell_lengths[1]);
        double const b_11 = -0.5 * cell_lengths[1] / (cell_lengths[0] + cell_lengths[1]);

        double const v1 = 2 * Kokkos::sqrt(3);

        double const vn1 = Kokkos::pow(2 + Kokkos::sqrt(3), n_cells_1)
                           - Kokkos::pow(2 - Kokkos::sqrt(3), n_cells_1);
        double const vn1_minus = Kokkos::pow(2 + Kokkos::sqrt(3), n_cells_1 - 1)
                                 - Kokkos::pow(2 - Kokkos::sqrt(3), n_cells_1 - 1);

        double const vn2 = Kokkos::pow(2 + Kokkos::sqrt(3), n_cells_2)
                           - Kokkos::pow(2 - Kokkos::sqrt(3), n_cells_2);
        double const vn2_minus = Kokkos::pow(2 + Kokkos::sqrt(3), n_cells_2 - 1)
                                 - Kokkos::pow(2 - Kokkos::sqrt(3), n_cells_2 - 1);

        double const denominator = vn1 * vn2 + vn1_minus * vn2 * a_11 + vn1 * vn2_minus * b_11;


        // Compute the coefficients a and b.
        m_coeff_deriv_patch_1 = Kokkos::pow(-1, n_cells_1 - 1) * v1 * a_11 * vn2 / denominator;
        m_coeff_deriv_patch_2 = Kokkos::pow(-1, n_cells_2 - 1) * v1 * b_11 * vn1 / denominator;


        // Compute the weights {omega_k}_{k = -Nl, ..., Nr} in c.
        // --- for k = 0
        m_coeff_values_patch_1(idx_1) = -3 / denominator
                                        * (a_11 / cell_lengths[1] * vn2 * (vn1 - vn1_minus)
                                           - b_11 / cell_lengths[0] * vn1 * (vn2 - vn2_minus));
        m_coeff_values_patch_2(idx_2) = -3 / denominator
                                        * (a_11 / cell_lengths[1] * vn2 * (vn1 - vn1_minus)
                                           - b_11 / cell_lengths[0] * vn1 * (vn2 - vn2_minus));

        // --- for k = 1, ..., Nr-1
        for (int k(1); k < n_cells_1; ++k) {
            idx_1 = get_increment_idx(idx_1, 1);

            const int n1_minus_k = n_cells_1 - k;
            double const vk_minus = Kokkos::pow(2 + Kokkos::sqrt(3), n1_minus_k - 1)
                                    - Kokkos::pow(2 - Kokkos::sqrt(3), n1_minus_k - 1);
            double const vk_plus = Kokkos::pow(2 + Kokkos::sqrt(3), n1_minus_k + 1)
                                   - Kokkos::pow(2 - Kokkos::sqrt(3), n1_minus_k + 1);

            m_coeff_values_patch_1(idx_1) = Kokkos::pow(-1, k + 1) * 3 / denominator * a_11
                                            / cell_lengths[1] * vn2 * (vk_plus - vk_minus);
        }

        // --- for k = Nr
        idx_1 = get_increment_idx(idx_1, 1);
        m_coeff_values_patch_1(idx_1) = Kokkos::pow(-1, n_cells_1 + 1) * 3 / denominator * a_11
                                        / cell_lengths[1] * vn2 * v1;

        // --- for k = -1, ..., -(Nl-1)
        for (int k(1); k < n_cells_2; ++k) {
            idx_2 = get_increment_idx(idx_2, 1);

            const int n2_minus_k = n_cells_2 - k;
            double const vk_minus = Kokkos::pow(2 + Kokkos::sqrt(3), n2_minus_k - 1)
                                    - Kokkos::pow(2 - Kokkos::sqrt(3), n2_minus_k - 1);
            double const vk_plus = Kokkos::pow(2 + Kokkos::sqrt(3), n2_minus_k + 1)
                                   - Kokkos::pow(2 - Kokkos::sqrt(3), n2_minus_k + 1);

            m_coeff_values_patch_2(idx_2) = Kokkos::pow(-1, k) * 3 / denominator * b_11
                                            / cell_lengths[0] * vn1 * (vk_plus - vk_minus);
        }

        // --- for k = -Nl
        idx_2 = get_increment_idx(idx_2, 1);
        m_coeff_values_patch_2(idx_2)
                = Kokkos::pow(-1, n_cells_2) * 3 / denominator * b_11 / cell_lengths[0] * vn1 * v1;
    };


    /**
     * @brief Get the lengths of the cell of in the right side and the left side of the given index.
     */
    template <class Idx1D>
    std::array<double, 2> get_cell_lengths(Idx1D const& idx) const
    {
        static_assert(
                std::is_same_v<Idx1D, Idx1D_1> || std::is_same_v<Idx1D, Idx1D_2>,
                "Wrong type of index given.");

        // Tag for the given patch.
        using IdxRange1D = std::
                conditional_t<std::is_same_v<Idx1D, Idx1D_1>, IdxRange1DPerp_1, IdxRange1DPerp_2>;

        // Tags for the other patch.
        using OIdx1D = std::conditional_t<std::is_same_v<Idx1D, Idx1D_1>, Idx1D_2, Idx1D_1>;
        using OIdxRange1D = std::
                conditional_t<std::is_same_v<Idx1D, Idx1D_1>, IdxRange1DPerp_2, IdxRange1DPerp_1>;

        using GlobalIdxRange = IdxRange<EdgePerpGrid1, EdgePerpGrid2>;

        Extremity extremity = std::is_same_v<Idx1D, Idx1D_1> ? m_extremity_1 : m_extremity_2;
        Extremity other_extremity = std::is_same_v<Idx1D, Idx1D_1> ? m_extremity_2 : m_extremity_1;

        IdxRange1D idx_range_1d(GlobalIdxRange(m_idx_range_perp_1, m_idx_range_perp_2));
        OIdxRange1D other_idx_range_1d(GlobalIdxRange(m_idx_range_perp_1, m_idx_range_perp_2));

        // Please, do not provide a index on the boundary except for the interface.
        assert(!(idx == idx_range_1d.front() && extremity == BACK));
        assert(!(idx == idx_range_1d.back() && extremity == FRONT));

        std::array<double, 2> cell_lengths;

        Idx1D const idx_minus((idx - Idx1D(1)).value());
        Idx1D const idx_plus((idx - Idx1D(-1)).value());

        double const delta_coord_plus = abs(ddc::coordinate(idx_plus) - ddc::coordinate(idx));
        double const delta_coord_minus = abs(ddc::coordinate(idx) - ddc::coordinate(idx_minus));

        bool const is_same_orientation_as_global
                = (std::is_same_v<Idx1D, Idx1D_1> && m_extremity_1 == FRONT)
                  || (std::is_same_v<Idx1D, Idx1D_2> && m_extremity_2 == BACK);

        // If the given index corresponds to the index at the interface,
        Idx1D const idx_interface
                = (extremity == FRONT) ? idx_range_1d.front() : idx_range_1d.back();
        if (idx == idx_interface) {
            Idx1D const idx_incremented = get_increment_idx(idx, 1);

            OIdx1D const other_idx = (other_extremity == FRONT) ? other_idx_range_1d.front()
                                                                : other_idx_range_1d.back();
            OIdx1D const other_idx_incremented = get_increment_idx(other_idx, 1);

            cell_lengths[0]
                    = abs(ddc::coordinate(other_idx_incremented) - ddc::coordinate(other_idx));
            cell_lengths[1] = abs(ddc::coordinate(idx_incremented) - ddc::coordinate(idx));
        }
        // If given index is on the patch 1 or on the patch 2,
        else if (is_same_orientation_as_global) {
            // . | →   or  → | .
            cell_lengths[0] = delta_coord_minus;
            cell_lengths[1] = delta_coord_plus;
        } else {
            //  . | ←  or  ← | .
            cell_lengths[0] = delta_coord_plus;
            cell_lengths[1] = delta_coord_minus;
        }

        if ((cell_lengths[0] == 0) || (cell_lengths[1] == 0)) {
            Kokkos::abort("[abort] Ill-defined: the length of the cells in"
                          "the cell_lengths array must be not zero.");
        }
        return cell_lengths;
    }


    /// @brief Compute the coefficients gammas.
    std::array<double, 3> get_gammas(std::array<double, 2> const& cell_lengths) const
    {
        double const length_left = cell_lengths[0];
        double const length_right = cell_lengths[1];

        double gamma_1 = 3. / 2. * 1. / (length_left + length_right) * length_right / length_left;
        double gamma_2 = 3. / 2. * 1. / (length_left + length_right)
                         * (length_left / length_right - length_right / length_left);
        double gamma_3 = -3. / 2. * 1. / (length_left + length_right) * length_left / length_right;

        return std::array<double, 3>({gamma_1, gamma_2, gamma_3});
    }

    /// @brief Compute the coefficient alpha.
    double get_alpha(std::array<double, 2> const& cell_lengths) const
    {
        double const length_left = cell_lengths[0];
        double const length_right = cell_lengths[1];
        return -0.5 * length_left / (length_right + length_left);
    }

    /// @brief Compute the coefficient beta.
    double get_beta(std::array<double, 2> const& cell_lengths) const
    {
        double const length_left = cell_lengths[0];
        double const length_right = cell_lengths[1];
        return -0.5 * length_right / (length_right + length_left);
    }


    /**
     * @brief Increment a given index.
     * We can to go far from the Interface. 
     * So, the increment increases the value of the index if the Interface is 
     * at the .front() of the index range.
     * E.g.
     *     0                  N > 0
     *     | ------------------>
     * Interface     increment(idx) = idx+1 
     * 
     * The increment decreases the value of the index if the Interface is 
     * at the .back() of the index range.
     * E.g.
     *    N > 0                  0
     *     | ------------------>
     * Interface     increment(idx) = idx-1 
     * 
     * The operator can add/sustract to "idx_0" an integer "increment".
     */
    template <class Idx>
    Idx get_increment_idx(Idx const& idx_0, int const increment) const
    {
        Extremity extremity = std::is_same_v<Idx, Idx1D_1> ? m_extremity_1 : m_extremity_2;

        Idx idx_incremented;
        if (extremity == FRONT) {
            return idx_incremented = Idx((idx_0 - Idx(-increment)).value());
        } else {
            return idx_incremented = Idx((idx_0 - Idx(increment)).value());
        }
    }
};