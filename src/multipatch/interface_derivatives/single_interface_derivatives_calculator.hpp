// SPDX-License-Identifier: MIT

#pragma once

#include <ddc/ddc.hpp>

#include "ddc_aliases.hpp"
#include "edge.hpp"
#include "edge_transformation.hpp"
#include "types.hpp"
#include "geometry_descriptors.hpp"


template <class T>
inline constexpr bool enable_single_derivative_calculator = false;

template <class T>
inline constexpr bool is_single_derivative_calculator_v
        = enable_single_derivative_calculator<std::remove_const_t<std::remove_reference_t<T>>>;

/**
 * @brief Compute the derivative of an equivalent global spline 
 * at the interface between two patches. 
 * 
 * For a given Interface, this operator computes the coefficients 
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
 *      The weights in the linear combination depend on the meshes of the two patches. 
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
 *     |   1   |   2   |
 *     |_______|_______|
 *          + ← → +
 *          - → ← -
 * Orientation of an equivalent global mesh : →
 * 
 * All the formulae and more details are given in the README.md. 
 * 
 * @tparam Interface The interface between two patches where we want 
 * to compute the derivatives.  
 * 
 * @warning The applied method only works for interpolation points located on 
 * the break points. In the case where "ddc::BoundCond::GREVILLE" is specified,
 * additional interpolation points are also placed in the first or last cell
 * of the patch.
 * Please be sure to initialise the discrete space of your Grid on the break points
 * (especially in the non-uniform case). 
 */
template <class InterfaceType>
class SingleInterfaceDerivativesCalculator
{
    static_assert(
            (!std::is_same_v<typename InterfaceType::Edge1, OutsideEdge>)&&(
                    !std::is_same_v<typename InterfaceType::Edge2, OutsideEdge>),
            "The interface cannot be an interface with the outside domain.");

    using EdgePerpGrid1 = typename InterfaceType::Edge1::perpendicular_grid;
    using EdgePerpGrid2 = typename InterfaceType::Edge2::perpendicular_grid;

    using EdgeParGrid1 = typename InterfaceType::Edge1::parallel_grid;
    using EdgeParGrid2 = typename InterfaceType::Edge2::parallel_grid;

    using Patch1 = typename InterfaceType::Edge1::associated_patch;
    using Patch2 = typename InterfaceType::Edge2::associated_patch;

    using IdxRange2D_1 = typename Patch1::IdxRange12;
    using IdxRange2D_2 = typename Patch2::IdxRange12;

    using IdxRange1DPerp_1 = IdxRange<EdgePerpGrid1>;
    using IdxRange1DPerp_2 = IdxRange<EdgePerpGrid2>;
    using IdxRange1DPar_1 = IdxRange<EdgeParGrid1>;
    using IdxRange1DPar_2 = IdxRange<EdgeParGrid2>;

    using Idx1D_1 = Idx<EdgePerpGrid1>;
    using Idx1D_2 = Idx<EdgePerpGrid2>;

    using BSplinesPerp1 = std::conditional_t<
            std::is_same_v<EdgePerpGrid1, typename Patch1::Grid1>,
            typename Patch1::BSplines1,
            typename Patch1::BSplines2>;
    using BSplinesPerp2 = std::conditional_t<
            std::is_same_v<EdgePerpGrid2, typename Patch2::Grid1>,
            typename Patch2::BSplines1,
            typename Patch2::BSplines2>;

    using GridBreakPt1 = std::conditional_t<
            BSplinesPerp1::is_uniform(),
            ddc::UniformBsplinesKnots<BSplinesPerp1>,
            ddc::NonUniformBsplinesKnots<BSplinesPerp1>>;
    using GridBreakPt2 = std::conditional_t<
            BSplinesPerp2::is_uniform(),
            ddc::UniformBsplinesKnots<BSplinesPerp2>,
            ddc::NonUniformBsplinesKnots<BSplinesPerp2>>;
    using IdxRangeBreakPt1 = IdxRange<GridBreakPt1>;
    using IdxRangeBreakPt2 = IdxRange<GridBreakPt2>;

public:
    /// @brief Interface between the two involved patches.
    using associated_interface = InterfaceType;

private:
    const bool m_is_cell_bound_1_with_extra_interpol_pt;
    const bool m_is_cell_bound_2_with_extra_interpol_pt;

    IdxRange1DPerp_1 const m_idx_range_perp_1;
    IdxRange1DPerp_2 const m_idx_range_perp_2;

    static Extremity constexpr m_extremity_1 = InterfaceType::Edge1::extremity;
    static Extremity constexpr m_extremity_2 = InterfaceType::Edge2::extremity;

    double m_coeff_deriv_patch_2;
    double m_coeff_deriv_patch_1;

    host_t<DFieldMem<IdxRange1DPerp_1>> m_weights_patch_1_alloc;
    host_t<DFieldMem<IdxRange1DPerp_2>> m_weights_patch_2_alloc;

    host_t<DField<IdxRange1DPerp_1>> m_weights_patch_1;
    host_t<DField<IdxRange1DPerp_2>> m_weights_patch_2;

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
     * If the interpolation points are uniform, it computes the coefficients
     * with an explicit formula. Otherwise, it uses the recursive formula 
     * (see README.md).
     *  
     * @param idx_range_1d_1 1D index range perpendicular to the Interface, 
     * on the patch 1. 
     * @param idx_range_1d_2 1D index range perpendicular to the Interface, 
     * on the patch 2. 
     * @param Bound1 The boundary condition type on the opposite edge of the interface on 
     * the patch 1. By default, the value is set to ddc::BoundCond::HERMITE. If 
     * ddc::BoundCond::GREVILLE is given, a treatment will be applied to consider the 
     * additional interpolation point. Giving ddc::BoundCond::PERIODIC does not make sense. 
     * @param Bound2 The boundary condition type on the opposite edge of the interface on 
     * the patch 2. By default, the value is set to ddc::BoundCond::HERMITE. If 
     * ddc::BoundCond::GREVILLE is given, a treatment will be applied to consider the 
     * additional interpolation point. Giving ddc::BoundCond::PERIODIC does not make sense.  
     */
    SingleInterfaceDerivativesCalculator(
            IdxRange1DPerp_1 const& idx_range_1d_1,
            IdxRange1DPerp_2 const& idx_range_1d_2,
            ddc::BoundCond const& Bound1 = ddc::BoundCond::HERMITE,
            ddc::BoundCond const& Bound2 = ddc::BoundCond::HERMITE)
        : m_is_cell_bound_1_with_extra_interpol_pt(Bound1 == ddc::BoundCond::GREVILLE)
        , m_is_cell_bound_2_with_extra_interpol_pt(Bound2 == ddc::BoundCond::GREVILLE)
        , m_idx_range_perp_1(idx_range_1d_1)
        , m_idx_range_perp_2(idx_range_1d_2)
        , m_weights_patch_1_alloc(m_idx_range_perp_1)
        , m_weights_patch_2_alloc(m_idx_range_perp_2)
        , m_weights_patch_1(m_weights_patch_1_alloc)
        , m_weights_patch_2(m_weights_patch_2_alloc)
    {
        assert(Bound1 != ddc::BoundCond::PERIODIC);
        assert(Bound2 != ddc::BoundCond::PERIODIC);

        // Two interpolation points have to be added if the derivatives are not closure condition.
        if (m_is_cell_bound_1_with_extra_interpol_pt) {
            assert(m_idx_range_perp_1.size() == ddc::discrete_space<BSplinesPerp1>().ncells() + 2);
        } else {
            assert(m_idx_range_perp_1.size() <= ddc::discrete_space<BSplinesPerp1>().ncells() + 1);
        }
        if (m_is_cell_bound_2_with_extra_interpol_pt) {
            assert(m_idx_range_perp_2.size() == ddc::discrete_space<BSplinesPerp2>().ncells() + 2);
        } else {
            assert(m_idx_range_perp_2.size() <= ddc::discrete_space<BSplinesPerp2>().ncells() + 1);
        }

        // The break points have to be interpolation points.
        check_break_points_are_interpolation_points<BSplinesPerp1, GridBreakPt1>(
                m_is_cell_bound_1_with_extra_interpol_pt,
                m_idx_range_perp_1,
                m_extremity_1);
        check_break_points_are_interpolation_points<BSplinesPerp2, GridBreakPt2>(
                m_is_cell_bound_2_with_extra_interpol_pt,
                m_idx_range_perp_2,
                m_extremity_2);


        // The additional interpolation points have to be in the boundary cells.
        if (m_is_cell_bound_1_with_extra_interpol_pt) {
            check_additional_interpolation_points_location<
                    BSplinesPerp1,
                    GridBreakPt1>(m_idx_range_perp_1, m_extremity_1);
        }
        if (m_is_cell_bound_2_with_extra_interpol_pt) {
            check_additional_interpolation_points_location<
                    BSplinesPerp2,
                    GridBreakPt2>(m_idx_range_perp_2, m_extremity_2);
        }

        if constexpr ((ddc::is_uniform_point_sampling_v<EdgePerpGrid1>)&&(
                              ddc::is_uniform_point_sampling_v<EdgePerpGrid2>)) {
            set_coefficients_uniform_per_patch_case();
        } else {
            set_coefficients_non_uniform_case();
        }
    }


    /**
     * @brief Instantiate SingleInterfaceDerivativesCalculator. 
     * See @ref SingleInterfaceDerivativesCalculatorInstantiator.
     * @param idx_range_a Index range on one patch. 
     * @param idx_range_b Index range on the other patch. 
     * @param Bound1 The boundary condition type on the opposite edge of the interface on 
     * the patch 1. By default, the value is set to ddc::BoundCond::HERMITE. 
     * @param Bound2 The boundary condition type on the opposite edge of the interface on 
     * the patch 2. By default, the value is set to ddc::BoundCond::HERMITE. 
     */
    template <class IdxRangeA, class IdxRangeB>
    SingleInterfaceDerivativesCalculator(
            IdxRangeA const& idx_range_a,
            IdxRangeB const& idx_range_b,
            ddc::BoundCond const& Bound1 = ddc::BoundCond::HERMITE,
            ddc::BoundCond const& Bound2 = ddc::BoundCond::HERMITE)
        : SingleInterfaceDerivativesCalculator(
                IdxRange1DPerp_1(
                        ddc::cartesian_prod_t<IdxRangeA, IdxRangeB>(idx_range_a, idx_range_b)),
                IdxRange1DPerp_2(
                        ddc::cartesian_prod_t<IdxRangeA, IdxRangeB>(idx_range_a, idx_range_b)),
                Bound1,
                Bound2)
    {
    }

    SingleInterfaceDerivativesCalculator(
            IdxRange1DPerp_1 const& idx_range_1d_1,
            IdxRange1DPerp_2 const& idx_range_1d_2,
            int const number_taken_cells)
        : SingleInterfaceDerivativesCalculator(
                (m_extremity_1 == Extremiy::FRONT)
                        ? idx_range_1d_1.take_first(IdxStep<EdgePerpGrid1>(number_taken_cells + 1))
                        : idx_range_1d_1.take_last(IdxStep<EdgePerpGrid1>(number_taken_cells + 1)),
                (m_extremity_2 == Extremiy::FRONT)
                        ? idx_range_1d_2.take_first(IdxStep<EdgePerpGrid2>(number_taken_cells + 1))
                        : idx_range_1d_2.take_last(IdxStep<EdgePerpGrid2>(number_taken_cells + 1)),
                ddc::BoundCond::HERMITE,
                ddc::BoundCond::HERMITE)
    {
    }

    template <class IdxRangeA, class IdxRangeB>
    SingleInterfaceDerivativesCalculator(
            IdxRangeA const& idx_range_a,
            IdxRangeB const& idx_range_b,
            int const number_taken_cells)
        : SingleInterfaceDerivativesCalculator(
                IdxRange1DPerp_1(
                        ddc::cartesian_prod_t<IdxRangeA, IdxRangeB>(idx_range_a, idx_range_b)),
                IdxRange1DPerp_2(
                        ddc::cartesian_prod_t<IdxRangeA, IdxRangeB>(idx_range_a, idx_range_b)),
                number_taken_cells)
    {
    }

    static Extremity constexpr m_extremity_1 = InterfaceType::Edge1::extremity;
    static Extremity constexpr m_extremity_2 = InterfaceType::Edge2::extremity;

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
     * @brief Get the coefficient (a) or (b) in front of the derivative on the given patch of the given Interface.
     * @return The value of the coefficient (a) or (b). 
     */
    template <class Patch>
    double get_coeff_deriv_on_patch() const
    {
        static_assert(
                (std::is_same_v<Patch, Patch1>) || (std::is_same_v<Patch, Patch2>),
                "The given Patch template parameter is not one of the Patch of the given "
                "Interface.");
        if constexpr (std::is_same_v<Patch, Patch1>) {
            return m_coeff_deriv_patch_1;
        } else {
            return m_coeff_deriv_patch_2;
        }
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
        Idx1D_1 interface_idx_1 = get_extremity_idx(m_extremity_1, m_idx_range_perp_1);
        Idx1D_2 interface_idx_2 = get_extremity_idx(m_extremity_2, m_idx_range_perp_2);
        assert(abs(function_1(interface_idx_1) - function_2(interface_idx_2)) < 1e-13);

        double coeff_values = ddc::host_transform_reduce(
                m_idx_range_perp_1,
                0.0,
                ddc::reducer::sum<double>(),
                [&](Idx1D_1 const& idx) { return function_1(idx) * m_weights_patch_1(idx); });

        // To avoid counting twice the value at the interface.
        IdxRange1DPerp_2 idx_range_perp_2_without_interface
                = (m_extremity_2 == FRONT)
                          ? m_idx_range_perp_2.remove_first(IdxStep<EdgePerpGrid2>(1))
                          : m_idx_range_perp_2.remove_last(IdxStep<EdgePerpGrid2>(1));
        coeff_values += ddc::host_transform_reduce(
                idx_range_perp_2_without_interface,
                0.0,
                ddc::reducer::sum<double>(),
                [&](Idx1D_2 const& idx) { return function_2(idx) * m_weights_patch_2(idx); });
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
    inline double get_function_coefficients(
            DConstField<IdxRange1DPerp_2, Kokkos::HostSpace, Layout2> const& function_2,
            DConstField<IdxRange1DPerp_1, Kokkos::HostSpace, Layout1> const& function_1) const
    {
        return get_function_coefficients(function_1, function_2);
    }


private:
    template <typename BSplinesPerp, typename GridBreakPt, typename EdgePerpGrid>
    void check_break_points_are_interpolation_points(
            bool const is_cell_bound_with_extra_interpol_pt,
            IdxRange<EdgePerpGrid> const& idx_range_perp,
            Extremity const extremity) const
    {
        // Get the index range for all the break points on the patch.
        IdxRange<GridBreakPt> idx_range_full_break_points
                = ddc::discrete_space<BSplinesPerp>().break_point_domain();
        IdxRange<GridBreakPt> idx_range_break_points = idx_range_full_break_points;
        // Select the break points where the interpolation points are defined.
        if (!is_cell_bound_with_extra_interpol_pt) {
            IdxStep<GridBreakPt> idx_step(idx_range_perp.size());
            if (extremity == Extremity::FRONT) {
                idx_range_break_points = idx_range_full_break_points.take_first(idx_step);
            } else {
                idx_range_break_points = idx_range_full_break_points.take_last(idx_step);
            }
        }
        // Each break point should be an interpolation point.
        ddc::host_for_each(idx_range_break_points, [&](Idx<GridBreakPt> const idx_break) {
            double const break_point = ddc::coordinate(idx_break);
            bool is_an_interpolation_pt = false;
            ddc::host_for_each(idx_range_perp, [&](Idx<EdgePerpGrid> const idx_interpol) {
                double const interpolation_point = ddc::coordinate(idx_interpol);
                if (abs(break_point - interpolation_point) < 1e-15) {
                    is_an_interpolation_pt = true;
                }
            });
            if (!is_an_interpolation_pt) {
                throw std::runtime_error(
                        "[abort] The break points have to be interpolation points. The break point "
                        + std::to_string(break_point) + " at the "
                        + std::to_string((idx_break - idx_range_break_points.front()).value())
                        + "th index was not found in the interpolation point grid given.");
            }
        });
    }

    template <typename BSplinesPerp, typename GridBreakPt, typename EdgePerpGrid>
    void check_additional_interpolation_points_location(
            IdxRange<EdgePerpGrid> const& idx_range_perp,
            Extremity const extremity)
    {
        IdxRange<GridBreakPt> idx_range_full_break_points
                = ddc::discrete_space<BSplinesPerp>().break_point_domain();

        double added_interpolation_pt;
        double break_coord_min;
        double break_coord_max;
        if (extremity == Extremity::FRONT) {
            Idx<EdgePerpGrid> idx = idx_range_perp.back() - IdxStep<EdgePerpGrid>(1);
            Idx<GridBreakPt> idx_break
                    = idx_range_full_break_points.back() - IdxStep<GridBreakPt>(1);
            added_interpolation_pt = double(ddc::coordinate(idx));
            break_coord_min = double(ddc::coordinate(idx_break));
            break_coord_max = double(ddc::coordinate(idx_range_full_break_points.back()));
        } else {
            Idx<EdgePerpGrid> idx = idx_range_perp.front() + IdxStep<EdgePerpGrid>(1);
            Idx<GridBreakPt> idx_break
                    = idx_range_full_break_points.front() + IdxStep<GridBreakPt>(1);
            added_interpolation_pt = double(ddc::coordinate(idx));
            break_coord_min = double(ddc::coordinate(idx_range_full_break_points.front()));
            break_coord_max = double(ddc::coordinate(idx_break));
        }

        if (!((break_coord_min < added_interpolation_pt)
              && (added_interpolation_pt < break_coord_max))) {
            throw std::runtime_error(
                    "[abort] The additional interpolation points have to be placed in the "
                    "first or last cell of the patch. The point "
                    + std::to_string(added_interpolation_pt) + " is not placed between "
                    + std::to_string(break_coord_min) + " and " + std::to_string(break_coord_max)
                    + ".");
        }
    }

    // ===========================================================================================
    // NON-UNIFORM INTERPOLATION POINTS                                                          |
    // ===========================================================================================
    /**
     * @brief Compute the coefficients a, b and the weights omega applying the recursive
     * formula. 
     */
    void set_coefficients_non_uniform_case()
    {
        // Memory allocation ---------------------------------------------------------------------
        // Define weight fields for number of cells equal to n, n-1 and n-2 in the recursion.
        host_t<DFieldMem<IdxRange1DPerp_1>> weights_1_alloc_minus_2(m_idx_range_perp_1);
        host_t<DFieldMem<IdxRange1DPerp_1>> weights_1_alloc_minus_1(m_idx_range_perp_1);
        host_t<DFieldMem<IdxRange1DPerp_1>> weights_1_alloc_minus_0(m_idx_range_perp_1);

        host_t<DFieldMem<IdxRange1DPerp_2>> weights_2_alloc_minus_2(m_idx_range_perp_2);
        host_t<DFieldMem<IdxRange1DPerp_2>> weights_2_alloc_minus_1(m_idx_range_perp_2);
        host_t<DFieldMem<IdxRange1DPerp_2>> weights_2_alloc_minus_0(m_idx_range_perp_2);

        host_t<DField<IdxRange1DPerp_1>> weights_1_minus_2(weights_1_alloc_minus_2);
        host_t<DField<IdxRange1DPerp_1>> weights_1_minus_1(weights_1_alloc_minus_1);
        host_t<DField<IdxRange1DPerp_1>> weights_1_minus_0(weights_1_alloc_minus_0);

        host_t<DField<IdxRange1DPerp_2>> weights_2_minus_2(weights_2_alloc_minus_2);
        host_t<DField<IdxRange1DPerp_2>> weights_2_minus_1(weights_2_alloc_minus_1);
        host_t<DField<IdxRange1DPerp_2>> weights_2_minus_0(weights_2_alloc_minus_0);

        int const n_points_1 = m_idx_range_perp_1.size();
        int const n_points_2 = m_idx_range_perp_2.size();

        // Please provide at least 2 cells.
        if (n_points_1 < 2) {
            throw std::runtime_error("Please provide at least 2 cells in the Patch1.");
        }
        if (n_points_2 < 2) {
            throw std::runtime_error("Please provide at least 2 cells in the Patch2.");
        }

        // Initialise fields
        ddc::parallel_fill(weights_1_minus_2, 0.0);
        ddc::parallel_fill(weights_1_minus_1, 0.0);
        ddc::parallel_fill(weights_1_minus_0, 0.0);
        ddc::parallel_fill(weights_2_minus_2, 0.0);
        ddc::parallel_fill(weights_2_minus_1, 0.0);
        ddc::parallel_fill(weights_2_minus_0, 0.0);

        // Store the 3 last steps for the recursion in an array.
        std::array<host_t<DField<IdxRange1DPerp_1>>, 3> weights_patch1
                = {weights_1_minus_2, weights_1_minus_1, weights_1_minus_0};
        std::array<host_t<DField<IdxRange1DPerp_2>>, 3> weights_patch2
                = {weights_2_minus_2, weights_2_minus_1, weights_2_minus_0};

        std::array<double, 3> coeff_deriv_patch1;
        std::array<double, 3> coeff_deriv_patch2;

        // For the local coefficients
        double alpha_i;
        double beta_i;
        std::array<double, 3> gammas_i;

        // Define a step to increment the indices to move away from the interface.
        IdxStep<EdgePerpGrid1> idx_step_1
                = m_extremity_1 == FRONT ? IdxStep<EdgePerpGrid1>(1) : IdxStep<EdgePerpGrid1>(-1);
        IdxStep<EdgePerpGrid2> idx_step_2
                = m_extremity_2 == FRONT ? IdxStep<EdgePerpGrid2>(1) : IdxStep<EdgePerpGrid2>(-1);

        // Initialisation recursivity: set (a_{n-2},b_{n-2},c_{n-2}) and (a_{n-1},b_{n-1},c_{n-1}) -----------------------------
        Idx1D_1 interface_idx_1 = get_extremity_idx(m_extremity_1, m_idx_range_perp_1);
        Idx1D_2 interface_idx_2 = get_extremity_idx(m_extremity_2, m_idx_range_perp_2);

        // Set (a_{n-2},b_{n-2},c_{n-2}) ---
        Idx1D_1 idx_1_minus;
        Idx1D_1 idx_1 = interface_idx_1;
        Idx1D_1 idx_1_plus = interface_idx_1 + idx_step_1;

        Idx1D_2 idx_2_minus;
        Idx1D_2 idx_2 = interface_idx_2;
        Idx1D_2 idx_2_plus = interface_idx_2 + idx_step_2;

        auto [cell_length_left, cell_length_right] = get_cell_lengths(idx_2);

        alpha_i = get_alpha(cell_length_left, cell_length_right);
        beta_i = get_beta(cell_length_left, cell_length_right);
        gammas_i = get_gammas(cell_length_left, cell_length_right);

        coeff_deriv_patch1[0] = beta_i; // b_{n-2}
        coeff_deriv_patch2[0] = alpha_i; // a_{n-2}

        // c_{n-2}
        weights_patch1[0](idx_1_plus) = gammas_i[0];
        weights_patch1[0](idx_1) = gammas_i[1];
        weights_patch2[0](idx_2) = gammas_i[1];
        weights_patch2[0](idx_2_plus) = gammas_i[2];

        // Set (a_{n-1},b_{n-1},c_{n-1}) ---
        idx_2_minus = interface_idx_2;
        idx_2 = interface_idx_2 + idx_step_2;
        idx_2_plus = interface_idx_2 + 2 * idx_step_2;

        std::tie(cell_length_left, cell_length_right) = get_cell_lengths(idx_2);

        alpha_i = get_alpha(cell_length_left, cell_length_right);
        beta_i = get_beta(cell_length_left, cell_length_right);
        gammas_i = get_gammas(cell_length_left, cell_length_right);

        coeff_deriv_patch2[1]
                = alpha_i * coeff_deriv_patch2[0] / (1 - coeff_deriv_patch2[0] * beta_i); // a_{n-1}
        coeff_deriv_patch1[1]
                = coeff_deriv_patch1[0] / (1 - coeff_deriv_patch2[0] * beta_i); // b_{n-1}

        // c_{n-1} = (c_{n-2} + a_{n-2} * gamma_i) / (1 - a_{n-2} * beta_i);
        ddc::host_for_each(m_idx_range_perp_1, [&](Idx1D_1 const& idx) {
            weights_patch1[1](idx)
                    = 1. / (1 - coeff_deriv_patch2[0] * beta_i) * weights_patch1[0](idx);
        });
        ddc::host_for_each(m_idx_range_perp_2, [&](Idx1D_2 const& idx) {
            weights_patch2[1](idx)
                    = 1. / (1 - coeff_deriv_patch2[0] * beta_i) * weights_patch2[0](idx);
        });
        double factor = 1. / (1 - coeff_deriv_patch2[0] * beta_i) * coeff_deriv_patch2[0];
        weights_patch2[1](idx_2_plus) += factor * gammas_i[2];
        weights_patch2[1](idx_2) += factor * gammas_i[1];
        weights_patch2[1](idx_2_minus) += factor * gammas_i[0];
        weights_patch1[1](idx_1) += factor * gammas_i[0];

        // Set (a_{n},b_{n},c_{n}) ---
        coeff_deriv_patch1[2] = coeff_deriv_patch1[1];
        coeff_deriv_patch2[2] = coeff_deriv_patch2[1];
        ddc::parallel_deepcopy(weights_patch1[2], get_const_field(weights_patch1[1]));
        ddc::parallel_deepcopy(weights_patch2[2], get_const_field(weights_patch2[1]));


        // Computing coefficients (c, a, b) on patch 2 -------------------------------------------
        // If m_is_cell_bound_2_with_extra_interpol_pt is true, the last cell is treated differently.
        int const n2_cells
                = m_is_cell_bound_2_with_extra_interpol_pt ? n_points_2 - 3 : n_points_2 - 1;

        recursion(
                n2_cells,
                interface_idx_2,
                coeff_deriv_patch2,
                coeff_deriv_patch1,
                weights_patch1,
                weights_patch2);

        // Correction if we use an additional interpolation point as closure in the last cell.
        if (m_is_cell_bound_2_with_extra_interpol_pt) {
            correction_boundary(
                    n_points_2,
                    interface_idx_2,
                    coeff_deriv_patch2,
                    coeff_deriv_patch1,
                    weights_patch1,
                    weights_patch2);
        }

        // Initialisation second recursivity: set (a_{n-2},b_{n-2},c_{n-2}) and (a_{n-1},b_{n-1},c_{n-1}) ----------------------
        // Set (a_{n-2},b_{n-2},c_{n-2}) ---
        idx_1_minus = interface_idx_1;
        idx_1 = interface_idx_1 + idx_step_1;
        idx_1_plus = interface_idx_1 + 2 * idx_step_1;

        std::tie(cell_length_left, cell_length_right) = get_cell_lengths(idx_1);

        alpha_i = get_alpha(cell_length_left, cell_length_right);
        beta_i = get_beta(cell_length_left, cell_length_right);
        gammas_i = get_gammas(cell_length_left, cell_length_right);

        coeff_deriv_patch1[0] = coeff_deriv_patch1[2]; // a_{n-2}
        coeff_deriv_patch2[0] = coeff_deriv_patch2[2]; // b_{n-2}
        ddc::parallel_deepcopy(weights_patch1[0], get_const_field(weights_patch1[2])); // c_{n-2}
        ddc::parallel_deepcopy(weights_patch2[0], get_const_field(weights_patch2[2])); // c_{n-2}

        // Set (a_{n-1},b_{n-1},c_{n-1}) ---
        coeff_deriv_patch2[1]
                = coeff_deriv_patch2[0] / (1 - coeff_deriv_patch1[0] * alpha_i); // a_{n-1}
        coeff_deriv_patch1[1]
                = coeff_deriv_patch1[0] * beta_i / (1 - coeff_deriv_patch1[0] * alpha_i); // b_{n-1}

        // c_{n-1} = (c_{n-2} + b_{n-2} * gamma_i) / (1 - b_{n-2} * alpha_i);
        ddc::host_for_each(m_idx_range_perp_1, [&](Idx1D_1 const& idx) {
            weights_patch1[1](idx)
                    = 1. / (1 - coeff_deriv_patch1[0] * alpha_i) * weights_patch1[0](idx);
        });
        ddc::host_for_each(m_idx_range_perp_2, [&](Idx1D_2 const& idx) {
            weights_patch2[1](idx)
                    = 1. / (1 - coeff_deriv_patch1[0] * alpha_i) * weights_patch2[0](idx);
        });
        factor = 1. / (1 - coeff_deriv_patch1[0] * alpha_i) * coeff_deriv_patch1[0];
        weights_patch1[1](idx_1_plus) += factor * gammas_i[0];
        weights_patch1[1](idx_1) += factor * gammas_i[1];
        weights_patch1[1](idx_1_minus) += factor * gammas_i[2];
        weights_patch2[1](interface_idx_2) += factor * gammas_i[2];

        // Computing coefficients (c, a, b) on patch 1 -------------------------------------------
        // If m_is_cell_bound_1_with_extra_interpol_pt is true, the last cell is treated differently.
        int const n1_cells
                = m_is_cell_bound_1_with_extra_interpol_pt ? n_points_1 - 3 : n_points_1 - 1;

        recursion(
                n1_cells,
                interface_idx_1,
                coeff_deriv_patch1,
                coeff_deriv_patch2,
                weights_patch1,
                weights_patch2);

        // Correction if we use an additional interpolation point as closure in the last cell.
        if (m_is_cell_bound_1_with_extra_interpol_pt) {
            correction_boundary(
                    n_points_1,
                    interface_idx_1,
                    coeff_deriv_patch1,
                    coeff_deriv_patch2,
                    weights_patch1,
                    weights_patch2);
        }

        // Store the final values.
        m_coeff_deriv_patch_1 = coeff_deriv_patch1[2];
        m_coeff_deriv_patch_2 = coeff_deriv_patch2[2];
        ddc::parallel_deepcopy(m_weights_patch_1, get_const_field(weights_patch1[2]));
        ddc::parallel_deepcopy(m_weights_patch_2, get_const_field(weights_patch2[2]));
    }


    /**
     * @brief Recursion to compute the coefficients a, b and the weights omega. 
     * 
     * @param[in] n_cells Number of the cells of the patch where the recursion is made.
     * @param[in] interface_idx Index at the interface. 
     * 
     * The following parameters are array with three elements corresponding to the 
     * values at the step n-1, n and n+1 in the recursion. 
     * @param[in,out] coeff_deriv_patch Coefficient in front of the derivative of the patch. 
     *  E.g. Patch1: coeff_deriv_patch = b        |  Patch2: coeff_deriv_patch = a
     * @param[in,out] coeff_deriv_other_patch Coefficient in front of the derivative of the other patch. 
     *  E.g. Patch1: coeff_deriv_other_patch = a  |  Patch2: coeff_deriv_other_patch = b
     * @param[in,out] weights_patch1 Weights in front of the function values on patch 1. 
     * @param[in,out] weights_patch2 Weights in front of the function values on patch 2. 
     */
    template <class Grid1D>
    void recursion(
            int const n_cells,
            Idx<Grid1D> const interface_idx,
            std::array<double, 3>& coeff_deriv_patch,
            std::array<double, 3>& coeff_deriv_other_patch,
            std::array<host_t<DField<IdxRange1DPerp_1>>, 3>& weights_patch1,
            std::array<host_t<DField<IdxRange1DPerp_2>>, 3>& weights_patch2)
    {
        using Idx1D = Idx<Grid1D>;
        constexpr bool is_on_patch1 = std::is_same_v<Idx1D, Idx1D_1>;

        Extremity constexpr extremity = is_on_patch1 ? m_extremity_1 : m_extremity_2;
        IdxStep<Grid1D> idx_step = extremity == FRONT ? IdxStep<Grid1D>(1) : IdxStep<Grid1D>(-1);

        for (int n(2); n < n_cells; n++) {
            Idx1D const idx_minus = interface_idx + (n - 1) * idx_step;
            Idx1D const idx = interface_idx + n * idx_step;
            Idx1D const idx_plus = interface_idx + (n + 1) * idx_step;

            auto [cell_length_left, cell_length_right] = get_cell_lengths(idx);

            // Local coefficient (alpha_i/beta_i) furthest from the interface
            double const coeff_furthest = is_on_patch1
                                                  ? get_beta(cell_length_left, cell_length_right)
                                                  : get_alpha(cell_length_left, cell_length_right);
            // Local coefficient (alpha_i/beta_i) closest to the interface
            double const coeff_closest = is_on_patch1
                                                 ? get_alpha(cell_length_left, cell_length_right)
                                                 : get_beta(cell_length_left, cell_length_right);
            std::array<double, 3> const gammas_i = get_gammas(cell_length_left, cell_length_right);

            double const denom = 1 - coeff_closest * coeff_deriv_patch[1] / coeff_deriv_patch[0];

            // Set (a_{n},b_{n},c_{n}) ---
            coeff_deriv_patch[2] = coeff_furthest * coeff_deriv_patch[1] / denom;
            coeff_deriv_other_patch[2]
                    = (coeff_deriv_other_patch[1]
                       - coeff_closest * coeff_deriv_patch[1] / coeff_deriv_patch[0]
                                 * coeff_deriv_other_patch[0])
                      / denom;

            // Patch1: c_{n} = (c_{n-1} - alpha_i * b_{n-1} / b_{n-2} * c_{n-2} + b_{n-1} * gamma_i) / denom;
            // Patch2: c_{n} = (c_{n-1} - beta_i * a_{n-1} / a_{n-2} * c_{n-2} + a_{n-1} * gamma_i) / denom;
            double coeff_closer_deriv = coeff_closest * coeff_deriv_patch[1] / coeff_deriv_patch[0];
            ddc::host_for_each(m_idx_range_perp_1, [&](Idx1D_1 const& idx_perp) {
                weights_patch1[2](idx_perp)
                        = 1. / denom
                          * (weights_patch1[1](idx_perp)
                             - coeff_closer_deriv * weights_patch1[0](idx_perp));
            });
            ddc::host_for_each(m_idx_range_perp_2, [&](Idx1D_2 const& idx_perp) {
                weights_patch2[2](idx_perp)
                        = 1. / denom
                          * (weights_patch2[1](idx_perp)
                             - coeff_closer_deriv * weights_patch2[0](idx_perp));
            });
            double const factor = 1. / denom * coeff_deriv_patch[1];
            if constexpr (is_on_patch1) {
                weights_patch1[2](idx_plus) += factor * gammas_i[0];
                weights_patch1[2](idx) += factor * gammas_i[1];
                weights_patch1[2](idx_minus) += factor * gammas_i[2];
            } else {
                weights_patch2[2](idx_plus) += factor * gammas_i[2];
                weights_patch2[2](idx) += factor * gammas_i[1];
                weights_patch2[2](idx_minus) += factor * gammas_i[0];
            }

            // Update coefficients.
            coeff_deriv_patch[0] = coeff_deriv_patch[1];
            coeff_deriv_patch[1] = coeff_deriv_patch[2];
            coeff_deriv_other_patch[0] = coeff_deriv_other_patch[1];
            coeff_deriv_other_patch[1] = coeff_deriv_other_patch[2];

            ddc::parallel_deepcopy(weights_patch1[0], get_const_field(weights_patch1[1]));
            ddc::parallel_deepcopy(weights_patch2[0], get_const_field(weights_patch2[1]));

            ddc::parallel_deepcopy(weights_patch1[1], get_const_field(weights_patch1[2]));
            ddc::parallel_deepcopy(weights_patch2[1], get_const_field(weights_patch2[2]));
        }
    }

    /**
     * @brief Correction of the local coefficients on the boundary for additional interpolation 
     * points in the middle of a cell.  
     * 
     * @param[in] n_points Number of interpolation points on the patch where the correction is made.
     * @param[in] interface_idx Index at the interface. 
     * 
     * The following parameters are array with three elements corresponding to the 
     * values at the step n-1, n and n+1 in the recursion. 
     * @param[in,out] coeff_deriv_patch Coefficient in front of the derivative of the patch. 
     *  E.g. Patch1: coeff_deriv_patch = b        |  Patch2: coeff_deriv_patch = a
     * @param[in,out] coeff_deriv_other_patch Coefficient in front of the derivative of the other patch. 
     *  E.g. Patch1: coeff_deriv_other_patch = a  |  Patch2: coeff_deriv_other_patch = b
     * @param[in,out] weights_patch1 Weights in front of the function values on patch 1. 
     * @param[in,out] weights_patch2 Weights in front of the function values on patch 2. 
     */
    template <class Grid1D>
    void correction_boundary(
            int const n_points,
            Idx<Grid1D> const interface_idx,
            std::array<double, 3>& coeff_deriv_patch,
            std::array<double, 3>& coeff_deriv_other_patch,
            std::array<host_t<DField<IdxRange1DPerp_1>>, 3>& weights_patch1,
            std::array<host_t<DField<IdxRange1DPerp_2>>, 3>& weights_patch2)
    {
        using Idx1D = Idx<Grid1D>;
        constexpr bool is_on_patch1 = std::is_same_v<Idx1D, Idx1D_1>;

        Extremity constexpr extremity = is_on_patch1 ? m_extremity_1 : m_extremity_2;
        IdxStep<Grid1D> idx_step = extremity == FRONT ? IdxStep<Grid1D>(1) : IdxStep<Grid1D>(-1);

        int const n = n_points - 3;
        Idx1D const idx_minus = interface_idx + (n - 1) * idx_step;
        Idx1D const idx = interface_idx + n * idx_step;
        Idx1D const idx_plus = interface_idx + (n + 2) * idx_step;
        // Additional interpolation point.
        Idx1D const idx_star = interface_idx + (n + 1) * idx_step;

        auto [cell_length_left, cell_length_right] = get_cell_lengths(idx);
        auto [cell_length_left_star, cell_length_right_star] = get_cell_lengths(idx_star);

        if constexpr (is_on_patch1) {
            cell_length_left += cell_length_left_star;
        } else {
            cell_length_right += cell_length_right_star;
        }

        // Local coefficient (alpha_i/beta_i) furthest from the interface
        double coeff_furthest = is_on_patch1 ? get_beta(cell_length_left, cell_length_right)
                                             : get_alpha(cell_length_left, cell_length_right);
        // Local coefficient (alpha_i/beta_i) closest to the interface
        double coeff_closest = is_on_patch1 ? get_alpha(cell_length_left, cell_length_right)
                                            : get_beta(cell_length_left, cell_length_right);
        std::array<double, 3> const gammas_i = get_gammas(cell_length_left, cell_length_right);

        // Correction of gammas_i, alpha_i and beta_i
        double x;
        if constexpr (is_on_patch1) {
            // (x^-_star - x^-_{i-1}) / Delta x^-_i
            x = cell_length_left_star / cell_length_left;
        } else {
            // (x^+_star - x^+_i) / Delta x^+_i
            x = cell_length_left_star / cell_length_right;
        }
        double const H0 = (1 - x) * (1 - x) * (1 + 2 * x);
        double const H1 = x * x * (3 - 2 * x);
        double const K0 = (1 - x) * (1 - x) * x;
        double const K1 = x * x * (x - 1);

        double gamma_plus;
        double gamma_star;
        double gamma_idx;
        double gamma_minus;
        if constexpr (is_on_patch1) {
            double denom_gamma = 1. / (1 + coeff_furthest * K1 / K0);

            coeff_closest = denom_gamma * coeff_closest;

            gamma_plus = denom_gamma * (gammas_i[0] - coeff_furthest * H0 / K0 / cell_length_left);
            gamma_star = denom_gamma * coeff_furthest / K0 / cell_length_left;
            gamma_idx = denom_gamma * (gammas_i[1] - coeff_furthest * H1 / K0 / cell_length_left);
            gamma_minus = denom_gamma * gammas_i[2];
        } else {
            double denom_gamma = 1. / (1 + coeff_furthest * K0 / K1);

            coeff_closest = denom_gamma * coeff_closest;

            gamma_plus = denom_gamma * (gammas_i[2] - coeff_furthest * H1 / K1 / cell_length_right);
            gamma_star = denom_gamma * coeff_furthest / K1 / cell_length_right;
            gamma_idx = denom_gamma * (gammas_i[1] - coeff_furthest * H0 / K1 / cell_length_right);
            gamma_minus = denom_gamma * gammas_i[0];
        }

        // No more boundary derivative, the local coefficient is set to zero.
        coeff_furthest = 0;

        double denom = 1 - coeff_closest * coeff_deriv_patch[1] / coeff_deriv_patch[0];

        // Set (a_{n},b_{n},c_{n}) ---
        coeff_deriv_patch[2] = coeff_furthest * coeff_deriv_patch[1] / denom;
        coeff_deriv_other_patch[2] = (coeff_deriv_other_patch[1]
                                      - coeff_closest * coeff_deriv_patch[1] / coeff_deriv_patch[0]
                                                * coeff_deriv_other_patch[0])
                                     / denom;

        // Patch1: c_{n} = (c_{n-1} - alpha_i * b_{n-1} / b_{n-2} * c_{n-2} + b_{n-1} * gamma_i) / denom;
        // Patch2: c_{n} = (c_{n-1} - beta_i * a_{n-1} / a_{n-2} * c_{n-2} + a_{n-1} * gamma_i) / denom;
        double coeff_closer_deriv = coeff_closest * coeff_deriv_patch[1] / coeff_deriv_patch[0];
        ddc::host_for_each(m_idx_range_perp_1, [&](Idx1D_1 const& idx_perp) {
            weights_patch1[2](idx_perp) = 1. / denom
                                          * (weights_patch1[1](idx_perp)
                                             - coeff_closer_deriv * weights_patch1[0](idx_perp));
        });
        ddc::host_for_each(m_idx_range_perp_2, [&](Idx1D_2 const& idx_perp) {
            weights_patch2[2](idx_perp) = 1. / denom
                                          * (weights_patch2[1](idx_perp)
                                             - coeff_closer_deriv * weights_patch2[0](idx_perp));
        });
        double const factor = 1. / denom * coeff_deriv_patch[1];
        if constexpr (is_on_patch1) {
            weights_patch1[2](idx_plus) += factor * gamma_plus;
            weights_patch1[2](idx_star) += factor * gamma_star;
            weights_patch1[2](idx) += factor * gamma_idx;
            weights_patch1[2](idx_minus) += factor * gamma_minus;
        } else {
            weights_patch2[2](idx_plus) += factor * gamma_plus;
            weights_patch2[2](idx_star) += factor * gamma_star;
            weights_patch2[2](idx) += factor * gamma_idx;
            weights_patch2[2](idx_minus) += factor * gamma_minus;
        }
    }


    // ===========================================================================================
    // UNIFORM INTERPOLATION POINTS                                                              |
    // ===========================================================================================
    /**
     * @brief Compute the coefficients a, b and the weights omega applying the explicit
     * formula for the uniform case. 
     */
    void set_coefficients_uniform_per_patch_case()
    {
        int const n_cells_1 = m_idx_range_perp_1.size() - 1;
        int const n_cells_2 = m_idx_range_perp_2.size() - 1;

        Idx1D_1 interface_idx_1 = get_extremity_idx(m_extremity_1, m_idx_range_perp_1);
        Idx1D_2 interface_idx_2 = get_extremity_idx(m_extremity_2, m_idx_range_perp_2);

        Idx1D_1 idx_1 = interface_idx_1;
        Idx1D_2 idx_2 = interface_idx_2;

        // Index increment to move away from the interface.
        IdxStep<EdgePerpGrid1> idx_step_1
                = m_extremity_1 == FRONT ? IdxStep<EdgePerpGrid1>(1) : IdxStep<EdgePerpGrid1>(-1);
        IdxStep<EdgePerpGrid2> idx_step_2
                = m_extremity_2 == FRONT ? IdxStep<EdgePerpGrid2>(1) : IdxStep<EdgePerpGrid2>(-1);

        auto [cell_length_left, cell_length_right] = get_cell_lengths(idx_1);

        double const a_11 = -0.5 * cell_length_left / (cell_length_left + cell_length_right);
        double const b_11 = -0.5 * cell_length_right / (cell_length_left + cell_length_right);

        double const u1 = 2 * Kokkos::sqrt(3);

        double const un1 = Kokkos::pow(2 + Kokkos::sqrt(3), n_cells_1)
                           - Kokkos::pow(2 - Kokkos::sqrt(3), n_cells_1);
        double const un1_minus = Kokkos::pow(2 + Kokkos::sqrt(3), n_cells_1 - 1)
                                 - Kokkos::pow(2 - Kokkos::sqrt(3), n_cells_1 - 1);

        double const un2 = Kokkos::pow(2 + Kokkos::sqrt(3), n_cells_2)
                           - Kokkos::pow(2 - Kokkos::sqrt(3), n_cells_2);
        double const un2_minus = Kokkos::pow(2 + Kokkos::sqrt(3), n_cells_2 - 1)
                                 - Kokkos::pow(2 - Kokkos::sqrt(3), n_cells_2 - 1);

        double const denominator = un1 * un2 + un1 * un2_minus * a_11 + un1_minus * un2 * b_11;


        // Compute the coefficients a and b.
        m_coeff_deriv_patch_1
                = Kokkos::pow(-1, (n_cells_1 - 1) % 2) * u1 * b_11 * un2 / denominator;
        m_coeff_deriv_patch_2
                = Kokkos::pow(-1, (n_cells_2 - 1) % 2) * u1 * a_11 * un1 / denominator;


        // Compute the weights {omega_k}_{k = -NL, ..., NR} in c.
        double const factor_a = 3 * a_11 / cell_length_right / denominator;
        double const factor_b = 3 * b_11 / cell_length_left / denominator;

        // --- for k = 0
        m_weights_patch_1(idx_1)
                = factor_a * un1 * (un2 - un2_minus) - factor_b * un2 * (un1 - un1_minus);
        m_weights_patch_2(idx_2) = m_weights_patch_1(idx_1);

        // --- for k = 1, ..., NR-1
        for (int k(1); k < n_cells_2; ++k) {
            idx_2 = idx_2 + idx_step_2;
            const int n2_minus_k = n_cells_2 - k;
            double const vk_minus = Kokkos::pow(2 + Kokkos::sqrt(3), n2_minus_k - 1)
                                    - Kokkos::pow(2 - Kokkos::sqrt(3), n2_minus_k - 1);
            double const vk_plus = Kokkos::pow(2 + Kokkos::sqrt(3), n2_minus_k + 1)
                                   - Kokkos::pow(2 - Kokkos::sqrt(3), n2_minus_k + 1);

            m_weights_patch_2(idx_2)
                    = Kokkos::pow(-1, k % 2) * factor_a * un1 * (vk_plus - vk_minus);
        }

        // --- for k = NR
        idx_2 = idx_2 + idx_step_2;
        m_weights_patch_2(idx_2) = Kokkos::pow(-1, n_cells_2 % 2) * factor_a * un1 * u1;

        // --- for k = -1, ..., -(NL-1)
        for (int k(1); k < n_cells_1; ++k) {
            idx_1 = idx_1 + idx_step_1;
            const int n1_minus_k = n_cells_1 - k;
            double const vk_minus = Kokkos::pow(2 + Kokkos::sqrt(3), n1_minus_k - 1)
                                    - Kokkos::pow(2 - Kokkos::sqrt(3), n1_minus_k - 1);
            double const vk_plus = Kokkos::pow(2 + Kokkos::sqrt(3), n1_minus_k + 1)
                                   - Kokkos::pow(2 - Kokkos::sqrt(3), n1_minus_k + 1);

            m_weights_patch_1(idx_1)
                    = Kokkos::pow(-1, (k + 1) % 2) * factor_b * un2 * (vk_plus - vk_minus);
        }

        // --- for k = -NL
        idx_1 = idx_1 + idx_step_1;
        m_weights_patch_1(idx_1) = Kokkos::pow(-1, (n_cells_1 + 1) % 2) * factor_b * un2 * u1;
    };


    /**
     * @brief Get the lengths of the cell on the right side and the left side of the given index.
     */
    template <class Grid1D>
    std::tuple<double, double> get_cell_lengths(Idx<Grid1D> const& idx) const
    {
        static_assert(
                std::is_same_v<Grid1D, EdgePerpGrid1> || std::is_same_v<Grid1D, EdgePerpGrid2>,
                "Wrong type of index given.");

        bool constexpr is_on_patch1 = std::is_same_v<Grid1D, EdgePerpGrid1>;

        // Tag for the grid on the other patch.
        using OGrid1D = std::conditional_t<is_on_patch1, EdgePerpGrid2, EdgePerpGrid1>;

        using IdxRangeGlobal = IdxRange<EdgePerpGrid1, EdgePerpGrid2>;

        Extremity constexpr extremity = is_on_patch1 ? m_extremity_1 : m_extremity_2;
        Extremity constexpr other_extremity = is_on_patch1 ? m_extremity_2 : m_extremity_1;

        IdxRange<Grid1D> idx_range_1d(IdxRangeGlobal(m_idx_range_perp_1, m_idx_range_perp_2));
        IdxRange<OGrid1D> other_idx_range_1d(
                IdxRangeGlobal(m_idx_range_perp_1, m_idx_range_perp_2));

        // Please, do not provide a index on the boundary except for the interface.
        assert(!(idx == idx_range_1d.front() && extremity == BACK));
        assert(!(idx == idx_range_1d.back() && extremity == FRONT));

        double length_left;
        double length_right;

        Idx<Grid1D> const idx_minus = idx - IdxStep<Grid1D>(1);
        Idx<Grid1D> const idx_plus = idx + IdxStep<Grid1D>(1);

        double const delta_coord_plus = abs(ddc::coordinate(idx_plus) - ddc::coordinate(idx));
        double const delta_coord_minus = abs(ddc::coordinate(idx) - ddc::coordinate(idx_minus));

        bool constexpr is_same_orientation_as_global = (is_on_patch1 && m_extremity_1 == BACK)
                                                       || (!is_on_patch1 && m_extremity_2 == FRONT);

        // If the given index corresponds to the index at the interface,
        Idx<Grid1D> const idx_interface = get_extremity_idx(extremity, idx_range_1d);
        if (idx == idx_interface) {
            IdxStep<Grid1D> idx_step
                    = extremity == FRONT ? IdxStep<Grid1D>(1) : IdxStep<Grid1D>(-1);
            Idx<Grid1D> const idx_incremented = idx + idx_step;

            IdxStep<OGrid1D> other_idx_step
                    = other_extremity == FRONT ? IdxStep<OGrid1D>(1) : IdxStep<OGrid1D>(-1);
            Idx<OGrid1D> const other_idx = get_extremity_idx(other_extremity, other_idx_range_1d);
            Idx<OGrid1D> const other_idx_incremented = other_idx + other_idx_step;

            length_left = abs(ddc::coordinate(other_idx_incremented) - ddc::coordinate(other_idx));
            length_right = abs(ddc::coordinate(idx_incremented) - ddc::coordinate(idx));
        }
        // If given index is on the patch 1 or on the patch 2,
        else if (is_same_orientation_as_global) {
            // . | →   or  → | .
            length_left = delta_coord_minus;
            length_right = delta_coord_plus;
        } else {
            //  . | ←  or  ← | .
            length_left = delta_coord_plus;
            length_right = delta_coord_minus;
        }

        if ((length_left == 0) || (length_right == 0)) {
            throw std::runtime_error(
                    "[abort] Ill-defined: the length of the cells must be not zero.");
        }
        return std::make_tuple(length_left, length_right);
    }


    /// @brief Compute the coefficients gammas_i.
    std::array<double, 3> get_gammas(double const cell_length_left, double const cell_length_right)
            const
    {
        double const factor = 3. / 2. / (cell_length_left + cell_length_right);
        double const R_on_L = cell_length_right / cell_length_left;
        double const L_on_R = cell_length_left / cell_length_right;

        double gamma_i_minus = -factor * R_on_L;
        double gamma_i = factor * (R_on_L - L_on_R);
        double gamma_i_plus = factor * L_on_R;

        return std::array<double, 3>({gamma_i_minus, gamma_i, gamma_i_plus});
    }

    /// @brief Compute the coefficient alpha_i.
    double get_alpha(double const cell_length_left, double const cell_length_right) const
    {
        return -0.5 * cell_length_left / (cell_length_right + cell_length_left);
    }

    /// @brief Compute the coefficient beta_i.
    double get_beta(double const cell_length_left, double const cell_length_right) const
    {
        return -0.5 * cell_length_right / (cell_length_right + cell_length_left);
    }


    template <class Grid1D>
    KOKKOS_INLINE_FUNCTION Idx<Grid1D> get_extremity_idx(
            Extremity const extremity,
            IdxRange<Grid1D> const& idx_range) const
    {
        return (extremity == FRONT) ? idx_range.front() : idx_range.back();
    }
};



template <class InterfaceType>
inline constexpr bool enable_single_derivative_calculator<
        SingleInterfaceDerivativesCalculator<InterfaceType>> = true;
