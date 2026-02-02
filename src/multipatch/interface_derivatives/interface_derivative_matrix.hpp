// SPDX-License-Identifier: MIT

#pragma once

#include <typeinfo>

#include <ddc/ddc.hpp>

#include "connectivity_details.hpp"
#include "ddc_aliases.hpp"
#include "edge.hpp"
#include "edge_transformation.hpp"
#include "interface.hpp"
#include "multipatch_field.hpp"
#include "multipatch_type.hpp"
#include "single_interface_derivatives_calculator_collection.hpp"
#include "types.hpp"
#include "view.hpp"

template <class Connectivity, class Grid1D, class PatchSeq, class DerivativesCalculatorCollection>
class InterfaceDerivativeMatrix;

/**
  * @brief Class to compute the interface derivatives along a given direction 
  * on all the interfaces with an approximation formula.
  * 
  * This operator is implemented for conforming equivalent global meshes. 
  * 
  * When we call the operator .solve_deriv() or .solve_cross_deriv(), the operator loops
  * over the interface to compute the vector C_trunc to determine all the interface derivatives 
  * and update the DerivField given in input. 
  * 
  * The vector C_trunc is given by the coefficients (see README),
  * @f$ \{\sum_{k = - N_{reduc}}^{N_{reduc}} \omega_{k, N_{reduc}, N_{reduc}}^{i_I} f_{i_I+k} \}_I @f$
  * 
  * with 
  *     * @f$ N_{reduc} @f$ the number of chosen cells for the approximation, 
  *     * @f$ i_I @f$ the index of the grid at the Ith interface, 
  *     * @f$ f_{i_I+k} @f$ the function values if we compute the first derivatives, or the first 
  * derivatives if we compute the cross-derivatives. 
  * 
  * @tparam Connectivity A MultipatchConnectivity class describing all the patch connections.
  * @tparam Grid1D A given direction.
  * @tparam Patches List of patches containing all the involved patches in the given direction. 
  * @tparam DerivativesCalculatorCollection A SingleInterfaceDerivativesCalculatorCollection
  * that stores the SingleInterfaceDerivativesCalculator needed to compute the interface derivatives. 
  */
template <class Connectivity, class Grid1D, class... Patches, class DerivativesCalculatorCollection>
class InterfaceDerivativeMatrix<
        Connectivity,
        Grid1D,
        ddc::detail::TypeSeq<Patches...>,
        DerivativesCalculatorCollection>
{
    /*
        All the interfaces given as input to the MultipatchConnectivity class.
         We expect the parameters defined on these interfaces.
    */
    using all_interface_collection = typename Connectivity::interface_collection;
    // All the sorted interfaces with the correct orientation.
    using interface_sorted_collection =
            typename Connectivity::template get_all_interfaces_along_direction_t<Grid1D>;

    // All the patches of the geometry.
    using all_patches = typename Connectivity::all_patches;

    static_assert(
            is_single_derivative_calculator_collection_v<DerivativesCalculatorCollection>,
            "Please provide a SingleInterfaceDerivativesCalculatorCollection type.");

    static constexpr std::size_t number_of_interfaces
            = ddc::type_seq_size_v<interface_sorted_collection>;

    using Interface0 = ddc::type_seq_element_t<0, interface_sorted_collection>;
    using InterfaceN
            = ddc::type_seq_element_t<number_of_interfaces - 1, interface_sorted_collection>;

    static_assert(
            (((!std::is_same_v<typename Interface0::Edge1, OutsideEdge>)
              || (!std::is_same_v<typename Interface0::Edge2, OutsideEdge>))
             == ((!std::is_same_v<typename InterfaceN::Edge1, OutsideEdge>)
                 || (!std::is_same_v<typename InterfaceN::Edge2, OutsideEdge>))),
            "If one global boundary has an outside edge, the other global boundary should have "
            "too.");

    static constexpr bool is_periodic
            = (!std::is_same_v<typename Interface0::Edge1, OutsideEdge>)&&(
                    !std::is_same_v<typename Interface0::Edge2, OutsideEdge>);
    /*
        Remove all the interfaces with an OutsideEdge.
        Rely on get_all_interfaces_along_direction_t to order the interfaces.
    */
    using outer_interface_collection = std::conditional_t<
            is_periodic,
            ddc::detail::TypeSeq<>,
            ddc::detail::TypeSeq<Interface0, InterfaceN>>;

    // All the interfaces in the Grid1D direction without the interfaces with OutsideEdge.
    using inner_interface_collection
            = ddc::type_seq_remove_t<interface_sorted_collection, outer_interface_collection>;

    // Number of inner interfaces. It fixes the number of elements in the matrix and vectors.
    static constexpr std::size_t n_inner_interfaces
            = ddc::type_seq_size_v<inner_interface_collection>;

    using FirstInterface = ddc::type_seq_element_t<0, inner_interface_collection>;

    // First patch of the sorted interfaces.
    using FirstPatch = std::conditional_t<
            is_periodic,
            find_patch_t<Grid1D, all_patches>,
            typename FirstInterface::Edge2::associated_patch>;

    // Sequence of 1D grids in the direction of Grid1D over the patches.
    using Grid1DSeq = collect_grids_on_dim_t<
            find_patch_t<Grid1D, all_patches>,
            Grid1D,
            interface_sorted_collection>;

    static_assert(
            ddc::type_seq_size_v<Grid1DSeq> == n_inner_interfaces + !is_periodic,
            "The number of 1D grids and inner interfaces should match.");

    static_assert(
            ddc::type_seq_size_v<Grid1DSeq> <= sizeof...(Patches),
            "The number of given patches is insufficient.");

    static_assert(
            ddc::type_seq_contains_v<
                    Grid1DSeq,
                    ddc::type_seq_merge_t<
                            ddc::detail::TypeSeq<typename Patches::Grid1...>,
                            ddc::detail::TypeSeq<typename Patches::Grid2...>>>,
            "A patch is missing in the given patch list. The list has to contain all the patches "
            "in the given Grid1D direction.");

    // Tag to identify we are computing the first derivatives.
    struct eval_deriv
    {
    };

    // Tag to identify we are computing the cross-derivatives.
    struct eval_cross_deriv
    {
    };

public:
    /// @brief First patch of the collection of patches in the direction of the given Grid1D.
    using first_patch = FirstPatch;

private:
    MultipatchType<IdxRangeOnPatch, Patches...> const& m_idx_ranges;

    DerivativesCalculatorCollection const& m_derivatives_calculators;

public:
    /**
     * @brief Instantiate InterfaceDerivativeMatrix. 
     *  
     * @param idx_ranges MultipatchType collection of index ranges defined on the given list of patches. 
     * @param derivatives_calculators SingleInterfaceDerivativesCalculatorCollection containing all the 
     *          interface derivative calculator for each interface in the given Grid1D direction. 
     */
    InterfaceDerivativeMatrix(
            MultipatchType<IdxRangeOnPatch, Patches...> const& idx_ranges,
            DerivativesCalculatorCollection const& derivatives_calculators)
        : m_idx_ranges(idx_ranges)
        , m_derivatives_calculators(derivatives_calculators)
    {
    }

    /**
     * @brief Compute the vector C_trunc to determine all the interface derivatives in the Grid1D
     * direction at a given index.
     * 
     * It uses the function values to compute the first-derivatives perpendicular to the interfaces. 
     *
     * @tparam IdxPar Type for the given index. It can be on the first or the second grid of the first patch.
     * @param[in,out] functions_and_derivs Collection fields with its first derivatives and cross-derivatives.
     * @param[in] idx_par Index of the line in the Grid1D direction where we compute all the interface derivatives.
     */
    template <class IdxPar>
    void solve_deriv(
            MultipatchField<DerivFieldOnPatch_host, Patches...> functions_and_derivs,
            IdxPar const& idx_par)
    {
        solve<eval_deriv>(functions_and_derivs, idx_par);
    }

    /**
    * @brief Compute the vector C_trunc to determine all the interface derivatives in the Grid1D
     * direction at all the indices of the first patch.
     * 
     * It uses the function values to compute the first-derivatives perpendicular to the interfaces. 
     * 
     * @param[in,out] functions_and_derivs Collection fields with its first derivatives and cross-derivatives 
     *              we want to compute.
     * @warning The global mesh is supposed to be conforming. 
     */
    void solve_deriv(MultipatchField<DerivFieldOnPatch_host, Patches...> functions_and_derivs)
    {
        constexpr bool is_first_idx_range_perp_on_dim1
                = ddc::in_tags_v<typename FirstPatch::Grid1, Grid1DSeq>;
        using IdxRangeParFirstType = std::conditional_t<
                is_first_idx_range_perp_on_dim1,
                typename FirstPatch::IdxRange2,
                typename FirstPatch::IdxRange1>;

        IdxRangeParFirstType idx_range_par_first(m_idx_ranges.template get<FirstPatch>());

        ddc::host_for_each(
                idx_range_par_first,
                [&](typename IdxRangeParFirstType::discrete_element_type const& idx_par) {
                    solve<eval_deriv>(functions_and_derivs, idx_par);
                });
    }

    /**
     * @brief Compute the vector C_trunc to determine all the interface cross-derivatives in the Grid1D
     * direction at a given index.
     * 
     * It uses the first-derivatives s to compute the cross-derivatives in a perpendicular direction to the interfaces. 
     *
     * @tparam IdxPar Type for the given index. It can be on the first or the second grid of the first patch.
     * @param[in,out] functions_and_derivs Collection fields with its first derivatives and cross-derivatives.
     * @param[in] idx_par Index of the line in the Grid1D direction where we compute all the interface cross-derivatives.
     */
    template <class IdxPar>
    void solve_cross_deriv(
            MultipatchField<DerivFieldOnPatch_host, Patches...> functions_and_derivs,
            IdxPar const& idx_par)
    {
        solve<eval_cross_deriv>(functions_and_derivs, idx_par);
    }

    /**
     * @brief Compute the vector C_trunc to determine all the corner cross-derivatives in the Grid1D
     * direction at all the last and the first indices of the first patch.
     * 
     * It uses the first-derivatives s to compute the cross-derivatives in a perpendicular direction to the interfaces. 
     *
     * @param[in,out] functions_and_derivs Collection fields with its first derivatives and cross-derivatives.
     * @warning The global mesh is supposed to be conforming. 
     */
    void solve_cross_deriv(MultipatchField<DerivFieldOnPatch_host, Patches...> functions_and_derivs)
    {
        constexpr bool is_first_idx_range_perp_on_dim1
                = ddc::in_tags_v<typename FirstPatch::Grid1, Grid1DSeq>;
        using IdxRangeParFirstType = std::conditional_t<
                is_first_idx_range_perp_on_dim1,
                typename FirstPatch::IdxRange2,
                typename FirstPatch::IdxRange1>;

        IdxRangeParFirstType idx_range_par_first(m_idx_ranges.template get<FirstPatch>());

        solve<eval_cross_deriv>(functions_and_derivs, idx_range_par_first.front());
        solve<eval_cross_deriv>(functions_and_derivs, idx_range_par_first.back());
    }


private:
    template <typename eval_type, class IdxPar>
    void solve(
            MultipatchField<DerivFieldOnPatch_host, Patches...> functions_and_derivs,
            IdxPar const& idx_par)
    {
        // Check the input types.
        static_assert(
                (std::is_same_v<IdxPar, typename FirstPatch::Idx1>)
                        || (std::is_same_v<IdxPar, typename FirstPatch::Idx2>),
                "The given index has to be a 1D index defined on the first patch. See "
                "InterfaceDerivativeMatrix<...>::first_patch.");

        // Update derivatives.
        update_derivatives<eval_type>(
                functions_and_derivs,
                idx_par,
                std::make_integer_sequence<std::size_t, n_inner_interfaces> {});
    }


    /// @brief Once the interface derivatives computed, fill in the correct derivative fields.
    template <typename eval_type, class GridPar1D, std::size_t... I>
    void update_derivatives(
            MultipatchField<DerivFieldOnPatch_host, Patches...> functions_and_derivs,
            Idx<GridPar1D> const& slice_idx,
            std::integer_sequence<std::size_t, I...>)
    {
        static_assert(
                (std::is_same_v<GridPar1D, typename FirstPatch::Grid1>)
                || (std::is_same_v<GridPar1D, typename FirstPatch::Grid2>)
                || (std::is_same_v<GridPar1D, ddc::Deriv<typename FirstPatch::Dim1>>)
                || (std::is_same_v<GridPar1D, ddc::Deriv<typename FirstPatch::Dim2>>));

        // Locally transform the index into an integer to avoid type issues.
        IdxRange<GridPar1D> idx_range_1d_first_patch(m_idx_ranges.template get<FirstPatch>());
        int slice_idx_value = (slice_idx - idx_range_1d_first_patch.front()).value();
        (update_derivatives_at_interface<I, eval_type>(functions_and_derivs, slice_idx_value), ...);
    }


    /// @brief Associate the Ith derivative values to the correct first derivative.
    template <
            std::size_t I,
            typename eval_type,
            std::enable_if_t<std::is_same_v<eval_type, eval_deriv>, bool> = true>
    void update_derivatives_at_interface(
            MultipatchField<DerivFieldOnPatch_host, Patches...> functions_and_derivs,
            int& slice_previous_idx_value)
    {
        using InterfaceI = ddc::type_seq_element_t<I, inner_interface_collection>;
        // Equivalent interface defined in the main file, i.e. given to MultipatchConnectivity.
        using EquivalentInterfaceI
                = find_associated_interface_t<typename InterfaceI::Edge1, all_interface_collection>;

        // The orientation Patch1|Patch2 of the interface matches with the sorted 1D grid sequence.
        constexpr bool is_good_orientation = std::is_same_v<
                typename EquivalentInterfaceI::Edge1::perpendicular_grid,
                ddc::type_seq_element_t<I, Grid1DSeq>>;

        using Patch_1 = typename EquivalentInterfaceI::Edge1::associated_patch;
        using Patch_2 = typename EquivalentInterfaceI::Edge2::associated_patch;

        constexpr Extremity extremity_1 = EquivalentInterfaceI::Edge1::extremity;
        constexpr Extremity extremity_2 = EquivalentInterfaceI::Edge2::extremity;

        using GridPerp1 = typename EquivalentInterfaceI::Edge1::perpendicular_grid;
        using GridPerp2 = typename EquivalentInterfaceI::Edge2::perpendicular_grid;

        using GridPar1 = typename EquivalentInterfaceI::Edge1::parallel_grid;
        using GridPar2 = typename EquivalentInterfaceI::Edge2::parallel_grid;

        using DerivPerp1 = typename ddc::Deriv<typename GridPerp1::continuous_dimension_type>;
        using DerivPerp2 = typename ddc::Deriv<typename GridPerp2::continuous_dimension_type>;

        // Type for the slice_previous_idx_value index.
        using IdxPatchSlice = std::conditional_t<is_good_orientation, Patch_1, Patch_2>;
        Idx<GridPar1> idx_slice_1;
        Idx<GridPar2> idx_slice_2;
        set_slice_indexes<
                EquivalentInterfaceI,
                IdxPatchSlice>(slice_previous_idx_value, idx_slice_1, idx_slice_2);

        // Get the fields of the left and right patch of the interface.
        DerivFieldOnPatch_host<Patch_1> function_and_derivs_1
                = functions_and_derivs.template get<Patch_1>();
        DerivFieldOnPatch_host<Patch_2> function_and_derivs_2
                = functions_and_derivs.template get<Patch_2>();

        // Use the function values to compute the first derivatives.
        DField<typename Patch_1::IdxRange12, Kokkos::HostSpace, Kokkos::layout_stride> function_1
                = function_and_derivs_1.get_values_field();
        DField<typename Patch_2::IdxRange12, Kokkos::HostSpace, Kokkos::layout_stride> function_2
                = function_and_derivs_2.get_values_field();

        // Compute the coefficient c_I for the interface I.
        double const interface_deriv
                = m_derivatives_calculators.template get<EquivalentInterfaceI>()
                          .get_function_coefficients(
                                  get_const_field(function_1[idx_slice_1]),
                                  get_const_field(function_2[idx_slice_2]));

        // Get the correct index ranges and indices for the slices.
        IdxRange<GridPerp1> idx_range_perp_1(m_idx_ranges.template get<Patch_1>());
        IdxRange<GridPerp2> idx_range_perp_2(m_idx_ranges.template get<Patch_2>());

        Idx<GridPerp1> idx_deriv_1 = get_idx_interface(idx_range_perp_1, extremity_1);
        Idx<GridPerp2> idx_deriv_2 = get_idx_interface(idx_range_perp_2, extremity_2);

        Idx<DerivPerp1, GridPerp1> idx_slice_deriv_1(Idx<DerivPerp1>(1), idx_deriv_1);
        Idx<DerivPerp2, GridPerp2> idx_slice_deriv_2(Idx<DerivPerp2>(1), idx_deriv_2);

        constexpr bool change_sign_1 = (extremity_1 == Extremity::FRONT);
        constexpr bool change_sign_2 = (extremity_2 == Extremity::BACK);
        constexpr int sign_1 = !change_sign_1 - change_sign_1;
        constexpr int sign_2 = !change_sign_2 - change_sign_2;

        // Update the derivatives.
        function_and_derivs_1(idx_slice_deriv_1, idx_slice_1) = interface_deriv * sign_1;
        function_and_derivs_2(idx_slice_deriv_2, idx_slice_2) = interface_deriv * sign_2;
    }

    /// @brief Associate the Ith derivative values to the correct cross-derivative.
    template <
            std::size_t I,
            typename eval_type,
            std::enable_if_t<std::is_same_v<eval_type, eval_cross_deriv>, bool> = true>
    void update_derivatives_at_interface(
            MultipatchField<DerivFieldOnPatch_host, Patches...> functions_and_derivs,
            int& slice_previous_idx_value)
    {
        using InterfaceI = ddc::type_seq_element_t<I, inner_interface_collection>;
        // Equivalent interface defined in the main file, i.e. given to MultipatchConnectivity.
        using EquivalentInterfaceI
                = find_associated_interface_t<typename InterfaceI::Edge1, all_interface_collection>;

        // The orientation Patch1|Patch2 of the interface matches with the sorted 1D grid sequence.
        constexpr bool is_good_orientation = std::is_same_v<
                typename EquivalentInterfaceI::Edge1::perpendicular_grid,
                ddc::type_seq_element_t<I, Grid1DSeq>>;

        using Patch_1 = typename EquivalentInterfaceI::Edge1::associated_patch;
        using Patch_2 = typename EquivalentInterfaceI::Edge2::associated_patch;

        constexpr Extremity extremity_1 = EquivalentInterfaceI::Edge1::extremity;
        constexpr Extremity extremity_2 = EquivalentInterfaceI::Edge2::extremity;

        using GridPerp1 = typename EquivalentInterfaceI::Edge1::perpendicular_grid;
        using GridPerp2 = typename EquivalentInterfaceI::Edge2::perpendicular_grid;

        using GridPar1 = typename EquivalentInterfaceI::Edge1::parallel_grid;
        using GridPar2 = typename EquivalentInterfaceI::Edge2::parallel_grid;

        using DerivPar1 = ddc::Deriv<typename GridPar1::continuous_dimension_type>;
        using DerivPar2 = ddc::Deriv<typename GridPar2::continuous_dimension_type>;

        using Grid1_1 = typename Patch_1::Grid1;
        using Grid2_1 = typename Patch_1::Grid2;
        using Grid1_2 = typename Patch_2::Grid1;
        using Grid2_2 = typename Patch_2::Grid2;

        using Deriv1_1 = ddc::Deriv<typename Patch_1::Dim1>;
        using Deriv2_1 = ddc::Deriv<typename Patch_1::Dim2>;
        using Deriv1_2 = ddc::Deriv<typename Patch_2::Dim1>;
        using Deriv2_2 = ddc::Deriv<typename Patch_2::Dim2>;

        // Get the fields of the left and right patch of the interface.
        DerivFieldOnPatch_host<Patch_1> function_and_derivs_1
                = functions_and_derivs.template get<Patch_1>();
        DerivFieldOnPatch_host<Patch_2> function_and_derivs_2
                = functions_and_derivs.template get<Patch_2>();

        IdxRange<GridPar1> idx_range_par_1(m_idx_ranges.template get<Patch_1>());
        IdxRange<GridPar2> idx_range_par_2(m_idx_ranges.template get<Patch_2>());

        IdxRange<GridPerp1> idx_range_perp_1(m_idx_ranges.template get<Patch_1>());
        IdxRange<GridPerp2> idx_range_perp_2(m_idx_ranges.template get<Patch_2>());

        // Type for the slice_previous_idx_value index.
        using IdxPatchSlice = std::conditional_t<is_good_orientation, Patch_1, Patch_2>;
        Idx<GridPar1> idx_slice_1;
        Idx<GridPar2> idx_slice_2;
        set_slice_indexes<
                EquivalentInterfaceI,
                IdxPatchSlice>(slice_previous_idx_value, idx_slice_1, idx_slice_2);

        // The slice indices has to be a point at a corner, i.e. index range boundaries.
        assert((idx_slice_1 == idx_range_par_1.front()) || (idx_slice_1 == idx_range_par_1.back()));
        assert((idx_slice_2 == idx_range_par_2.front()) || (idx_slice_2 == idx_range_par_2.back()));

        Idx<DerivPar1, GridPar1> idx_slice_deriv_1(Idx<DerivPar1>(1), idx_slice_1);
        Idx<DerivPar2, GridPar2> idx_slice_deriv_2(Idx<DerivPar2>(1), idx_slice_2);

        DField<IdxRange<GridPerp1>, Kokkos::HostSpace, Kokkos::layout_stride> derivs_1
                = function_and_derivs_1[idx_slice_deriv_1];
        DField<IdxRange<GridPerp2>, Kokkos::HostSpace, Kokkos::layout_stride> derivs_2
                = function_and_derivs_2[idx_slice_deriv_2];

        // Compute the coefficient c_I for the interface I.
        double const corner_derivative
                = m_derivatives_calculators.template get<EquivalentInterfaceI>()
                          .get_derivatives_coefficients(
                                  get_const_field(derivs_1),
                                  get_const_field(derivs_2));

        // Get the correct indices for the slices.
        Idx<GridPerp1> idx_deriv_1 = get_idx_interface(idx_range_perp_1, extremity_1);
        Idx<GridPerp2> idx_deriv_2 = get_idx_interface(idx_range_perp_2, extremity_2);

        // If the parallel grids of the two patches agree.
        constexpr bool are_same_direction_par = EquivalentInterfaceI::orientations_agree;

        constexpr bool change_sign_1 = (extremity_1 == Extremity::FRONT);
        constexpr bool change_sign_2
                = ((extremity_2 == Extremity::BACK) && are_same_direction_par)
                  || ((extremity_2 == Extremity::FRONT) && !are_same_direction_par);

        constexpr int sign_1 = !change_sign_1 - change_sign_1;
        constexpr int sign_2 = !change_sign_2 - change_sign_2;

        // Update the cross-derivatives.
        Idx<Deriv1_1, Grid1_1, Deriv2_1, Grid2_1>
                idx_cross_deriv1(Idx<Deriv1_1>(1), Idx<Deriv2_1>(1), idx_slice_1, idx_deriv_1);
        function_and_derivs_1(idx_cross_deriv1) = corner_derivative * sign_1;

        Idx<Deriv1_2, Grid1_2, Deriv2_2, Grid2_2>
                idx_cross_deriv2(Idx<Deriv1_2>(1), Idx<Deriv2_2>(1), idx_slice_2, idx_deriv_2);
        function_and_derivs_2(idx_cross_deriv2) = corner_derivative * sign_2;
    }


    /**
     * @brief Get the indices on Patch1 and Patch2 of the given interface.
     * The index on Patch1 corresponds to the value given in input
     * (if the interface is well oriented, otherwise it is Patch2).
     * The index on Patch2 is the equivalent index of the one on Patch1.
     * It also updates the value given in input with the value of the index
     * on Patch2 for the next interface 
     * (if the next interface is well oriented, otherwise it is Patch1).
     * 
     * Example for well oriented interfaces:
     * 
     *  |           Patch1 | Patch2      Patch1 | Patch2        |  
     *  |            idx_1 > idx_2        idx_1 > idx_2         |
     *  |              ^   |   v            ^   |   v           |
     *  |idx_val   idx_val | idx_val    idx_val | idx_val       |
     */
    template <
            typename InterfaceI,
            class IdxPatchSlice,
            class ParallGrid1 = typename InterfaceI::Edge1::parallel_grid,
            class ParallGrid2 = typename InterfaceI::Edge2::parallel_grid>
    void set_slice_indexes(
            int& slice_idx_value,
            Idx<ParallGrid1>& slice_idx_1,
            Idx<ParallGrid2>& slice_idx_2)
    {
        using Patch_1 = typename InterfaceI::Edge1::associated_patch;
        using Patch_2 = typename InterfaceI::Edge2::associated_patch;

        IdxRange<ParallGrid1> const idx_range_parall_1(m_idx_ranges.template get<Patch_1>());
        IdxRange<ParallGrid2> const idx_range_parall_2(m_idx_ranges.template get<Patch_2>());

        EdgeTransformation<InterfaceI> index_converter(idx_range_parall_1, idx_range_parall_2);

        if constexpr (std::is_same_v<Patch_1, IdxPatchSlice>) {
            slice_idx_1 = Idx<ParallGrid1>(slice_idx_value);
            slice_idx_2 = index_converter(slice_idx_1);
            slice_idx_value = (slice_idx_2 - idx_range_parall_2.front()).value();
        } else {
            slice_idx_2 = Idx<ParallGrid2>(slice_idx_value);
            slice_idx_1 = index_converter(slice_idx_2);
            slice_idx_value = (slice_idx_1 - idx_range_parall_1.front()).value();
        }
    }

    /**
     * @brief Get the index of the given index range slice corresponding to the 
     * extremity at the interface of the 1D grid. 
     */
    template <class Grid1Or2>
    Idx<Grid1Or2> get_idx_interface(IdxRange<Grid1Or2> const idx_range_slice, Extremity extremity)
            const
    {
        return (extremity == Extremity::FRONT) ? idx_range_slice.front() : idx_range_slice.back();
    }
};
