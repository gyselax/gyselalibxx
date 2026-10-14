// SPDX-License-Identifier: MIT

#pragma once

#include <typeinfo>

#include <ddc/ddc.hpp>

#include "connectivity_details.hpp"
#include "ddc_aliases.hpp"
#include "edge.hpp"
#include "edge_transformation.hpp"
#include "interface.hpp"
#include "interface_derivative_coefficients_collection.hpp"
#include "multipatch_field.hpp"
#include "multipatch_type.hpp"
#include "types.hpp"
#include "view.hpp"

template <class Connectivity, class Grid1D, class PatchSeq, class DerivCoeffCollection>
class InterfacesDerivativeCalculator;

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
  * @tparam DerivCoeffCollection A InterfaceDerivCoeffsCollection
  * that stores the InterfaceDerivCoeffs needed to compute the interface derivatives. 
  */
template <class Connectivity, class Grid1D, class... Patches, class DerivCoeffCollection>
class InterfacesDerivativeCalculator<
        Connectivity,
        Grid1D,
        ddc::detail::TypeSeq<Patches...>,
        DerivCoeffCollection>
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
            is_interface_derivative_coefficients_collection_v<DerivCoeffCollection>,
            "Please provide a InterfaceDerivCoeffsCollection type.");

    static constexpr std::size_t number_of_interfaces
            = ddc::type_seq_size_v<interface_sorted_collection>;

    using InterfaceBoundFirst = ddc::type_seq_element_t<0, interface_sorted_collection>;
    using InterfaceBoundLast
            = ddc::type_seq_element_t<number_of_interfaces - 1, interface_sorted_collection>;

    static_assert(
            (((!std::is_same_v<typename InterfaceBoundFirst::Edge1, OutsideEdge>)
              || (!std::is_same_v<typename InterfaceBoundFirst::Edge2, OutsideEdge>))
             == ((!std::is_same_v<typename InterfaceBoundLast::Edge1, OutsideEdge>)
                 || (!std::is_same_v<typename InterfaceBoundLast::Edge2, OutsideEdge>))),
            "If one global boundary has an outside edge, the other global boundary should have "
            "too.");

    static constexpr bool is_periodic
            = (!std::is_same_v<typename InterfaceBoundFirst::Edge1, OutsideEdge>)&&(
                    !std::is_same_v<typename InterfaceBoundFirst::Edge2, OutsideEdge>);
    /*
        Remove all the interfaces with an OutsideEdge.
        Rely on get_all_interfaces_along_direction_t to order the interfaces.
    */
    using outer_interface_collection = std::conditional_t<
            is_periodic,
            ddc::detail::TypeSeq<>,
            ddc::detail::TypeSeq<InterfaceBoundFirst, InterfaceBoundLast>>;

    // All the interfaces in the Grid1D direction without the interfaces with OutsideEdge.
    using inner_interface_collection
            = ddc::type_seq_remove_t<interface_sorted_collection, outer_interface_collection>;

    // Number of inner interfaces. It fixes the number of elements in the matrix and vectors.
    static constexpr std::size_t n_inner_interfaces
            = ddc::type_seq_size_v<inner_interface_collection>;

    using InterfaceInnerFirst = ddc::type_seq_element_t<0, inner_interface_collection>;

    // First patch of the sorted interfaces.
    using PatchFirst = std::conditional_t<
            is_periodic,
            find_patch_t<Grid1D, all_patches>,
            typename InterfaceInnerFirst::Edge2::associated_patch>;

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

private:
    template <std::size_t I>
    struct InterfaceInfo
    {
        using interface = ddc::type_seq_element_t<I, inner_interface_collection>;
        /// Equivalent interface defined in the main file, i.e. given to MultipatchConnectivity.
        using stored_interface
                = find_associated_interface_t<typename interface::Edge1, all_interface_collection>;

        /// The orientation Patch1|Patch2 of the interface matches with the sorted 1D grid sequence.
        static constexpr bool is_good_orientation = std::is_same_v<
                typename stored_interface::Edge1::perpendicular_grid,
                ddc::type_seq_element_t<I, Grid1DSeq>>;

        using patch1 = typename stored_interface::Edge1::associated_patch;
        using patch2 = typename stored_interface::Edge2::associated_patch;

        static constexpr Extremity extremity_1 = stored_interface::Edge1::extremity;
        static constexpr Extremity extremity_2 = stored_interface::Edge2::extremity;

        using grid_perp1 = typename stored_interface::Edge1::perpendicular_grid;
        using grid_perp2 = typename stored_interface::Edge2::perpendicular_grid;

        using grid_par1 = typename stored_interface::Edge1::parallel_grid;
        using grid_par2 = typename stored_interface::Edge2::parallel_grid;

        using deriv_perp1 = typename ddc::Deriv<typename grid_perp1::continuous_dimension_type>;
        using deriv_perp2 = typename ddc::Deriv<typename grid_perp2::continuous_dimension_type>;

        using deriv_par1 = typename ddc::Deriv<typename grid_par1::continuous_dimension_type>;
        using deriv_par2 = typename ddc::Deriv<typename grid_par2::continuous_dimension_type>;

        using idx_cross_deriv1
                = Idx<ddc::Deriv<typename patch1::Dim1>,
                      typename patch1::Grid1,
                      ddc::Deriv<typename patch1::Dim2>,
                      typename patch1::Grid2>;
        using idx_cross_deriv2
                = Idx<ddc::Deriv<typename patch2::Dim1>,
                      typename patch2::Grid1,
                      ddc::Deriv<typename patch2::Dim2>,
                      typename patch2::Grid2>;

        using idx_deriv_perp1 = Idx<deriv_perp1, typename patch1::Grid1, typename patch1::Grid2>;
        using idx_deriv_perp2 = Idx<deriv_perp2, typename patch2::Grid1, typename patch2::Grid2>;
    };

public:
    /// @brief First patch of the collection of patches in the direction of the given Grid1D.
    using first_patch = PatchFirst;

private:
    MultipatchType<IdxRangeOnPatch, Patches...> const& m_idx_ranges;

    DerivCoeffCollection const& m_deriv_coeffs_collect;

public:
    /**
     * @brief Instantiate InterfacesDerivativeCalculator. 
     *  
     * @param idx_ranges MultipatchType collection of index ranges defined on the given list of patches. 
     * @param deriv_coeffs_collect InterfaceDerivCoeffsCollection containing all the 
     *          interface derivative calculator for each interface in the given Grid1D direction. 
     */
    InterfacesDerivativeCalculator(
            MultipatchType<IdxRangeOnPatch, Patches...> const& idx_ranges,
            DerivCoeffCollection const& deriv_coeffs_collect)
        : m_idx_ranges(idx_ranges)
        , m_deriv_coeffs_collect(deriv_coeffs_collect)
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
                = ddc::in_tags_v<typename PatchFirst::Grid1, Grid1DSeq>;
        using IdxRangeParFirstType = std::conditional_t<
                is_first_idx_range_perp_on_dim1,
                typename PatchFirst::IdxRange2,
                typename PatchFirst::IdxRange1>;

        IdxRangeParFirstType idx_range_par_first(m_idx_ranges.template get<PatchFirst>());

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
                = ddc::in_tags_v<typename PatchFirst::Grid1, Grid1DSeq>;
        using IdxRangeParFirstType = std::conditional_t<
                is_first_idx_range_perp_on_dim1,
                typename PatchFirst::IdxRange2,
                typename PatchFirst::IdxRange1>;

        IdxRangeParFirstType idx_range_par_first(m_idx_ranges.template get<PatchFirst>());

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
                (std::is_same_v<IdxPar, typename PatchFirst::Idx1>)
                        || (std::is_same_v<IdxPar, typename PatchFirst::Idx2>),
                "The given index has to be a 1D index defined on the first patch. See "
                "InterfacesDerivativeCalculator<...>::first_patch.");

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
                (std::is_same_v<GridPar1D, typename PatchFirst::Grid1>)
                || (std::is_same_v<GridPar1D, typename PatchFirst::Grid2>)
                || (std::is_same_v<GridPar1D, ddc::Deriv<typename PatchFirst::Dim1>>)
                || (std::is_same_v<GridPar1D, ddc::Deriv<typename PatchFirst::Dim2>>));

        // Locally transform the index into an integer to avoid type issues.
        IdxRange<GridPar1D> idx_range_1d_first_patch(m_idx_ranges.template get<PatchFirst>());
        int slice_idx_value = (slice_idx - idx_range_1d_first_patch.front()).value();
        (update_derivatives_at_interface<I, eval_type>(functions_and_derivs, slice_idx_value), ...);
    }


    /// @brief Associate the Ith derivative values to the correct first derivative or cross-derivative.
    template <std::size_t I, typename eval_type>
    void update_derivatives_at_interface(
            MultipatchField<DerivFieldOnPatch_host, Patches...> functions_and_derivs,
            int& slice_previous_idx_value)
    {
        using Patch_1 = typename InterfaceInfo<I>::patch1;
        using Patch_2 = typename InterfaceInfo<I>::patch2;

        IdxRange<typename InterfaceInfo<I>::grid_perp1> idx_range_perp_1(
                m_idx_ranges.template get<Patch_1>());
        IdxRange<typename InterfaceInfo<I>::grid_perp2> idx_range_perp_2(
                m_idx_ranges.template get<Patch_2>());

        // Type for the slice_previous_idx_value index.
        Idx<typename InterfaceInfo<I>::grid_par1> idx_slice_1;
        Idx<typename InterfaceInfo<I>::grid_par2> idx_slice_2;
        set_slice_indexes<I>(slice_previous_idx_value, idx_slice_1, idx_slice_2);

        double const reconstructed_deriv = get_reconstructed_deriv<I, eval_type>(
                idx_slice_1,
                idx_slice_2,
                functions_and_derivs.template get<Patch_1>(),
                functions_and_derivs.template get<Patch_2>(),
                m_deriv_coeffs_collect.template get<typename InterfaceInfo<I>::stored_interface>());

        // Get the correct indices for the slices.
        Idx<typename InterfaceInfo<I>::grid_perp1> idx_deriv_1
                = get_idx_interface(idx_range_perp_1, InterfaceInfo<I>::extremity_1);
        Idx<typename InterfaceInfo<I>::grid_perp2> idx_deriv_2
                = get_idx_interface(idx_range_perp_2, InterfaceInfo<I>::extremity_2);

        std::array<const int, 2> signs = get_signs<I, eval_type>();

        // Update the derivatives or cross-derivatives.
        Idx<typename InterfaceInfo<I>::deriv_perp1> i_dperp1(1);
        Idx<typename InterfaceInfo<I>::deriv_perp2> i_dperp2(1);
        if constexpr (std::is_same_v<eval_type, eval_deriv>) {
            typename InterfaceInfo<I>::idx_deriv_perp1
                    idx_deriv1(i_dperp1, idx_deriv_1, idx_slice_1);
            typename InterfaceInfo<I>::idx_deriv_perp2
                    idx_deriv2(i_dperp2, idx_deriv_2, idx_slice_2);

            functions_and_derivs.template get<Patch_1>()(idx_deriv1)
                    = reconstructed_deriv * signs[0];
            functions_and_derivs.template get<Patch_2>()(idx_deriv2)
                    = reconstructed_deriv * signs[1];
        } else {
            Idx<typename InterfaceInfo<I>::deriv_par1> i_dpar1(1);
            Idx<typename InterfaceInfo<I>::deriv_par2> i_dpar2(1);

            typename InterfaceInfo<I>::idx_cross_deriv1
                    idx_deriv1(i_dperp1, i_dpar1, idx_deriv_1, idx_slice_1);
            typename InterfaceInfo<I>::idx_cross_deriv2
                    idx_deriv2(i_dperp2, i_dpar2, idx_deriv_2, idx_slice_2);

            functions_and_derivs.template get<Patch_1>()(idx_deriv1)
                    = reconstructed_deriv * signs[0];
            functions_and_derivs.template get<Patch_2>()(idx_deriv2)
                    = reconstructed_deriv * signs[1];
        }
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
    template <std::size_t I>
    void set_slice_indexes(
            int& slice_idx_value,
            Idx<typename InterfaceInfo<I>::grid_par1>& slice_idx_1,
            Idx<typename InterfaceInfo<I>::grid_par2>& slice_idx_2)
    {
        using Patch_1 = typename InterfaceInfo<I>::patch1;
        using Patch_2 = typename InterfaceInfo<I>::patch2;

        using ParallGrid1 = typename InterfaceInfo<I>::grid_par1;
        using ParallGrid2 = typename InterfaceInfo<I>::grid_par2;

        using PatchIdxSlice
                = std::conditional_t<InterfaceInfo<I>::is_good_orientation, Patch_1, Patch_2>;

        IdxRange<ParallGrid1> const idx_range_parall_1(m_idx_ranges.template get<Patch_1>());
        IdxRange<ParallGrid2> const idx_range_parall_2(m_idx_ranges.template get<Patch_2>());

        EdgeTransformation<typename InterfaceInfo<I>::stored_interface>
                index_converter(idx_range_parall_1, idx_range_parall_2);

        if constexpr (std::is_same_v<Patch_1, PatchIdxSlice>) {
            slice_idx_1 = Idx<ParallGrid1>(slice_idx_value);
            slice_idx_2 = index_converter(slice_idx_1);
            slice_idx_value = (slice_idx_2 - idx_range_parall_2.front()).value();
        } else {
            slice_idx_2 = Idx<ParallGrid2>(slice_idx_value);
            slice_idx_1 = index_converter(slice_idx_2);
            slice_idx_value = (slice_idx_1 - idx_range_parall_1.front()).value();
        }
    }

    template <
            std::size_t I,
            typename eval_type,
            class InterfaceDerivCoeffs,
            std::enable_if_t<std::is_same_v<eval_type, eval_deriv>, bool> = true>
    double const get_reconstructed_deriv(
            Idx<typename InterfaceInfo<I>::grid_par1> const& idx_slice_1,
            Idx<typename InterfaceInfo<I>::grid_par2> const& idx_slice_2,
            DerivFieldOnPatch_host<typename InterfaceInfo<I>::patch1> functions_and_derivs_1,
            DerivFieldOnPatch_host<typename InterfaceInfo<I>::patch2> functions_and_derivs_2,
            InterfaceDerivCoeffs const& deriv_coeffs)
    {
        // Use the function values to compute the first derivatives.
        DField<typename InterfaceInfo<I>::patch1::IdxRange12,
               Kokkos::HostSpace,
               Kokkos::layout_stride>
                function_1 = functions_and_derivs_1.get_values_field();
        DField<typename InterfaceInfo<I>::patch2::IdxRange12,
               Kokkos::HostSpace,
               Kokkos::layout_stride>
                function_2 = functions_and_derivs_2.get_values_field();

        // Compute the coefficient c_I for the interface I.
        double const interface_deriv = deriv_coeffs.get_approx_deriv(
                get_const_field(function_1[idx_slice_1]),
                get_const_field(function_2[idx_slice_2]));
        return interface_deriv;
    }

    template <
            std::size_t I,
            typename eval_type,
            class InterfaceDerivCoeffs,
            std::enable_if_t<std::is_same_v<eval_type, eval_cross_deriv>, bool> = true>
    double const get_reconstructed_deriv(
            Idx<typename InterfaceInfo<I>::grid_par1> const& idx_slice_1,
            Idx<typename InterfaceInfo<I>::grid_par2> const& idx_slice_2,
            DerivFieldOnPatch_host<typename InterfaceInfo<I>::patch1> functions_and_derivs_1,
            DerivFieldOnPatch_host<typename InterfaceInfo<I>::patch2> functions_and_derivs_2,
            InterfaceDerivCoeffs const& deriv_coeffs)
    {
        IdxRange<typename InterfaceInfo<I>::grid_par1> idx_range_par_1(
                m_idx_ranges.template get<typename InterfaceInfo<I>::patch1>());
        IdxRange<typename InterfaceInfo<I>::grid_par2> idx_range_par_2(
                m_idx_ranges.template get<typename InterfaceInfo<I>::patch2>());

        // The slice indices has to be a point at a corner, i.e. index range boundaries.
        assert((idx_slice_1 == idx_range_par_1.front()) || (idx_slice_1 == idx_range_par_1.back()));
        assert((idx_slice_2 == idx_range_par_2.front()) || (idx_slice_2 == idx_range_par_2.back()));

        Idx<typename InterfaceInfo<I>::deriv_par1, typename InterfaceInfo<I>::grid_par1>
                idx_slice_deriv_1(Idx<typename InterfaceInfo<I>::deriv_par1>(1), idx_slice_1);
        Idx<typename InterfaceInfo<I>::deriv_par2, typename InterfaceInfo<I>::grid_par2>
                idx_slice_deriv_2(Idx<typename InterfaceInfo<I>::deriv_par2>(1), idx_slice_2);

        DField<IdxRange<typename InterfaceInfo<I>::grid_perp1>,
               Kokkos::HostSpace,
               Kokkos::layout_stride>
                derivs_1 = functions_and_derivs_1[idx_slice_deriv_1];
        DField<IdxRange<typename InterfaceInfo<I>::grid_perp2>,
               Kokkos::HostSpace,
               Kokkos::layout_stride>
                derivs_2 = functions_and_derivs_2[idx_slice_deriv_2];

        // Compute the coefficient c_I for the interface I.
        double const corner_derivative = deriv_coeffs.get_approx_cross_deriv(
                get_const_field(derivs_1),
                get_const_field(derivs_2));
        return corner_derivative;
    }

    /**
     * @brief Get the index of the given index range slice corresponding to the 
     * extremity at the interface of the 1D grid. 
     */
    template <class Grid1Or2>
    static Idx<Grid1Or2> get_idx_interface(
            IdxRange<Grid1Or2> const idx_range_slice,
            Extremity extremity)
    {
        return (extremity == Extremity::FRONT) ? idx_range_slice.front() : idx_range_slice.back();
    }

    /**
     * @brief Get the signs of the correction on the derivatives. 
     * The sign of the computed interface derivative can be wrong in the case
     * the axes do not agree. 
     */
    template <std::size_t I, typename eval_type>
    static std::array<const int, 2> get_signs()
    {
        // If the parallel grids of the two patches agree.
        constexpr bool are_same_direction_par
                = std::is_same_v<eval_type, eval_deriv>
                          ? true
                          : InterfaceInfo<I>::interface::orientations_agree;
        constexpr Extremity extremity_1 = InterfaceInfo<I>::extremity_1;
        constexpr Extremity extremity_2 = InterfaceInfo<I>::extremity_2;

        constexpr bool change_sign_1 = (extremity_1 == Extremity::FRONT);
        constexpr bool change_sign_2
                = ((extremity_2 == Extremity::BACK) && are_same_direction_par)
                  || ((extremity_2 == Extremity::FRONT) && !are_same_direction_par);

        constexpr int sign_1 = !change_sign_1 - change_sign_1;
        constexpr int sign_2 = !change_sign_2 - change_sign_2;
        return std::array<const int, 2> {sign_1, sign_2};
    }
};
