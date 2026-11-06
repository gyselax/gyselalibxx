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

template <
        class Derived,
        class Connectivity,
        class Grid1D,
        class PatchSeq,
        class DerivCoeffCollection>
class DerivativeCalculatorCommon;

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
template <
        class Derived,
        class Connectivity,
        class Grid1D,
        class... Patches,
        class DerivCoeffCollection>
class DerivativeCalculatorCommon<
        Derived,
        Connectivity,
        Grid1D,
        ddc::TypeSeq<Patches...>,
        DerivCoeffCollection>
{
    /*
        All the interfaces given as input to the MultipatchConnectivity class.
         We expect the parameters defined on these interfaces.
    */
    using all_interface_collection = typename Connectivity::interface_collection;

protected:
    // All the sorted interfaces with the correct orientation.
    using interface_sorted_collection =
            typename Connectivity::template get_all_interfaces_along_direction_t<Grid1D>;

private:
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

protected:
    static constexpr bool is_periodic
            = (!std::is_same_v<typename InterfaceBoundFirst::Edge1, OutsideEdge>)&&(
                    !std::is_same_v<typename InterfaceBoundFirst::Edge2, OutsideEdge>);

private:
    /*
        Remove all the interfaces with an OutsideEdge.
        Rely on get_all_interfaces_along_direction_t to order the interfaces.
    */
    using outer_interface_collection = std::conditional_t<
            is_periodic,
            ddc::TypeSeq<>,
            ddc::TypeSeq<InterfaceBoundFirst, InterfaceBoundLast>>;

protected:
    // All the interfaces in the Grid1D direction without the interfaces with OutsideEdge.
    using inner_interface_collection
            = ddc::type_seq_remove_t<interface_sorted_collection, outer_interface_collection>;

    // Number of inner interfaces. It fixes the number of elements in the matrix and vectors.
    static constexpr std::size_t n_inner_interfaces
            = ddc::type_seq_size_v<inner_interface_collection>;

private:
    using InterfaceInnerFirst = ddc::type_seq_element_t<0, inner_interface_collection>;

    // First patch of the sorted interfaces.
    using PatchFirst = std::conditional_t<
            is_periodic,
            find_patch_t<Grid1D, all_patches>,
            typename InterfaceInnerFirst::Edge2::associated_patch>;

protected:
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
                            ddc::TypeSeq<typename Patches::Grid1...>,
                            ddc::TypeSeq<typename Patches::Grid2...>>>,
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

protected:
    // Interface with the correct patch order.
    template <std::size_t I>
    using interface = ddc::type_seq_element_t<I, inner_interface_collection>;

    /// Equivalent interface defined in the main file, i.e. given to MultipatchConnectivity.
    template <std::size_t I>
    using stored_interface
            = find_associated_interface_t<typename interface<I>::Edge1, all_interface_collection>;

    // The order of Patch1|Patch2 of the interface matches with the sorted 1D grid sequence.
    template <std::size_t I>
    static constexpr bool patch_order_matches_interfaces = std::is_same_v<
            typename stored_interface<I>::Edge1::perpendicular_grid,
            ddc::type_seq_element_t<I, Grid1DSeq>>;

    template <std::size_t I>
    using patch1 = typename stored_interface<I>::Edge1::associated_patch;
    template <std::size_t I>
    using patch2 = typename stored_interface<I>::Edge2::associated_patch;

    template <std::size_t I>
    static constexpr Extremity extremity_1 = stored_interface<I>::Edge1::extremity;
    template <std::size_t I>
    static constexpr Extremity extremity_2 = stored_interface<I>::Edge2::extremity;

    template <std::size_t I>
    using grid_perp1 = typename stored_interface<I>::Edge1::perpendicular_grid;
    template <std::size_t I>
    using grid_perp2 = typename stored_interface<I>::Edge2::perpendicular_grid;

    template <std::size_t I>
    using grid_par1 = typename stored_interface<I>::Edge1::parallel_grid;
    template <std::size_t I>
    using grid_par2 = typename stored_interface<I>::Edge2::parallel_grid;

    template <std::size_t I>
    using deriv_perp1 = typename ddc::Deriv<typename grid_perp1<I>::continuous_dimension_type>;
    template <std::size_t I>
    using deriv_perp2 = typename ddc::Deriv<typename grid_perp2<I>::continuous_dimension_type>;

    template <std::size_t I>
    using deriv_par1 = typename ddc::Deriv<typename grid_par1<I>::continuous_dimension_type>;
    template <std::size_t I>
    using deriv_par2 = typename ddc::Deriv<typename grid_par2<I>::continuous_dimension_type>;

    template <std::size_t I>
    using idx_cross_deriv1
            = Idx<ddc::Deriv<typename patch1<I>::Dim1>,
                  typename patch1<I>::Grid1,
                  ddc::Deriv<typename patch1<I>::Dim2>,
                  typename patch1<I>::Grid2>;
    template <std::size_t I>
    using idx_cross_deriv2
            = Idx<ddc::Deriv<typename patch2<I>::Dim1>,
                  typename patch2<I>::Grid1,
                  ddc::Deriv<typename patch2<I>::Dim2>,
                  typename patch2<I>::Grid2>;

    template <std::size_t I>
    using idx_deriv_perp1
            = Idx<deriv_perp1<I>, typename patch1<I>::Grid1, typename patch1<I>::Grid2>;
    template <std::size_t I>
    using idx_deriv_perp2
            = Idx<deriv_perp2<I>, typename patch2<I>::Grid1, typename patch2<I>::Grid2>;

public:
    /// @brief First patch of the collection of patches in the direction of the given Grid1D.
    using patch_first = PatchFirst;

protected:
    MultipatchType<IdxRangeOnPatch, Patches...> const& m_idx_ranges;

    DerivCoeffCollection const& m_deriv_coeffs_collect;

public:
    /**
     * @brief Instantiate DerivativeCalculatorCommon. 
     *  
     * @param idx_ranges MultipatchType collection of index ranges defined on the given list of patches. 
     * @param deriv_coeffs_collect InterfaceDerivCoeffsCollection containing all the 
     *          interface derivative calculator for each interface in the given Grid1D direction. 
     */
    DerivativeCalculatorCommon(
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
        derived().template solve<eval_deriv>(functions_and_derivs, idx_par);
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
                    derived().template solve<eval_deriv>(functions_and_derivs, idx_par);
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
        derived().template solve<eval_cross_deriv>(functions_and_derivs, idx_par);
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

        derived()
                .template solve<
                        eval_cross_deriv>(functions_and_derivs, idx_range_par_first.front());
        derived()
                .template solve<eval_cross_deriv>(functions_and_derivs, idx_range_par_first.back());
    }

private:
    Derived& derived()
    {
        return static_cast<Derived&>(*this);
    }

protected:
    template <
            std::size_t I,
            typename eval_type,
            std::enable_if_t<std::is_same_v<eval_type, eval_deriv>, bool> = true>
    double const get_approximated_deriv(
            Idx<grid_par1<I>> const& idx_to_slice_1,
            Idx<grid_par2<I>> const& idx_to_slice_2,
            MultipatchField<DerivFieldOnPatch_host, Patches...> functions_and_derivs)
    {
        // Use the function values to compute the first derivatives.
        DField<typename patch1<I>::IdxRange12, Kokkos::HostSpace, Kokkos::layout_stride> function_1
                = functions_and_derivs.template get<patch1<I>>().get_values_field();
        DField<typename patch2<I>::IdxRange12, Kokkos::HostSpace, Kokkos::layout_stride> function_2
                = functions_and_derivs.template get<patch2<I>>().get_values_field();

        // Compute the coefficient c_I for the interface I.
        double const interface_deriv
                = m_deriv_coeffs_collect.template get<stored_interface<I>>().get_approx_deriv(
                        get_const_field(function_1[idx_to_slice_1]),
                        get_const_field(function_2[idx_to_slice_2]));
        return interface_deriv;
    }

    template <
            std::size_t I,
            typename eval_type,
            std::enable_if_t<std::is_same_v<eval_type, eval_cross_deriv>, bool> = true>
    double const get_approximated_deriv(
            Idx<grid_par1<I>> const& idx_to_slice_1,
            Idx<grid_par2<I>> const& idx_to_slice_2,
            MultipatchField<DerivFieldOnPatch_host, Patches...> functions_and_derivs)
    {
        IdxRange<grid_par1<I>> idx_range_par_1(m_idx_ranges.template get<patch1<I>>());
        IdxRange<grid_par2<I>> idx_range_par_2(m_idx_ranges.template get<patch2<I>>());

        // The slice indices has to be a point at a corner, i.e. index range boundaries.
        assert((idx_to_slice_1 == idx_range_par_1.front())
               || (idx_to_slice_1 == idx_range_par_1.back()));
        assert((idx_to_slice_2 == idx_range_par_2.front())
               || (idx_to_slice_2 == idx_range_par_2.back()));

        Idx<deriv_par1<I>, grid_par1<I>>
                idx_to_slice_deriv_1(Idx<deriv_par1<I>>(1), idx_to_slice_1);
        Idx<deriv_par2<I>, grid_par2<I>>
                idx_to_slice_deriv_2(Idx<deriv_par2<I>>(1), idx_to_slice_2);

        DField<IdxRange<grid_perp1<I>>, Kokkos::HostSpace, Kokkos::layout_stride> derivs_1
                = functions_and_derivs.template get<patch1<I>>()[idx_to_slice_deriv_1];
        DField<IdxRange<grid_perp2<I>>, Kokkos::HostSpace, Kokkos::layout_stride> derivs_2
                = functions_and_derivs.template get<patch2<I>>()[idx_to_slice_deriv_2];

        // Compute the coefficient c_I for the interface I.
        double const corner_derivative
                = m_deriv_coeffs_collect.template get<stored_interface<I>>().get_approx_cross_deriv(
                        get_const_field(derivs_1),
                        get_const_field(derivs_2));
        return corner_derivative;
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
            int& value_idx_to_slice,
            Idx<grid_par1<I>>& idx_to_slice_1,
            Idx<grid_par2<I>>& idx_to_slice_2)
    {
        // The Patch with the idx_to_slice value corresponding to the given value_idx_to_slice.
        using PatchIndexSlice
                = std::conditional_t<patch_order_matches_interfaces<I>, patch1<I>, patch2<I>>;

        IdxRange<grid_par1<I>> const idx_range_parall_1(m_idx_ranges.template get<patch1<I>>());
        IdxRange<grid_par2<I>> const idx_range_parall_2(m_idx_ranges.template get<patch2<I>>());

        EdgeTransformation<stored_interface<I>>
                index_converter(idx_range_parall_1, idx_range_parall_2);

        if constexpr (std::is_same_v<patch1<I>, PatchIndexSlice>) {
            idx_to_slice_1 = Idx<grid_par1<I>>(value_idx_to_slice);
            idx_to_slice_2 = index_converter(idx_to_slice_1);
            value_idx_to_slice = (idx_to_slice_2 - idx_range_parall_2.front()).value();
        } else {
            idx_to_slice_2 = Idx<grid_par2<I>>(value_idx_to_slice);
            idx_to_slice_1 = index_converter(idx_to_slice_2);
            value_idx_to_slice = (idx_to_slice_1 - idx_range_parall_1.front()).value();
        }
    }

    /**
     * @brief Get the index of the given index range slice corresponding to one of the 
     * extremities of the 1D grid. 
     * If is_on_interface is true, the extremity corresponds to the interface. 
     * If is_on_interface is false, the extremity corresponds to the other bound of the interface
     * (next or previous). 
     */
    template <
            bool is_on_interface,
            std::size_t I,
            class Patch,
            bool is_patch1 = std::is_same_v<Patch, patch1<I>>>
    std::conditional_t<is_patch1, Idx<grid_perp1<I>>, Idx<grid_perp2<I>>> get_idx_interface()
    {
        using IdxRangePerp
                = std::conditional_t<is_patch1, IdxRange<grid_perp1<I>>, IdxRange<grid_perp2<I>>>;
        IdxRangePerp idx_range_slice(m_idx_ranges.template get<Patch>());
        Extremity extremity = is_patch1 ? extremity_1<I> : extremity_2<I>;
        Extremity extremity_test = is_on_interface ? Extremity::FRONT : Extremity::BACK;
        return (extremity == extremity_test) ? idx_range_slice.front() : idx_range_slice.back();
    }
};
