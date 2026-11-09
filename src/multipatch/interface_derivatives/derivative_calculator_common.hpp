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
class DerivativeCalculatorCommon;

/**
  * @brief Abstract class to hold the interface derivative coefficient collection
  * and the index ranges.
  * This is the superclass for interface derivative calculators defined in this folder.
  * 
  * @tparam Connectivity A MultipatchConnectivity class describing all the patch connections.
  * @tparam Grid1D A given direction.
  * @tparam Patches List of patches containing all the involved patches in the given direction. 
  * @tparam DerivCoeffCollection A InterfaceDerivCoeffsCollection
  * that stores the InterfaceDerivCoeffs needed to compute the interface derivatives. 
  */
template <class Connectivity, class Grid1D, class... Patches, class DerivCoeffCollection>
class DerivativeCalculatorCommon<
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
    /// @brief Tag for all the sorted interfaces with the correct patch order.
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
    /// @brief Boolean to indicate the periodicty of the global domain.
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
    /// @brief Tag for all the interfaces in the Grid1D direction without the interfaces with OutsideEdge.
    using inner_interface_collection
            = ddc::type_seq_remove_t<interface_sorted_collection, outer_interface_collection>;

    /// @brief Number of inner interfaces (fixing the number of elements in the matrix and vectors).
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
    /// @brief Sequence of 1D grids in the direction of the given Grid1D over the patches.
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

protected:
    /// @brief Tag the interface with the correct patch order.
    template <std::size_t I>
    using interface = ddc::type_seq_element_t<I, inner_interface_collection>;

    /// @brief Tag the equivalent interface defined in the main file, i.e. given to MultipatchConnectivity.
    template <std::size_t I>
    using stored_interface
            = find_associated_interface_t<typename interface<I>::Edge1, all_interface_collection>;

    /// @brief Boolean if the order of Patch1|Patch2 of the interface matches with the sorted 1D grid sequence.
    template <std::size_t I>
    static constexpr bool patch_order_matches_interfaces = std::is_same_v<
            typename stored_interface<I>::Edge1::perpendicular_grid,
            ddc::type_seq_element_t<I, Grid1DSeq>>;

    /// @brief Tag the first patch of the stored interface.
    template <std::size_t I>
    using patch1 = typename stored_interface<I>::Edge1::associated_patch;
    /// @brief Tag the second patch of the stored interface.
    template <std::size_t I>
    using patch2 = typename stored_interface<I>::Edge2::associated_patch;

    /// @brief Tag the connection extremity of the first patch of the interface.
    template <std::size_t I>
    static constexpr Extremity extremity_1 = stored_interface<I>::Edge1::extremity;
    /// @brief Tag the connection extremity of the second patch of the interface.
    template <std::size_t I>
    static constexpr Extremity extremity_2 = stored_interface<I>::Edge2::extremity;

    /// @brief Tag the perpendicular grid of the first patch to the interface.
    template <std::size_t I>
    using grid_perp1 = typename stored_interface<I>::Edge1::perpendicular_grid;
    /// @brief Tag the perpendicular grid of the second patch to the interface.
    template <std::size_t I>
    using grid_perp2 = typename stored_interface<I>::Edge2::perpendicular_grid;

    /// @brief Tag the parallel grid of the first patch to the interface.
    template <std::size_t I>
    using grid_par1 = typename stored_interface<I>::Edge1::parallel_grid;
    /// @brief Tag the parallel grid of the second patch to the interface.
    template <std::size_t I>
    using grid_par2 = typename stored_interface<I>::Edge2::parallel_grid;

    /// @brief Tag the derivative dimension along perpendicular grid of the first patch.
    template <std::size_t I>
    using deriv_perp1 = typename ddc::Deriv<typename grid_perp1<I>::continuous_dimension_type>;
    /// @brief Tag the derivative dimension along perpendicular grid of the second patch.
    template <std::size_t I>
    using deriv_perp2 = typename ddc::Deriv<typename grid_perp2<I>::continuous_dimension_type>;

    /// @brief Tag the derivative dimension along parallel grid of the first patch.
    template <std::size_t I>
    using deriv_par1 = typename ddc::Deriv<typename grid_par1<I>::continuous_dimension_type>;
    /// @brief Tag the derivative dimension along parallel grid of the second patch.
    template <std::size_t I>
    using deriv_par2 = typename ddc::Deriv<typename grid_par2<I>::continuous_dimension_type>;

    /// @brief Tag the index type for cross-derivative on the first patch.
    template <std::size_t I>
    using idx_cross_deriv1
            = Idx<ddc::Deriv<typename patch1<I>::Dim1>,
                  typename patch1<I>::Grid1,
                  ddc::Deriv<typename patch1<I>::Dim2>,
                  typename patch1<I>::Grid2>;
    /// @brief Tag the index type for cross-derivatives on the second patch.
    template <std::size_t I>
    using idx_cross_deriv2
            = Idx<ddc::Deriv<typename patch2<I>::Dim1>,
                  typename patch2<I>::Grid1,
                  ddc::Deriv<typename patch2<I>::Dim2>,
                  typename patch2<I>::Grid2>;

    /// @brief Tag the index type for derivatives along the perpendicular grid on the first patch.
    template <std::size_t I>
    using idx_deriv_perp1
            = Idx<deriv_perp1<I>, typename patch1<I>::Grid1, typename patch1<I>::Grid2>;
    /// @brief Tag the index type for derivatives along the perpendicular grid on the second patch.
    template <std::size_t I>
    using idx_deriv_perp2
            = Idx<deriv_perp2<I>, typename patch2<I>::Grid1, typename patch2<I>::Grid2>;

public:
    /// @brief First patch of the collection of patches in the direction of the given Grid1D.
    using patch_first = PatchFirst;

    /// @brief Grid 1D given to the class indicating the direction of the mesh line along which the derivatives are computed.
    using grid_direction = Grid1D;

    /// @brief Patch of the given direction grid.
    using patch_direction = find_patch_t<Grid1D, all_patches>;

    /// @brief Derivative dimension of the direction patch perpendicular to the mesh line.
    using deriv_perp_dim = ddc::Deriv<typename Grid1D::continuous_dimension_type>;
    /// @brief Derivative dimension of the direction patch parallel to the mesh line.
    using deriv_par_dim = ddc::Deriv<std::conditional_t<
            std::is_same_v<typename patch_direction::Grid1, grid_direction>,
            typename patch_direction::Dim2,
            typename patch_direction::Dim1>>;

    /// @brief Index on the derivative dimension perpendicular to the mesh line.
    using idx_deriv_perp = Idx<deriv_perp_dim>;
    /// @brief Index on the two derivative dimensions of patch_direction.
    using idx_cross_deriv = Idx<deriv_perp_dim, deriv_par_dim>;

protected:
    /// @brief MultipatchType on IdxRangeOnPatch to collect the index ranges on the patches.
    MultipatchType<IdxRangeOnPatch, Patches...> const& m_idx_ranges;

    /// @brief InterfaceDerivCoeffsCollection to collect all the needed InterfaceDerivCoeffs.
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

protected:
    /**
     * @brief Compute the approximated derivative (or coefficient "c") at the interface along the mesh line
     * indexed at the given indices. 
     * 
     * @tparam I Index of the interface. 
     * @param[in] idx_to_slice_1 Index on the parallel grid of the first patch to select the mesh line along 
     *              which the derivative is computed. 
     * @param[in] idx_to_slice_2 Index on the parallel grid of the second patch to select the mesh line along 
     *              which the derivative is computed. 
     * @param[in] functions_and_derivs Collection fields with its first derivatives and cross-derivatives.
     * @return The approximated derivative at the interface. 
     */
    template <std::size_t I, bool is_first_deriv, std::enable_if_t<is_first_deriv, bool> = true>
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

    /**
     * @brief Compute the approximated cross-derivative (or coefficient "c") at the interface along the mesh line
     * indexed at the given indices. 
     * Work only for indices in the corners of the patches (where the cross-derivatives are supposed to be 
     * computed).
     * 
     * @tparam I Index of the interface. 
     * @param[in] idx_to_slice_1 Index on the parallel grid of the first patch to select the mesh line along 
     *              which the derivative is computed.
     * @param[in] idx_to_slice_2 Index on the parallel grid of the second patch to select the mesh line along 
     *              which the derivative is computed. 
     * @param[in] functions_and_derivs Collection fields with its first derivatives and cross-derivatives.
     * @return The approximated cross-derivative at the interface. 
     */
    template <std::size_t I, bool is_first_deriv, std::enable_if_t<!is_first_deriv, bool> = true>
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
     * 
     * @param[in,out] value_idx_to_slice Integer containing the value of the index on the parallel grid 
     *              of the last interface.
     * @param[out] idx_to_slice_1 Index on the parallel grid of the first patch to select the mesh line along 
     *              which the derivative will be computed.
     * @param[out] idx_to_slice_2 Index on the parallel grid of the second patch to select the mesh line along 
     *              which the derivative will be computed. 
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
