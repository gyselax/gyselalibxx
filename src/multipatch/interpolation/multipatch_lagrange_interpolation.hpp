// SPDX-License-Identifier: MIT

#pragma once
#include "ddc_aliases.hpp"
#include "lagrange_interpolation.hpp"
#include "multipatch_interpolation_2d.hpp"

namespace detail {

/**
 * @brief A helper to define the 2D Lagrange interpolation used on each patch of a
 * MultipatchLagrangeInterpolator2D.
 *
 * @tparam ExecSpace The Kokkos execution space used for computations.
 * @tparam DataType The type of the function values.
 * @tparam Basis1OnPatch A type template such that Basis1OnPatch<Patch> is the Lagrange basis
 *          along the first dimension of the patch.
 * @tparam Basis2OnPatch A type template such that Basis2OnPatch<Patch> is the Lagrange basis
 *          along the second dimension of the patch.
 * @tparam ExtrapRules1 The extrapolation rules along the first dimension of each patch.
 * @tparam ExtrapRules2 The extrapolation rules along the second dimension of each patch.
 */
template <
        class ExecSpace,
        class DataType,
        template <class>
        class Basis1OnPatch,
        template <class>
        class Basis2OnPatch,
        class ExtrapRules1,
        class ExtrapRules2>
struct LagrangeInterpolatorOnPatch
{
    /// @brief The 2D Lagrange interpolation on the given patch.
    template <class Patch>
    using type = ::LagrangeInterpolator<
            ExecSpace,
            DataType,
            IdxRange<Basis1OnPatch<Patch>, Basis2OnPatch<Patch>>,
            typename Patch::IdxRange12,
            ExtrapRules1,
            ExtrapRules2>;
};

} // namespace detail

/**
 * @brief A MultipatchInterpolation2D which uses a 2D Lagrange interpolation on every patch.
 *
 * On each patch, the Lagrange interpolation is defined on the bases
 * (Basis1OnPatch<Patch>, Basis2OnPatch<Patch>) and the grids (Grid1, Grid2) of the patch.
 * The same extrapolation rules are used on every patch. They are given as tags
 * (e.g. ExtrapolationRule::Periodic, ExtrapolationRule::Null_Null) which are resolved
 * separately on each patch.
 *
 * The execution space is the one of the patch locator (PatchLocator::exec_space).
 *
 * @tparam PatchLocator The operator used to find the patch where a coordinate is physically located.
 * @tparam MultipatchExtrapolationRule The extrapolation rule used outside of the global domain
 *          (e.g. ConstantExtrapolationRuleOnion or NullExtrapolationRule).
 * @tparam Basis1OnPatch A type template such that Basis1OnPatch<Patch> is the Lagrange basis
 *          (UniformLagrangeBasis or NonUniformLagrangeBasis) along the first dimension of the patch.
 * @tparam Basis2OnPatch A type template such that Basis2OnPatch<Patch> is the Lagrange basis
 *          along the second dimension of the patch.
 * @tparam ExtrapRules1 A ddc::TypeSeq<MinExtrapolationRule, MaxExtrapolationRule> describing the
 *          extrapolation along the first dimension of each patch (see LagrangeInterpolator).
 * @tparam ExtrapRules2 A ddc::TypeSeq<MinExtrapolationRule, MaxExtrapolationRule> describing the
 *          extrapolation along the second dimension of each patch (see LagrangeInterpolator).
 * @tparam DataType The type of the function values (default: double).
 */
template <
        class PatchLocator,
        class MultipatchExtrapolationRule,
        template <class>
        class Basis1OnPatch,
        template <class>
        class Basis2OnPatch,
        class ExtrapRules1,
        class ExtrapRules2,
        class DataType = double>
using MultipatchLagrangeInterpolator2D = MultipatchInterpolation2D<
        detail::LagrangeInterpolatorOnPatch<
                typename PatchLocator::exec_space,
                DataType,
                Basis1OnPatch,
                Basis2OnPatch,
                ExtrapRules1,
                ExtrapRules2>::template type,
        PatchLocator,
        MultipatchExtrapolationRule>;
