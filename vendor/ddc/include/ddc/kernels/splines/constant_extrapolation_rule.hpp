// Copyright (C) The DDC development team, see COPYRIGHT.md file
//
// SPDX-License-Identifier: MIT

#pragma once

#include "view.hpp"

namespace ddc {

template <class DimI, class... Dim>
struct ConstantExtrapolationRule
{
};

/**
 * @brief A functor for describing a spline boundary value by a constant extrapolation for 1D evaluator.
 *
 * To define the value of a function on B-splines out of the domain, we here use a constant
 * extrapolation on the edge.
 */
template <class DimI>
struct ConstantExtrapolationRule<DimI>
{
private:
    ddc::Coordinate<DimI> m_eval_pos;

public:
    /**
     * @brief Instantiate a ConstantExtrapolationRule.
     *
     * The boundary value will be the same as at the coordinate eval_pos given.
     *
     * @param[in] eval_pos
     * 			Coordinate inside the domain where we will evaluate each points outside the domain.
     */
    explicit ConstantExtrapolationRule(ddc::Coordinate<DimI> eval_pos) : m_eval_pos(eval_pos) {}

    /**
     * @brief Get the value of the function on B-splines at a coordinate outside the domain.
     *
     * @param[in] pos
     * 			The coordinate where we want to evaluate the function on B-splines.
     * @param[in] spline_coef
     *			The coefficients of the function on B-splines.
     *
     * @return A double with the value of the function on B-splines evaluated at the coordinate.
     */
    template <class CoordType, class BSplines, class Layout, class MemorySpace>
    KOKKOS_FUNCTION double operator()(
            CoordType,
            ddc::ChunkSpan<double const, ddc::DiscreteDomain<BSplines>, Layout, MemorySpace> const
                    spline_coef) const
    {
        static_assert(in_tags_v<DimI, to_type_seq_t<CoordType>>);

        std::array<double, BSplines::degree() + 1> vals;

        ddc::DiscreteElement<BSplines> idx
                = ddc::discrete_space<BSplines>().eval_basis(vals, m_eval_pos);

        double y = 0.0;
        for (std::size_t i = 0; i < BSplines::degree() + 1; ++i) {
            y += spline_coef(idx + i) * vals[i];
        }
        return y;
    }
};

/**
 * @brief A functor for describing a spline boundary value by a constant extrapolation for 2D evaluator.
 *
 * To define the value of a function on B-splines out of the domain, we here use a constant
 * extrapolation on the edge.
 */
template <class DimI, class DimNI>
struct ConstantExtrapolationRule<DimI, DimNI>
{
private:
    ddc::Coordinate<DimI> m_eval_pos;
    ddc::Coordinate<DimNI> m_eval_pos_not_interest_min;
    ddc::Coordinate<DimNI> m_eval_pos_not_interest_max;

public:
    /**
     * @brief Instantiate a ConstantExtrapolationRule.
     *
     * The boundary value will be the same as at the coordinate given in a dimension given.
     * The dimension of the input defines the dimension of the boundary condition.
     * The second and the third parameters are needed in case of non-periodic splines on the
     * dimension off-interest (the complementary dimension of the boundary condition),
     * because the evaluator can receive coordinates outside the domain in both dimension.
     *
     * @param[in] eval_pos_bc
     * 			Coordinate in the dimension given inside the domain where we will evaluate
     * 			each points outside the domain.
     * @param[in] eval_pos_not_interest_min
     * 			The minimum coordinate inside the domain on the complementary dimension of the boundary condition.
     * @param[in] eval_pos_not_interest_max
     * 			The maximum coordinate inside the domain on the complementary dimension of the boundary condition.
     */
    explicit ConstantExtrapolationRule(
            ddc::Coordinate<DimI> eval_pos,
            ddc::Coordinate<DimNI> eval_pos_not_interest_min,
            ddc::Coordinate<DimNI> eval_pos_not_interest_max)
        : m_eval_pos(eval_pos)
        , m_eval_pos_not_interest_min(eval_pos_not_interest_min)
        , m_eval_pos_not_interest_max(eval_pos_not_interest_max)
    {
    }

    /**
     * @brief Instantiate a ConstantExtrapolationRule.
     *
     * The boundary value will be the same as at the coordinate given in a dimension given.
     * The dimension of the input defines the dimension of the boundary condition.
     * No second and third parameters are needed in case of periodic splines on the
     * dimension off-interest (the complementary dimension of the boundary condition).
     *
     * @param[in] eval_pos_bc
     * 			Coordinate in the dimension given inside the domain where we will evaluate
     * 			each points outside the domain.
     */
    template <class DimNI_sfinae = DimNI, std::enable_if_t<DimNI_sfinae::PERIODIC, int> = 0>
    explicit ConstantExtrapolationRule(ddc::Coordinate<DimI> eval_pos)
        : m_eval_pos(eval_pos)
        , m_eval_pos_not_interest_min(0.)
        , m_eval_pos_not_interest_max(0.)
    {
    }

    /**
     * @brief Get the value of the function on B-splines at a coordinate outside the domain.
     *
     * In the dimension defined in the constructor Dim1 (or Dim2), it sets the coordinate pos_1 (or pos_2)
     * given at the m_eval_pos coordinate if it is outside the domain.
     * If the coordinate on the complementary dimension of the boundary condition dimension ddc::select<DimNI>(coord_extrap) is
     * outside the domain, then it also sets the coordinate at eval_pos_not_interest_min
     * (if ddc::select<DimNI>(coord_extrap) @f$ < @f$ eval_pos_not_interest_min) or
     * at eval_pos_not_interest_max (if ddc::select<DimNI>(coord_extrap) @f$ > @f$ eval_pos_not_interest_max).
     *
     * @param[in] coord_extrap
     * 			The coordinates where we want to evaluate the function on B-splines
     * @param[in] spline_coef
     *			The coefficients of the function on B-splines.
     *
     *@return A double with the value of the function on B-splines evaluated at the coordinate.
     */
    template <class CoordType, class BSplines1, class BSplines2, class Layout, class MemorySpace>
    KOKKOS_FUNCTION double operator()(
            CoordType coord_extrap,
            ddc::ChunkSpan<
                    double const,
                    ddc::DiscreteDomain<BSplines1, BSplines2>,
                    Layout,
                    MemorySpace> const spline_coef) const
    {
        static_assert(
                in_tags_v<
                        DimI,
                        to_type_seq_t<CoordType>> && in_tags_v<DimNI, to_type_seq_t<CoordType>>);

        ddc::Coordinate<DimI, DimNI> eval_pos;
        if constexpr (DimNI::PERIODIC) {
            eval_pos = ddc::Coordinate<DimI, DimNI>(m_eval_pos, ddc::select<DimNI>(coord_extrap));
        } else {
            eval_pos = ddc::Coordinate<DimI, DimNI>(
                    m_eval_pos,
                    Kokkos::
                            clamp(ddc::select<DimNI>(coord_extrap),
                                  m_eval_pos_not_interest_min,
                                  m_eval_pos_not_interest_max));
        }

        std::array<double, BSplines1::degree() + 1> vals1;
        std::array<double, BSplines2::degree() + 1> vals2;

        ddc::DiscreteElement<BSplines1> idx1
                = ddc::discrete_space<BSplines1>()
                          .eval_basis(vals1, ddc::select<typename BSplines1::tag_type>(eval_pos));
        ddc::DiscreteElement<BSplines2> idx2
                = ddc::discrete_space<BSplines2>()
                          .eval_basis(vals2, ddc::select<typename BSplines2::tag_type>(eval_pos));

        double y = 0.0;
        for (std::size_t i = 0; i < BSplines1::degree() + 1; ++i) {
            for (std::size_t j = 0; j < BSplines2::degree() + 1; ++j) {
                y += spline_coef(idx1 + i, idx2 + j) * vals1[i] * vals2[j];
            }
        }

        return y;
    }
};
} // namespace ddc
