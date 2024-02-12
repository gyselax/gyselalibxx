#pragma once

#include "view.hpp"

namespace ddc {

template <class DimI, class... Dim>
struct ConstantExtrapolationRule
{
};

template <class DimI>
struct ConstantExtrapolationRule<DimI>
{
private:
    ddc::Coordinate<DimI> m_eval_pos;

public:
    explicit ConstantExtrapolationRule(ddc::Coordinate<DimI> eval_pos) : m_eval_pos(eval_pos) {}

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

template <class DimI, class DimNI>
struct ConstantExtrapolationRule<DimI, DimNI>
{
private:
    ddc::Coordinate<DimI> m_eval_pos;
    ddc::Coordinate<DimNI> m_eval_pos_not_interest_min;
    ddc::Coordinate<DimNI> m_eval_pos_not_interest_max;

public:
    explicit ConstantExtrapolationRule(
            ddc::Coordinate<DimI> eval_pos,
            ddc::Coordinate<DimNI> eval_pos_not_interest_min,
            ddc::Coordinate<DimNI> eval_pos_not_interest_max)
        : m_eval_pos(eval_pos)
        , m_eval_pos_not_interest_min(eval_pos_not_interest_min)
        , m_eval_pos_not_interest_max(eval_pos_not_interest_max)
    {
    }

    template <class DimNI_sfinae = DimNI, std::enable_if_t<DimNI_sfinae::PERIODIC, int> = 0>
    explicit ConstantExtrapolationRule(ddc::Coordinate<DimI> eval_pos)
        : m_eval_pos(eval_pos)
        , m_eval_pos_not_interest_min(0.)
        , m_eval_pos_not_interest_max(0.)
    {
    }

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
