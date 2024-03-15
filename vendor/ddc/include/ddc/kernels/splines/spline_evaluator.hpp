#pragma once

#include <array>

#include <ddc/ddc.hpp>

#include "Kokkos_Macros.hpp"
#include "view.hpp"

namespace ddc {

template <
        class ExecSpace,
        class MemorySpace,
        class BSplinesType,
        class InterpolationMesh,
        class LeftExtrapolationRule,
        class RightExtrapolationRule,
        class... IDimX>
class SplineEvaluator
{
private:
    // Tags to determine what to evaluate
    struct eval_type
    {
    };

    struct eval_deriv_type
    {
    };

    using tag_type = typename BSplinesType::tag_type;

public:
    using exec_space = ExecSpace;

    using memory_space = MemorySpace;

    using bsplines_type = BSplinesType;

    using interpolation_mesh_type = InterpolationMesh;

    using interpolation_domain_type = ddc::DiscreteDomain<interpolation_mesh_type>;

    using vals_domain_type = ddc::DiscreteDomain<IDimX...>;

    using bsplines_domain_type = ddc::DiscreteDomain<bsplines_type>;

    using batch_domain_type =
            typename ddc::detail::convert_type_seq_to_discrete_domain<ddc::type_seq_remove_t<
                    ddc::detail::TypeSeq<IDimX...>,
                    ddc::detail::TypeSeq<interpolation_mesh_type>>>;

    template <typename Tag>
    using spline_dim_type
            = std::conditional_t<std::is_same_v<Tag, interpolation_mesh_type>, bsplines_type, Tag>;

    using spline_domain_type =
            typename ddc::detail::convert_type_seq_to_discrete_domain<ddc::type_seq_replace_t<
                    ddc::detail::TypeSeq<IDimX...>,
                    ddc::detail::TypeSeq<interpolation_mesh_type>,
                    ddc::detail::TypeSeq<bsplines_type>>>;


private:
    const spline_domain_type m_spline_domain;

    LeftExtrapolationRule m_left_extrap_rule;

    RightExtrapolationRule m_right_bc;

public:
    static_assert(
            std::is_same_v<LeftExtrapolationRule,
                            typename ddc::PeriodicExtrapolationRule<
                                    tag_type>> == bsplines_type::is_periodic()
                    && std::is_same_v<
                               RightExtrapolationRule,
                               typename ddc::PeriodicExtrapolationRule<
                                       tag_type>> == bsplines_type::is_periodic(),
            "PeriodicExtrapolationRule has to be used if and only if dimension is periodic");
    static_assert(
            std::is_invocable_r_v<
                    double,
                    LeftExtrapolationRule,
                    ddc::Coordinate<tag_type>,
                    ddc::ChunkSpan<
                            double const,
                            bsplines_domain_type,
                            std::experimental::layout_right,
                            memory_space>>,
            "LeftExtrapolationRule::operator() has to be callable with usual arguments.");
    static_assert(
            std::is_invocable_r_v<
                    double,
                    RightExtrapolationRule,
                    ddc::Coordinate<tag_type>,
                    ddc::ChunkSpan<
                            double const,
                            bsplines_domain_type,
                            std::experimental::layout_right,
                            memory_space>>,
            "RightExtrapolationRule::operator() has to be callable with usual arguments.");

    explicit SplineEvaluator(
            spline_domain_type const& spline_domain,
            LeftExtrapolationRule const& left_extrap_rule,
            RightExtrapolationRule const& right_extrap_rule)
        : m_spline_domain(spline_domain)
        , m_left_extrap_rule(left_extrap_rule)
        , m_right_bc(right_extrap_rule)
    {
    }

    SplineEvaluator(SplineEvaluator const& x) = default;

    SplineEvaluator(SplineEvaluator&& x) = default;

    ~SplineEvaluator() = default;

    SplineEvaluator& operator=(SplineEvaluator const& x) = default;

    SplineEvaluator& operator=(SplineEvaluator&& x) = default;



    KOKKOS_FUNCTION spline_domain_type spline_domain() const noexcept
    {
        return m_spline_domain;
    }

    KOKKOS_FUNCTION bsplines_domain_type bsplines_domain() const noexcept // TODO : clarify name
    {
        return ddc::discrete_space<bsplines_type>().full_domain();
    }

    KOKKOS_FUNCTION batch_domain_type batch_domain() const noexcept
    {
        return ddc::remove_dims_of(spline_domain(), bsplines_domain());
    }

    template <class Layout, class... CoordsDims>
    KOKKOS_FUNCTION double operator()(
            ddc::Coordinate<CoordsDims...> const& coord_eval,
            ddc::ChunkSpan<double const, bsplines_domain_type, Layout, memory_space> const
                    spline_coef) const
    {
        return eval(coord_eval, spline_coef);
    }

    template <class Layout1, class Layout2, class Layout3, class... CoordsDims>
    void operator()(
            ddc::ChunkSpan<double, vals_domain_type, Layout1, memory_space> const spline_eval,
            ddc::ChunkSpan<
                    ddc::Coordinate<CoordsDims...> const,
                    vals_domain_type,
                    Layout2,
                    memory_space> const coords_eval,
            ddc::ChunkSpan<double const, spline_domain_type, Layout3, memory_space> const
                    spline_coef) const
    {
        interpolation_domain_type const interpolation_domain(spline_eval.domain());
        ddc::parallel_for_each(
                exec_space(),
                batch_domain(),
                KOKKOS_CLASS_LAMBDA(typename batch_domain_type::discrete_element_type const j) {
                    const auto spline_eval_1D = spline_eval[j];
                    const auto coords_eval_1D = coords_eval[j];
                    const auto spline_coef_1D = spline_coef[j];
                    for (auto const i : interpolation_domain) {
                        spline_eval_1D(i) = eval(coords_eval_1D(i), spline_coef_1D);
                    }
                });
    }

    template <class Layout, class... CoordsDims>
    KOKKOS_FUNCTION double deriv(
            ddc::Coordinate<CoordsDims...> const& coord_eval,
            ddc::ChunkSpan<double const, bsplines_domain_type, Layout, memory_space> const
                    spline_coef) const
    {
        return eval_no_bc<eval_deriv_type>(coord_eval, spline_coef);
    }

    template <class Layout1, class Layout2, class Layout3, class... CoordsDims>
    void deriv(
            ddc::ChunkSpan<double, vals_domain_type, Layout1, memory_space> const spline_eval,
            ddc::ChunkSpan<
                    ddc::Coordinate<CoordsDims...> const,
                    vals_domain_type,
                    Layout2,
                    memory_space> const coords_eval,
            ddc::ChunkSpan<double const, spline_domain_type, Layout3, memory_space> const
                    spline_coef) const
    {
        interpolation_domain_type const interpolation_domain(spline_eval.domain());
        ddc::parallel_for_each(
                exec_space(),
                batch_domain(),
                KOKKOS_CLASS_LAMBDA(typename batch_domain_type::discrete_element_type const j) {
                    const auto spline_eval_1D = spline_eval[j];
                    const auto coords_eval_1D = coords_eval[j];
                    const auto spline_coef_1D = spline_coef[j];
                    for (auto const i : interpolation_domain) {
                        spline_eval_1D(i)
                                = eval_no_bc<eval_deriv_type>(coords_eval_1D(i), spline_coef_1D);
                    }
                });
    }

    template <class Layout1, class Layout2>
    void integrate(
            ddc::ChunkSpan<double, batch_domain_type, Layout1, memory_space> const integrals,
            ddc::ChunkSpan<double const, spline_domain_type, Layout2, memory_space> const
                    spline_coef) const
    {
        ddc::Chunk values_alloc(
                ddc::DiscreteDomain<bsplines_type>(spline_coef.domain()),
                ddc::KokkosAllocator<double, memory_space>());
        ddc::ChunkSpan values = values_alloc.span_view();
        Kokkos::parallel_for(
                Kokkos::RangePolicy<exec_space>(0, 1),
                KOKKOS_LAMBDA(int) { ddc::discrete_space<bsplines_type>().integrals(values); });

        ddc::parallel_for_each(
                exec_space(),
                batch_domain(),
                KOKKOS_LAMBDA(typename batch_domain_type::discrete_element_type const j) {
                    integrals(j) = 0;
                    for (typename bsplines_domain_type::discrete_element_type const i :
                         values.domain()) {
                        integrals(j) += spline_coef(i, j) * values(i);
                    }
                });
    }

private:
    template <class Layout, class... CoordsDims>
    KOKKOS_INLINE_FUNCTION double eval(
            ddc::Coordinate<CoordsDims...> const& coord_eval,
            ddc::ChunkSpan<double const, bsplines_domain_type, Layout, memory_space> const
                    spline_coef) const
    {
        ddc::Coordinate<typename interpolation_mesh_type::continuous_dimension_type>
                coord_eval_interpolation
                = ddc::select<typename interpolation_mesh_type::continuous_dimension_type>(
                        coord_eval);
        if constexpr (bsplines_type::is_periodic()) {
            if (coord_eval_interpolation < ddc::discrete_space<bsplines_type>().rmin()
                || coord_eval_interpolation > ddc::discrete_space<bsplines_type>().rmax()) {
                coord_eval_interpolation -= Kokkos::floor(
                                                    (coord_eval_interpolation
                                                     - ddc::discrete_space<bsplines_type>().rmin())
                                                    / ddc::discrete_space<bsplines_type>().length())
                                            * ddc::discrete_space<bsplines_type>().length();
            }
        } else {
            if (coord_eval_interpolation < ddc::discrete_space<bsplines_type>().rmin()) {
                return m_left_extrap_rule(coord_eval_interpolation, spline_coef);
            }
            if (coord_eval_interpolation > ddc::discrete_space<bsplines_type>().rmax()) {
                return m_right_bc(coord_eval_interpolation, spline_coef);
            }
        }
        return eval_no_bc<eval_type>(coord_eval_interpolation, spline_coef);
    }

    template <class EvalType, class Layout, class... CoordsDims>
    KOKKOS_INLINE_FUNCTION double eval_no_bc(
            ddc::Coordinate<CoordsDims...> const& coord_eval,
            ddc::ChunkSpan<double const, bsplines_domain_type, Layout, memory_space> const
                    spline_coef) const
    {
        static_assert(
                std::is_same_v<EvalType, eval_type> || std::is_same_v<EvalType, eval_deriv_type>);
        ddc::DiscreteElement<bsplines_type> jmin;
        std::array<double, bsplines_type::degree() + 1> vals;
        ddc::Coordinate<typename interpolation_mesh_type::continuous_dimension_type>
                coord_eval_interpolation
                = ddc::select<typename interpolation_mesh_type::continuous_dimension_type>(
                        coord_eval);
        if constexpr (std::is_same_v<EvalType, eval_type>) {
            jmin = ddc::discrete_space<bsplines_type>().eval_basis(vals, coord_eval_interpolation);
        } else if constexpr (std::is_same_v<EvalType, eval_deriv_type>) {
            jmin = ddc::discrete_space<bsplines_type>().eval_deriv(vals, coord_eval_interpolation);
        }
        double y = 0.0;
        for (std::size_t i = 0; i < bsplines_type::degree() + 1; ++i) {
            y += spline_coef(ddc::DiscreteElement<bsplines_type>(jmin + i)) * vals[i];
        }
        return y;
    }
};
} // namespace ddc
