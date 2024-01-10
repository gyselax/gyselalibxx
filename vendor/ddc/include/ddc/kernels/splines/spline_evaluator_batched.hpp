#pragma once

#include <array>

#include <ddc/ddc.hpp>

#include "Kokkos_Macros.hpp"
#include "spline_boundary_value.hpp"
#include "spline_evaluator.hpp"
#include "view.hpp"

namespace ddc {

template <class SplineEvaluator, class... IDimX>
class SplineEvaluatorBatched
{
private:
    // Tags to determine what to evaluate
    struct eval_type
    {
    };

    struct eval_deriv_type
    {
    };

    using tag_type = typename SplineEvaluator::tag_type;

public:
    using exec_space = typename SplineEvaluator::exec_space;

    using memory_space = typename SplineEvaluator::memory_space;

    using bsplines_type = typename SplineEvaluator::bsplines_type;

    using evaluator_type = SplineEvaluator;

    using interpolation_mesh_type = typename SplineEvaluator::mesh_type;

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
    SplineEvaluator spline_evaluator;
    const spline_domain_type m_spline_domain; // Necessary ?


public:
    SplineEvaluatorBatched() = delete;

    explicit SplineEvaluatorBatched(
            spline_domain_type const& spline_domain,
            SplineBoundaryValue<bsplines_type> const& left_bc,
            SplineBoundaryValue<bsplines_type> const& right_bc)
        : spline_evaluator(left_bc, right_bc)
        , m_spline_domain(spline_domain) // Necessary ?
    {
    }

    SplineEvaluatorBatched(SplineEvaluatorBatched const& x) = default;

    SplineEvaluatorBatched(SplineEvaluatorBatched&& x) = default;

    ~SplineEvaluatorBatched() = default;

    SplineEvaluatorBatched& operator=(SplineEvaluatorBatched const& x) = default;

    SplineEvaluatorBatched& operator=(SplineEvaluatorBatched&& x) = default;



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
    double operator()(
            ddc::Coordinate<CoordsDims...> const& coord_eval,
            ddc::ChunkSpan<double const, spline_domain_type, Layout, memory_space> const
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
        ddc::for_each(
                ddc::policies::policy(exec_space()),
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
    double deriv(
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
        ddc::for_each(
                ddc::policies::policy(exec_space()),
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

        ddc::for_each(
                ddc::policies::policy(exec_space()),
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
