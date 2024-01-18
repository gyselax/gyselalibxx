#pragma once

#include <array>

#include <ddc/ddc.hpp>

#include "Kokkos_Macros.hpp"
#include "spline_boundary_value.hpp"
#include "spline_evaluator.hpp"
#include "view.hpp"

namespace ddc {

template <
        class ExecSpace,
        class MemorySpace,
        class BSplinesType1,
        class BSplinesType2,
        class interpolation_mesh_type1,
        class interpolation_mesh_type2,
        class... IDimX>
class SplineEvaluator2DBatched
{
private:
    // Tags to determine what to evaluate
    struct eval_type
    {
    };

    struct eval_deriv_type
    {
    };

public:
    using exec_space = ExecSpace;

    using memory_space = MemorySpace;

    using evaluator_type1
            = ddc::SplineEvaluator<ExecSpace, MemorySpace, BSplinesType1, interpolation_mesh_type1>;
    using evaluator_type2
            = ddc::SplineEvaluator<ExecSpace, MemorySpace, BSplinesType2, interpolation_mesh_type2>;

    using bsplines_type1 = typename evaluator_type1::bsplines_type;
    using bsplines_type2 = typename evaluator_type2::bsplines_type;

    using interpolation_domain_type1 = ddc::DiscreteDomain<interpolation_mesh_type1>;
    using interpolation_domain_type2 = ddc::DiscreteDomain<interpolation_mesh_type2>;
    using interpolation_domain_type
            = ddc::DiscreteDomain<interpolation_mesh_type1, interpolation_mesh_type2>;

    using vals_domain_type = ddc::DiscreteDomain<IDimX...>;

    using bsplines_domain_type1 = ddc::DiscreteDomain<bsplines_type1>;
    using bsplines_domain_type2 = ddc::DiscreteDomain<bsplines_type2>;
    using bsplines_domain_type = ddc::DiscreteDomain<bsplines_type1, bsplines_type2>;

    using batch_domain_type =
            typename ddc::detail::convert_type_seq_to_discrete_domain<ddc::type_seq_remove_t<
                    ddc::detail::TypeSeq<IDimX...>,
                    ddc::detail::TypeSeq<interpolation_mesh_type1, interpolation_mesh_type2>>>;

    template <typename Tag>
    using spline_dim_type = std::conditional_t<
            std::is_same_v<Tag, interpolation_mesh_type1>,
            bsplines_type1,
            std::conditional_t<std::is_same_v<Tag, interpolation_mesh_type2>, bsplines_type2, Tag>>;

    using spline_domain_type =
            typename ddc::detail::convert_type_seq_to_discrete_domain<ddc::type_seq_replace_t<
                    ddc::detail::TypeSeq<IDimX...>,
                    ddc::detail::TypeSeq<interpolation_mesh_type1, interpolation_mesh_type2>,
                    ddc::detail::TypeSeq<bsplines_type1, bsplines_type2>>>;


private:
    spline_domain_type m_spline_domain;


public:
    explicit SplineEvaluator2DBatched(
            spline_domain_type const& spline_domain,
            [[maybe_unused]] SplineBoundaryValue<bsplines_type1> const&
                    left1_bc, // Unused, to be restored in next MR
            [[maybe_unused]] SplineBoundaryValue<bsplines_type1> const& right1_bc,
            [[maybe_unused]] SplineBoundaryValue<bsplines_type2> const& left2_bc,
            [[maybe_unused]] SplineBoundaryValue<bsplines_type2> const& right2_bc)
        : m_spline_domain(spline_domain)
    {
    }

    SplineEvaluator2DBatched(SplineEvaluator2DBatched const& x) = default;

    SplineEvaluator2DBatched(SplineEvaluator2DBatched&& x) = default;

    ~SplineEvaluator2DBatched() = default;

    SplineEvaluator2DBatched& operator=(SplineEvaluator2DBatched const& x) = default;

    SplineEvaluator2DBatched& operator=(SplineEvaluator2DBatched&& x) = default;



    KOKKOS_FUNCTION spline_domain_type spline_domain() const noexcept
    {
        return m_spline_domain;
    }

    KOKKOS_FUNCTION bsplines_domain_type bsplines_domain() const noexcept // TODO : clarify name
    {
        return bsplines_domain_type(
                ddc::discrete_space<bsplines_type1>().full_domain(),
                ddc::discrete_space<bsplines_type2>().full_domain());
    }

    KOKKOS_FUNCTION batch_domain_type batch_domain() const noexcept
    {
        return ddc::remove_dims_of(spline_domain(), bsplines_domain());
    }

    template <class Layout, class... CoordsDims>
    KOKKOS_FUNCTION double operator()(
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
        interpolation_domain_type1 const interpolation_domain1(spline_eval.domain());
        interpolation_domain_type2 const interpolation_domain2(spline_eval.domain());
        ddc::for_each(
                ddc::policies::policy(exec_space()),
                batch_domain(),
                KOKKOS_CLASS_LAMBDA(typename batch_domain_type::discrete_element_type const j) {
                    const auto spline_eval_2D = spline_eval[j];
                    const auto coords_eval_2D = coords_eval[j];
                    const auto spline_coef_2D = spline_coef[j];
                    for (auto const i1 : interpolation_domain1) {
                        for (auto const i2 : interpolation_domain2) {
                            spline_eval_2D(i1, i2) = eval(coords_eval_2D(i1, i2), spline_coef_2D);
                        }
                    }
                });
    }

    template <class Layout, class... CoordsDims>
    KOKKOS_FUNCTION double deriv_dim_1(
            ddc::Coordinate<CoordsDims...> const& coord_eval,
            ddc::ChunkSpan<double const, bsplines_domain_type, Layout, memory_space> const
                    spline_coef) const
    {
        return eval_no_bc<eval_deriv_type, eval_type>(coord_eval, spline_coef);
    }

    template <class Layout, class... CoordsDims>
    KOKKOS_FUNCTION double deriv_dim_2(
            ddc::Coordinate<CoordsDims...> const& coord_eval,
            ddc::ChunkSpan<double const, bsplines_domain_type, Layout, memory_space> const
                    spline_coef) const
    {
        return eval_no_bc<eval_type, eval_deriv_type>(coord_eval, spline_coef);
    }

    template <class Layout, class... CoordsDims>
    KOKKOS_FUNCTION double deriv_1_and_2(
            ddc::Coordinate<CoordsDims...> const& coord_eval,
            ddc::ChunkSpan<double const, bsplines_domain_type, Layout, memory_space> const
                    spline_coef) const
    {
        return eval_no_bc<eval_deriv_type, eval_deriv_type>(coord_eval, spline_coef);
    }

    template <class InterestDim, class Layout, class... CoordsDims>
    KOKKOS_FUNCTION double deriv(
            ddc::Coordinate<CoordsDims...> const& coord_eval,
            ddc::ChunkSpan<double const, bsplines_domain_type, Layout, memory_space> const
                    spline_coef) const
    {
        static_assert(
                std::is_same_v<
                        InterestDim,
                        typename interpolation_mesh_type1::
                                continuous_dimension_type> || std::is_same_v<InterestDim, typename interpolation_mesh_type2::continuous_dimension_type>);
        if constexpr (std::is_same_v<
                              InterestDim,
                              typename interpolation_mesh_type1::continuous_dimension_type>) {
            return deriv_dim_1(coord_eval, spline_coef);
        } else if constexpr (std::is_same_v<
                                     InterestDim,
                                     typename interpolation_mesh_type2::
                                             continuous_dimension_type>) {
            return deriv_dim_2(coord_eval, spline_coef);
        }
    }

    template <class InterestDim1, class InterestDim2, class Layout, class... CoordsDims>
    KOKKOS_FUNCTION double deriv2(
            ddc::Coordinate<CoordsDims...> const& coord_eval,
            ddc::ChunkSpan<double const, bsplines_domain_type, Layout, memory_space> const
                    spline_coef) const
    {
        static_assert(
                (std::is_same_v<
                         InterestDim1,
                         typename interpolation_mesh_type1::
                                 continuous_dimension_type> && std::is_same_v<InterestDim2, typename interpolation_mesh_type2::continuous_dimension_type>)
                || (std::is_same_v<
                            InterestDim2,
                            typename interpolation_mesh_type1::
                                    continuous_dimension_type> && std::is_same_v<InterestDim1, typename interpolation_mesh_type2::continuous_dimension_type>));
        return deriv_1_and_2(coord_eval, spline_coef);
    }

    template <class Layout1, class Layout2, class Layout3, class... CoordsDims>
    void deriv_dim_1(
            ddc::ChunkSpan<double, vals_domain_type, Layout1, memory_space> const spline_eval,
            ddc::ChunkSpan<
                    ddc::Coordinate<CoordsDims...> const,
                    vals_domain_type,
                    Layout2,
                    memory_space> const coords_eval,
            ddc::ChunkSpan<double const, spline_domain_type, Layout3, memory_space> const
                    spline_coef) const
    {
        interpolation_domain_type1 const interpolation_domain1(spline_eval.domain());
        interpolation_domain_type2 const interpolation_domain2(spline_eval.domain());
        ddc::for_each(
                ddc::policies::policy(exec_space()),
                batch_domain(),
                KOKKOS_CLASS_LAMBDA(typename batch_domain_type::discrete_element_type const j) {
                    const auto spline_eval_2D = spline_eval[j];
                    const auto coords_eval_2D = coords_eval[j];
                    const auto spline_coef_2D = spline_coef[j];
                    for (auto const i1 : interpolation_domain1) {
                        for (auto const i2 : interpolation_domain2) {
                            spline_eval_2D(i1, i2) = eval_no_bc<
                                    eval_deriv_type,
                                    eval_type>(coords_eval_2D(i1, i2), spline_coef_2D);
                        }
                    }
                });
    }

    template <class Layout1, class Layout2, class Layout3, class... CoordsDims>
    void deriv_dim_2(
            ddc::ChunkSpan<double, vals_domain_type, Layout1, memory_space> const spline_eval,
            ddc::ChunkSpan<
                    ddc::Coordinate<CoordsDims...> const,
                    vals_domain_type,
                    Layout2,
                    memory_space> const coords_eval,
            ddc::ChunkSpan<double const, spline_domain_type, Layout3, memory_space> const
                    spline_coef) const
    {
        interpolation_domain_type1 const interpolation_domain1(spline_eval.domain());
        interpolation_domain_type2 const interpolation_domain2(spline_eval.domain());
        ddc::for_each(
                ddc::policies::policy(exec_space()),
                batch_domain(),
                KOKKOS_CLASS_LAMBDA(typename batch_domain_type::discrete_element_type const j) {
                    const auto spline_eval_2D = spline_eval[j];
                    const auto coords_eval_2D = coords_eval[j];
                    const auto spline_coef_2D = spline_coef[j];
                    for (auto const i1 : interpolation_domain1) {
                        for (auto const i2 : interpolation_domain2) {
                            spline_eval_2D(i1, i2) = eval_no_bc<
                                    eval_type,
                                    eval_deriv_type>(coords_eval_2D(i1, i2), spline_coef_2D);
                        }
                    }
                });
    }

    template <class Layout1, class Layout2, class Layout3, class... CoordsDims>
    void deriv_1_and_2(
            ddc::ChunkSpan<double, vals_domain_type, Layout1, memory_space> const spline_eval,
            ddc::ChunkSpan<
                    ddc::Coordinate<CoordsDims...> const,
                    vals_domain_type,
                    Layout2,
                    memory_space> const coords_eval,
            ddc::ChunkSpan<double const, spline_domain_type, Layout3, memory_space> const
                    spline_coef) const
    {
        interpolation_domain_type1 const interpolation_domain1(spline_eval.domain());
        interpolation_domain_type2 const interpolation_domain2(spline_eval.domain());
        ddc::for_each(
                ddc::policies::policy(exec_space()),
                batch_domain(),
                KOKKOS_CLASS_LAMBDA(typename batch_domain_type::discrete_element_type const j) {
                    const auto spline_eval_2D = spline_eval[j];
                    const auto coords_eval_2D = coords_eval[j];
                    const auto spline_coef_2D = spline_coef[j];
                    for (auto const i1 : interpolation_domain1) {
                        for (auto const i2 : interpolation_domain2) {
                            spline_eval_2D(i1, i2) = eval_no_bc<
                                    eval_deriv_type,
                                    eval_deriv_type>(coords_eval_2D(i1, i2), spline_coef_2D);
                        }
                    }
                });
    }

    template <class InterestDim, class Layout1, class Layout2, class Layout3, class... CoordsDims>
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
        static_assert(
                std::is_same_v<
                        InterestDim,
                        typename interpolation_mesh_type1::
                                continuous_dimension_type> || std::is_same_v<InterestDim, typename interpolation_mesh_type2::continuous_dimension_type>);
        if constexpr (std::is_same_v<
                              InterestDim,
                              typename interpolation_mesh_type1::continuous_dimension_type>) {
            return deriv_dim_1(spline_eval, coords_eval, spline_coef);
        } else if constexpr (std::is_same_v<
                                     InterestDim,
                                     typename interpolation_mesh_type2::
                                             continuous_dimension_type>) {
            return deriv_dim_2(spline_eval, coords_eval, spline_coef);
        }
    }

    template <
            class InterestDim1,
            class InterestDim2,
            class Layout1,
            class Layout2,
            class Layout3,
            class... CoordsDims>
    void deriv2(
            ddc::ChunkSpan<double, vals_domain_type, Layout1, memory_space> const spline_eval,
            ddc::ChunkSpan<
                    ddc::Coordinate<CoordsDims...> const,
                    vals_domain_type,
                    Layout2,
                    memory_space> const coords_eval,
            ddc::ChunkSpan<double const, spline_domain_type, Layout3, memory_space> const
                    spline_coef) const
    {
        static_assert(
                (std::is_same_v<
                         InterestDim1,
                         typename interpolation_mesh_type1::
                                 continuous_dimension_type> && std::is_same_v<InterestDim2, typename interpolation_mesh_type2::continuous_dimension_type>)
                || (std::is_same_v<
                            InterestDim2,
                            typename interpolation_mesh_type1::
                                    continuous_dimension_type> && std::is_same_v<InterestDim1, typename interpolation_mesh_type2::continuous_dimension_type>));
        return deriv_1_and_2(spline_eval, coords_eval, spline_coef);
    }

    template <class Layout1, class Layout2>
    void integrate(
            ddc::ChunkSpan<double, batch_domain_type, Layout1, memory_space> const integrals,
            ddc::ChunkSpan<double const, spline_domain_type, Layout2, memory_space> const
                    spline_coef) const
    {
        ddc::Chunk values1_alloc(
                ddc::DiscreteDomain<bsplines_type1>(spline_coef.domain()),
                ddc::KokkosAllocator<double, memory_space>());
        ddc::ChunkSpan values1 = values1_alloc.span_view();
        ddc::Chunk values2_alloc(
                ddc::DiscreteDomain<bsplines_type2>(spline_coef.domain()),
                ddc::KokkosAllocator<double, memory_space>());
        ddc::ChunkSpan values2 = values2_alloc.span_view();
        Kokkos::parallel_for(
                Kokkos::RangePolicy<exec_space>(0, 1),
                KOKKOS_LAMBDA(int) {
                    ddc::discrete_space<bsplines_type1>().integrals(values1);
                    ddc::discrete_space<bsplines_type2>().integrals(values2);
                });

        ddc::for_each(
                ddc::policies::policy(exec_space()),
                batch_domain(),
                KOKKOS_LAMBDA(typename batch_domain_type::discrete_element_type const j) {
                    integrals(j) = 0;
                    for (typename bsplines_domain_type1::discrete_element_type const i1 :
                         values1.domain()) {
                        for (typename bsplines_domain_type2::discrete_element_type const i2 :
                             values2.domain()) {
                            integrals(j) += spline_coef(i1, i2, j) * values1(i1) * values2(i2);
                        }
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
        ddc::Coordinate<typename interpolation_mesh_type1::continuous_dimension_type>
                coord_eval_interpolation1(coord_eval);
        ddc::Coordinate<typename interpolation_mesh_type2::continuous_dimension_type>
                coord_eval_interpolation2(coord_eval);
        if constexpr (bsplines_type1::is_periodic()) {
            if (coord_eval_interpolation1 < ddc::discrete_space<bsplines_type1>().rmin()
                || coord_eval_interpolation1 > ddc::discrete_space<bsplines_type1>().rmax()) {
                coord_eval_interpolation1
                        -= Kokkos::floor(
                                   (coord_eval_interpolation1
                                    - ddc::discrete_space<bsplines_type1>().rmin())
                                   / ddc::discrete_space<bsplines_type1>().length())
                           * ddc::discrete_space<bsplines_type1>().length();
            }
        }
        if constexpr (bsplines_type2::is_periodic()) {
            if (coord_eval_interpolation2 < ddc::discrete_space<bsplines_type2>().rmin()
                || coord_eval_interpolation2 > ddc::discrete_space<bsplines_type2>().rmax()) {
                coord_eval_interpolation2
                        -= Kokkos::floor(
                                   (coord_eval_interpolation2
                                    - ddc::discrete_space<bsplines_type2>().rmin())
                                   / ddc::discrete_space<bsplines_type2>().length())
                           * ddc::discrete_space<bsplines_type2>().length();
            }
        }
        return eval_no_bc<eval_type, eval_type>(
                ddc::Coordinate<
                        typename interpolation_mesh_type1::continuous_dimension_type,
                        typename interpolation_mesh_type2::continuous_dimension_type>(
                        coord_eval_interpolation1,
                        coord_eval_interpolation2),
                spline_coef);
    }

    template <class EvalType1, class EvalType2, class Layout, class... CoordsDims>
    KOKKOS_INLINE_FUNCTION double eval_no_bc(
            ddc::Coordinate<CoordsDims...> const& coord_eval,
            ddc::ChunkSpan<double const, bsplines_domain_type, Layout, memory_space> const
                    spline_coef) const
    {
        static_assert(
                std::is_same_v<EvalType1, eval_type> || std::is_same_v<EvalType1, eval_deriv_type>);
        static_assert(
                std::is_same_v<EvalType2, eval_type> || std::is_same_v<EvalType2, eval_deriv_type>);
        ddc::DiscreteElement<bsplines_type1> jmin1;
        ddc::DiscreteElement<bsplines_type2> jmin2;

        std::array<double, bsplines_type1::degree() + 1> vals1;
        std::array<double, bsplines_type2::degree() + 1> vals2;
        ddc::Coordinate<typename interpolation_mesh_type1::continuous_dimension_type>
                coord_eval_interpolation1
                = ddc::select<typename interpolation_mesh_type1::continuous_dimension_type>(
                        coord_eval);
        ddc::Coordinate<typename interpolation_mesh_type2::continuous_dimension_type>
                coord_eval_interpolation2
                = ddc::select<typename interpolation_mesh_type2::continuous_dimension_type>(
                        coord_eval);

        if constexpr (std::is_same_v<EvalType1, eval_type>) {
            jmin1 = ddc::discrete_space<bsplines_type1>()
                            .eval_basis(vals1, coord_eval_interpolation1);
        } else if constexpr (std::is_same_v<EvalType1, eval_deriv_type>) {
            jmin1 = ddc::discrete_space<bsplines_type1>()
                            .eval_deriv(vals1, coord_eval_interpolation1);
        }
        if constexpr (std::is_same_v<EvalType2, eval_type>) {
            jmin2 = ddc::discrete_space<bsplines_type2>()
                            .eval_basis(vals2, coord_eval_interpolation2);
        } else if constexpr (std::is_same_v<EvalType2, eval_deriv_type>) {
            jmin2 = ddc::discrete_space<bsplines_type2>()
                            .eval_deriv(vals2, coord_eval_interpolation2);
        }

        double y = 0.0;
        for (std::size_t i = 0; i < bsplines_type1::degree() + 1; ++i) {
            for (std::size_t j = 0; j < bsplines_type2::degree() + 1; ++j) {
                y += spline_coef(ddc::DiscreteElement<
                                 bsplines_type1,
                                 bsplines_type2>(jmin1 + i, jmin2 + j))
                     * vals1[i] * vals2[j];
            }
        }
        return y;
    }
};
} // namespace ddc
