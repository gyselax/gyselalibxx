#pragma once
#include "ddc/chunk_span.hpp"
#include "ddc/deepcopy.hpp"
#include "ddc/discrete_domain.hpp"
#include "ddc/kokkos_allocator.hpp"

#include "deriv.hpp"
#include "spline_builder.hpp"

namespace ddc {
template <class SplineBuilder, class... IDimX>
class SplineBuilderBatched
{
private:
    using tag_type = typename SplineBuilder::bsplines_type::tag_type;

public:
    using exec_space = typename SplineBuilder::exec_space;

    using memory_space = typename SplineBuilder::memory_space;

    using bsplines_type = typename SplineBuilder::bsplines_type;

    using deriv_type = ddc::Deriv<tag_type>;

    using builder_type = SplineBuilder;

    using interpolation_mesh_type = typename SplineBuilder::mesh_type;

    using interpolation_domain_type = ddc::DiscreteDomain<interpolation_mesh_type>;

    using vals_domain_type = ddc::DiscreteDomain<IDimX...>;

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

    using spline_tr_domain_type =
            typename ddc::detail::convert_type_seq_to_discrete_domain<ddc::type_seq_merge_t<
                    ddc::detail::TypeSeq<bsplines_type>,
                    ddc::type_seq_remove_t<
                            ddc::detail::TypeSeq<IDimX...>,
                            ddc::detail::TypeSeq<interpolation_mesh_type>>>>;

    using derivs_domain_type =
            typename ddc::detail::convert_type_seq_to_discrete_domain<ddc::type_seq_replace_t<
                    ddc::detail::TypeSeq<IDimX...>,
                    ddc::detail::TypeSeq<interpolation_mesh_type>,
                    ddc::detail::TypeSeq<deriv_type>>>;

    static constexpr ddc::BoundCond BcXmin = SplineBuilder::s_bc_xmin;
    static constexpr ddc::BoundCond BcXmax = SplineBuilder::s_bc_xmax;

    static constexpr int s_nbc_xmin = builder_type::s_nbc_xmin;
    static constexpr int s_nbc_xmax = builder_type::s_nbc_xmax;

private:
    builder_type spline_builder;
    const vals_domain_type m_vals_domain;

public:
    explicit SplineBuilderBatched(
            vals_domain_type const& vals_domain,
            std::optional<int> cols_per_chunk = std::nullopt,
            std::optional<unsigned int> preconditionner_max_block_size = std::nullopt)
        : spline_builder(
                ddc::select<interpolation_mesh_type>(vals_domain),
                cols_per_chunk,
                preconditionner_max_block_size)
        , m_vals_domain(vals_domain)
    {
        static_assert(
                ((BcXmin == BoundCond::PERIODIC) == (BcXmax == BoundCond::PERIODIC)),
                "Incompatible boundary conditions");
    }

    SplineBuilderBatched(SplineBuilderBatched const& x) = delete;

    SplineBuilderBatched(SplineBuilderBatched&& x) = default;

    ~SplineBuilderBatched() = default;

    SplineBuilderBatched& operator=(SplineBuilderBatched const& x) = delete;

    SplineBuilderBatched& operator=(SplineBuilderBatched&& x) = default;

    vals_domain_type vals_domain() const noexcept
    {
        return m_vals_domain;
    }

    interpolation_domain_type interpolation_domain() const noexcept
    {
        return spline_builder.interpolation_domain();
    }

    batch_domain_type batch_domain() const noexcept
    {
        return ddc::remove_dims_of(vals_domain(), interpolation_domain());
    }

    ddc::DiscreteDomain<bsplines_type> bsplines_domain() const noexcept // TODO : clarify name
    {
        return ddc::discrete_space<bsplines_type>().full_domain();
    }

    spline_domain_type spline_domain() const noexcept
    {
        return ddc::replace_dim_of<
                interpolation_mesh_type,
                bsplines_type>(vals_domain(), bsplines_domain());
    }

    spline_tr_domain_type spline_tr_domain() const noexcept
    {
        return spline_tr_domain_type(bsplines_domain(), batch_domain());
    }

    int offset() const noexcept
    {
        return spline_builder.offset();
    }

    template <class Layout>
    void operator()(
            ddc::ChunkSpan<double, spline_domain_type, Layout, memory_space> spline,
            ddc::ChunkSpan<double const, vals_domain_type, Layout, memory_space> vals,
            std::optional<
                    ddc::ChunkSpan<double const, derivs_domain_type, Layout, memory_space>> const
                    derivs_xmin
            = std::nullopt,
            std::optional<
                    ddc::ChunkSpan<double const, derivs_domain_type, Layout, memory_space>> const
                    derivs_xmax
            = std::nullopt) const;
};

template <class SplineBuilder, class... IDimX>
template <class Layout>
void SplineBuilderBatched<SplineBuilder, IDimX...>::operator()(
        ddc::ChunkSpan<double, spline_domain_type, Layout, memory_space> spline,
        ddc::ChunkSpan<double const, vals_domain_type, Layout, memory_space> vals,
        std::optional<ddc::ChunkSpan<double const, derivs_domain_type, Layout, memory_space>> const
                derivs_xmin,
        std::optional<ddc::ChunkSpan<double const, derivs_domain_type, Layout, memory_space>> const
                derivs_xmax) const
{
    assert(vals.template extent<interpolation_mesh_type>()
           == ddc::discrete_space<bsplines_type>().nbasis() - spline_builder.s_nbe_xmin
                      - spline_builder.s_nbe_xmax);

    const bool odd = spline_builder.s_odd;

    const std::size_t nbc_xmin = spline_builder.s_nbc_xmin;
    const std::size_t nbc_xmax = spline_builder.s_nbc_xmax;

    assert((BcXmin == ddc::BoundCond::HERMITE)
           != (!derivs_xmin.has_value() || derivs_xmin->template extent<deriv_type>() == 0));
    assert((BcXmax == ddc::BoundCond::HERMITE)
           != (!derivs_xmax.has_value() || derivs_xmax->template extent<deriv_type>() == 0));
    if constexpr (BcXmin == BoundCond::HERMITE) {
        assert(ddc::DiscreteElement<deriv_type>(derivs_xmin->domain().front()).uid() == 1);
    }
    if constexpr (BcXmax == BoundCond::HERMITE) {
        assert(ddc::DiscreteElement<deriv_type>(derivs_xmax->domain().front()).uid() == 1);
    }

    // Hermite boundary conditions at xmin, if any
    // NOTE: For consistency with the linear system, the i-th derivative
    //       provided by the user must be multiplied by dx^i
    if constexpr (BcXmin == BoundCond::HERMITE) {
        assert(derivs_xmin->template extent<deriv_type>() == nbc_xmin);
        auto derivs_xmin_values = *derivs_xmin;
        auto const dx_proxy = spline_builder.dx();
        ddc::for_each(
                ddc::policies::policy(exec_space()),
                batch_domain(),
                KOKKOS_LAMBDA(typename batch_domain_type::discrete_element_type j) {
                    for (int i = nbc_xmin; i > 0; --i) {
                        spline(ddc::DiscreteElement<bsplines_type>(nbc_xmin - i), j)
                                = derivs_xmin_values(ddc::DiscreteElement<deriv_type>(i), j)
                                  * ddc::detail::ipow(dx_proxy, i + odd - 1);
                    }
                });
    }

    // TODO : Consider optimizing
    // Fill spline with vals (to work in spline afterward and preserve vals)
    auto const& offset_proxy = spline_builder.offset();
    auto const& interp_size_proxy = interpolation_domain().extents();
    auto const& nbasis_proxy = ddc::discrete_space<bsplines_type>().nbasis();
    ddc::for_each(
            ddc::policies::policy(exec_space()),
            batch_domain(),
            KOKKOS_LAMBDA(typename batch_domain_type::discrete_element_type j) {
                for (int i = nbc_xmin; i < nbc_xmin + offset_proxy; ++i) {
                    spline(ddc::DiscreteElement<bsplines_type>(i), j) = 0.0;
                }
                for (int i = 0; i < interp_size_proxy; ++i) {
                    spline(ddc::DiscreteElement<bsplines_type>(nbc_xmin + i + offset_proxy), j)
                            = vals(ddc::DiscreteElement<interpolation_mesh_type>(i), j);
                }
            });
    // Hermite boundary conditions at xmax, if any
    // NOTE: For consistency with the linear system, the i-th derivative
    //       provided by the user must be multiplied by dx^i
    if constexpr (BcXmax == BoundCond::HERMITE) {
        assert(derivs_xmax->template extent<deriv_type>() == nbc_xmax);
        auto derivs_xmax_values = *derivs_xmax;
        auto const dx_proxy = spline_builder.dx();
        ddc::for_each(
                ddc::policies::policy(exec_space()),
                batch_domain(),
                KOKKOS_LAMBDA(typename batch_domain_type::discrete_element_type j) {
                    for (int i = 0; i < nbc_xmax; ++i) {
                        spline(ddc::DiscreteElement<bsplines_type>(nbasis_proxy - nbc_xmax - i), j)
                                = derivs_xmax_values(ddc::DiscreteElement<deriv_type>(i + 1), j)
                                  * ddc::detail::ipow(dx_proxy, i + odd);
                    }
                });
    }

    // TODO : Consider optimizing
    // Allocate and fill a transposed version of spline in order to get dimension of interest as last dimension (optimal for GPU, necessary for Ginkgo)
    ddc::Chunk spline_tr_alloc(spline_tr_domain(), ddc::KokkosAllocator<double, memory_space>());
    ddc::ChunkSpan spline_tr = spline_tr_alloc.span_view();
    ddc::for_each(
            ddc::policies::policy(exec_space()),
            batch_domain(),
            KOKKOS_LAMBDA(typename batch_domain_type::discrete_element_type const j) {
                for (int i = 0; i < nbasis_proxy; i++) {
                    spline_tr(ddc::DiscreteElement<bsplines_type>(i), j)
                            = spline(ddc::DiscreteElement<bsplines_type>(i + offset_proxy), j);
                }
            });
    // Create a 2D Kokkos::View to manage spline_tr as a matrix
    Kokkos::View<double**, Kokkos::LayoutRight, exec_space> bcoef_section(
            spline_tr.data_handle(),
            ddc::discrete_space<bsplines_type>().nbasis(),
            batch_domain().size());
    // Compute spline coef
    spline_builder.matrix->solve_batch_inplace(bcoef_section);
    // Transpose back spline_tr in spline
    ddc::for_each(
            ddc::policies::policy(exec_space()),
            batch_domain(),
            KOKKOS_LAMBDA(typename batch_domain_type::discrete_element_type const j) {
                for (int i = 0; i < nbasis_proxy; i++) {
                    spline(ddc::DiscreteElement<bsplines_type>(i + offset_proxy), j)
                            = spline_tr(ddc::DiscreteElement<bsplines_type>(i), j);
                }
            });

    // Not sure yet of what this part do
    if (bsplines_type::is_periodic()) {
        ddc::for_each(
                ddc::policies::policy(exec_space()),
                batch_domain(),
                KOKKOS_LAMBDA(typename batch_domain_type::discrete_element_type const j) {
                    if (offset_proxy != 0) {
                        for (int i = 0; i < offset_proxy; ++i) {
                            spline(ddc::DiscreteElement<bsplines_type>(i), j) = spline(
                                    ddc::DiscreteElement<bsplines_type>(nbasis_proxy + i),
                                    j);
                        }
                        for (std::size_t i = offset_proxy; i < bsplines_type::degree(); ++i) {
                            spline(ddc::DiscreteElement<bsplines_type>(nbasis_proxy + i), j)
                                    = spline(ddc::DiscreteElement<bsplines_type>(i), j);
                        }
                    }
                    for (std::size_t i(0); i < bsplines_type::degree(); ++i) {
                        const ddc::DiscreteElement<bsplines_type> i_start(i);
                        const ddc::DiscreteElement<bsplines_type> i_end(nbasis_proxy + i);

                        spline(i_end, j) = spline(i_start, j);
                    }
                });
    }
}
} // namespace ddc
