#pragma once
#include "ddc/chunk_span.hpp"
#include "ddc/deepcopy.hpp"
#include "ddc/discrete_domain.hpp"
#include "ddc/kokkos_allocator.hpp"

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

    static constexpr ddc::BoundCond BcXmin = SplineBuilder::s_bc_xmin;
    static constexpr ddc::BoundCond BcXmax = SplineBuilder::s_bc_xmax;

private:
    builder_type spline_builder;
    const vals_domain_type m_vals_domain;

public:
    SplineBuilderBatched(
            vals_domain_type const& vals_domain,
            std::optional<int> cols_per_par_chunk = std::nullopt,
            std::optional<int> par_chunks_per_seq_chunk = std::nullopt,
            std::optional<unsigned int> preconditionner_max_block_size = std::nullopt)
        : spline_builder(
                ddc::select<interpolation_mesh_type>(vals_domain),
                cols_per_par_chunk,
                par_chunks_per_seq_chunk,
                preconditionner_max_block_size)
        , m_vals_domain(vals_domain)
    {
        static_assert(
                BcXmin == BoundCond::PERIODIC && BcXmax == BoundCond::PERIODIC,
                "Boundary conditions other than PERIODIC are not supported yet in "
                "SpSplineBuilderBatched");
    };

    SplineBuilderBatched(SplineBuilderBatched const& x) = delete;

    SplineBuilderBatched(SplineBuilderBatched&& x) = default;

    ~SplineBuilderBatched() = default;

    SplineBuilderBatched& operator=(SplineBuilderBatched const& x) = delete;

    SplineBuilderBatched& operator=(SplineBuilderBatched&& x) = default;

    vals_domain_type const vals_domain() const noexcept
    {
        return m_vals_domain;
    }

    interpolation_domain_type const interpolation_domain() const noexcept
    {
        return spline_builder.interpolation_domain();
    }

    batch_domain_type const batch_domain() const noexcept
    {
        return ddc::remove_dims_of(vals_domain(), interpolation_domain());
    }

    ddc::DiscreteDomain<bsplines_type> const bsplines_domain() const noexcept // TODO : clarify name
    {
        return ddc::discrete_space<bsplines_type>().full_domain();
    }

    spline_domain_type const spline_domain() const noexcept
    {
        return ddc::replace_dim_of<
                interpolation_mesh_type,
                bsplines_type>(vals_domain(), bsplines_domain());
    }

    spline_tr_domain_type const spline_tr_domain() const noexcept
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
            ddc::ChunkSpan<double, vals_domain_type, Layout, memory_space> vals) const;
};

template <class SplineBuilder, class... IDimX>
template <class Layout>
void SplineBuilderBatched<SplineBuilder, IDimX...>::operator()(
        ddc::ChunkSpan<double, spline_domain_type, Layout, memory_space> spline,
        ddc::ChunkSpan<double, vals_domain_type, Layout, memory_space> vals) const
{
    const std::size_t nbc_xmin = spline_builder.s_nbc_xmin;

    // TODO : Consider optimizing
    // Fill spline with vals (to work in spline afterward and preserve vals)
    auto const& offset_proxy = spline_builder.offset();
    auto const& interp_size_proxy = interpolation_domain().extents();
    auto const& nbasis_proxy = ddc::discrete_space<bsplines_type>().nbasis();
    ddc::for_each(
            ddc::policies::policy(exec_space()),
            batch_domain(),
            DDC_LAMBDA(typename batch_domain_type::discrete_element_type j) {
                for (int i = nbc_xmin; i < nbc_xmin + offset_proxy; ++i) {
                    spline(ddc::DiscreteElement<bsplines_type>(i), j) = 0.0;
                }
                for (int i = 0; i < interp_size_proxy; ++i) {
                    spline(ddc::DiscreteElement<bsplines_type>(nbc_xmin + i + offset_proxy), j)
                            = vals(ddc::DiscreteElement<interpolation_mesh_type>(i), j);
                }
            });

    // TODO : Consider optimizing
    // Allocate and fill a transposed version of spline in order to get dimension of interest as last dimension (optimal for GPU, necessary for Ginkgo)
    ddc::Chunk spline_tr_alloc(spline_tr_domain(), ddc::KokkosAllocator<double, memory_space>());
    ddc::ChunkSpan spline_tr = spline_tr_alloc.span_view();
    ddc::for_each(
            ddc::policies::policy(exec_space()),
            batch_domain(),
            DDC_LAMBDA(typename batch_domain_type::discrete_element_type const j) {
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
            DDC_LAMBDA(typename batch_domain_type::discrete_element_type const j) {
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
                DDC_LAMBDA(typename batch_domain_type::discrete_element_type const j) {
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
