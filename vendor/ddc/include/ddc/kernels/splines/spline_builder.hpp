#pragma once
#include "ddc/chunk_span.hpp"
#include "ddc/discrete_domain.hpp"
#include "ddc/kokkos_allocator.hpp"

#include "deriv.hpp"

namespace ddc {
enum class SplineSolver {
    GINKGO
}; // Only GINKGO available atm, other solvers will be implemented in the futur

constexpr bool is_spline_interpolation_mesh_uniform(
        bool const is_uniform,
        ddc::BoundCond const BcXmin,
        ddc::BoundCond const BcXmax,
        int degree)
{
    int N_BE_MIN = n_boundary_equations(BcXmin, degree);
    int N_BE_MAX = n_boundary_equations(BcXmax, degree);
    bool is_periodic = (BcXmin == ddc::BoundCond::PERIODIC) && (BcXmax == ddc::BoundCond::PERIODIC);
    return is_uniform && ((N_BE_MIN != 0 && N_BE_MAX != 0) || is_periodic);
}

template <
        class ExecSpace,
        class MemorySpace,
        class BSplines,
        class InterpolationMesh,
        ddc::BoundCond BcXmin,
        ddc::BoundCond BcXmax,
        SplineSolver Solver,
        class... IDimX>
class SplineBuilder
{
    static_assert(
            (BSplines::is_periodic() && (BcXmin == ddc::BoundCond::PERIODIC)
             && (BcXmax == ddc::BoundCond::PERIODIC))
            || (!BSplines::is_periodic() && (BcXmin != ddc::BoundCond::PERIODIC)
                && (BcXmax != ddc::BoundCond::PERIODIC)));
    static_assert(!BSplines::is_radial());

private:
    using tag_type = typename InterpolationMesh::continuous_dimension_type;

public:
    using exec_space = ExecSpace;

    using memory_space = MemorySpace;

    using interpolation_mesh_type = InterpolationMesh;

    using bsplines_type = BSplines;

    using deriv_type = ddc::Deriv<tag_type>;

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

    static constexpr bool s_odd = BSplines::degree() % 2;
    static constexpr int s_nbe_xmin = n_boundary_equations(BcXmin, BSplines::degree());
    static constexpr int s_nbe_xmax = n_boundary_equations(BcXmax, BSplines::degree());
    static constexpr int s_nbc_xmin = n_user_input(BcXmin, BSplines::degree());
    static constexpr int s_nbc_xmax = n_user_input(BcXmax, BSplines::degree());

    static constexpr ddc::BoundCond s_bc_xmin = BcXmin;
    static constexpr ddc::BoundCond s_bc_xmax = BcXmax;

private:
    vals_domain_type m_vals_domain;

    int m_offset;

    double m_dx; // average cell size for normalization of derivatives

    // interpolator specific
    std::unique_ptr<ddc::detail::Matrix> matrix;

public:
    int compute_offset(interpolation_domain_type const& interpolation_domain);

    explicit SplineBuilder(
            vals_domain_type const& vals_domain,
            std::optional<int> cols_per_chunk = std::nullopt,
            std::optional<unsigned int> preconditionner_max_block_size = std::nullopt)
        : m_vals_domain(vals_domain)
        , m_offset(compute_offset(interpolation_domain()))
        , m_dx((ddc::discrete_space<BSplines>().rmax() - ddc::discrete_space<BSplines>().rmin())
               / ddc::discrete_space<BSplines>().ncells())
    {
        static_assert(
                ((BcXmin == BoundCond::PERIODIC) == (BcXmax == BoundCond::PERIODIC)),
                "Incompatible boundary conditions");

        // Calculate block sizes
        int lower_block_size, upper_block_size;
        if constexpr (bsplines_type::is_uniform()) {
            compute_block_sizes_uniform(lower_block_size, upper_block_size);
        } else {
            compute_block_sizes_non_uniform(lower_block_size, upper_block_size);
        }
        allocate_matrix(
                lower_block_size,
                upper_block_size,
                cols_per_chunk,
                preconditionner_max_block_size);
    }

    SplineBuilder(SplineBuilder const& x) = delete;

    SplineBuilder(SplineBuilder&& x) = default;

    ~SplineBuilder() = default;

    SplineBuilder& operator=(SplineBuilder const& x) = delete;

    SplineBuilder& operator=(SplineBuilder&& x) = default;

    vals_domain_type vals_domain() const noexcept
    {
        return m_vals_domain;
    }

    interpolation_domain_type interpolation_domain() const noexcept
    {
        return interpolation_domain_type(vals_domain());
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

    derivs_domain_type derivs_xmin_domain() const noexcept
    {
        return ddc::replace_dim_of<interpolation_mesh_type, deriv_type>(
                vals_domain(),
                ddc::DiscreteDomain<deriv_type>(
                        ddc::DiscreteElement<deriv_type>(1),
                        ddc::DiscreteVector<deriv_type>(s_nbc_xmin)));
    }

    derivs_domain_type derivs_xmax_domain() const noexcept
    {
        return ddc::replace_dim_of<interpolation_mesh_type, deriv_type>(
                vals_domain(),
                ddc::DiscreteDomain<deriv_type>(
                        ddc::DiscreteElement<deriv_type>(1),
                        ddc::DiscreteVector<deriv_type>(s_nbc_xmax)));
    }

    /**
     * @brief Get the interpolation matrix.
     *
     * Get the interpolation matrix. This can be useful for debugging (as it allows
     * one to print the matrix) or for more complex quadrature schemes.
     *
     * @return A reference to the interpolation matrix.
     */
    const ddc::detail::Matrix& get_interpolation_matrix() const noexcept
    {
        return *matrix;
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

private:
    void compute_block_sizes_uniform(int& lower_block_size, int& upper_block_size) const;

    void compute_block_sizes_non_uniform(int& lower_block_size, int& upper_block_size) const;

    void allocate_matrix(
            int lower_block_size,
            int upper_block_size,
            std::optional<int> cols_per_chunk = std::nullopt,
            std::optional<unsigned int> preconditionner_max_block_size = std::nullopt);

    void build_matrix_system();
};

template <
        class ExecSpace,
        class MemorySpace,
        class BSplines,
        class InterpolationMesh,
        ddc::BoundCond BcXmin,
        ddc::BoundCond BcXmax,
        SplineSolver Solver,
        class... IDimX>
int SplineBuilder<
        ExecSpace,
        MemorySpace,
        BSplines,
        InterpolationMesh,
        BcXmin,
        BcXmax,
        Solver,
        IDimX...>::compute_offset(interpolation_domain_type const& interpolation_domain)
{
    int offset;
    if constexpr (bsplines_type::is_periodic()) {
        // Calculate offset so that the matrix is diagonally dominant
        std::array<double, bsplines_type::degree() + 1> values;
        ddc::DiscreteElement<interpolation_mesh_type> start(interpolation_domain.front());
        auto jmin = ddc::discrete_space<BSplines>()
                            .eval_basis(values, ddc::coordinate(start + BSplines::degree()));
        if constexpr (bsplines_type::degree() % 2 == 0) {
            offset = jmin.uid() - start.uid() + bsplines_type::degree() / 2 - BSplines::degree();
        } else {
            int const mid = bsplines_type::degree() / 2;
            offset = jmin.uid() - start.uid() + (values[mid] > values[mid + 1] ? mid : mid + 1)
                     - BSplines::degree();
        }
    } else {
        offset = 0;
    }
    return offset;
}

template <
        class ExecSpace,
        class MemorySpace,
        class BSplines,
        class InterpolationMesh,
        ddc::BoundCond BcXmin,
        ddc::BoundCond BcXmax,
        SplineSolver Solver,
        class... IDimX>
void SplineBuilder<
        ExecSpace,
        MemorySpace,
        BSplines,
        InterpolationMesh,
        BcXmin,
        BcXmax,
        Solver,
        IDimX...>::compute_block_sizes_uniform(int& lower_block_size, int& upper_block_size) const
{
    switch (BcXmin) {
    case ddc::BoundCond::PERIODIC:
        upper_block_size = (bsplines_type::degree()) / 2;
        break;
    case ddc::BoundCond::NATURAL:
    case ddc::BoundCond::HERMITE:
        upper_block_size = s_nbc_xmin;
        break;
    case ddc::BoundCond::GREVILLE:
        upper_block_size = bsplines_type::degree() - 1;
        break;
    default:
        throw std::runtime_error("ddc::BoundCond not handled");
    }
    switch (BcXmax) {
    case ddc::BoundCond::PERIODIC:
        lower_block_size = (bsplines_type::degree()) / 2;
        break;
    case ddc::BoundCond::NATURAL:
    case ddc::BoundCond::HERMITE:
        lower_block_size = s_nbc_xmax;
        break;
    case ddc::BoundCond::GREVILLE:
        lower_block_size = bsplines_type::degree() - 1;
        break;
    default:
        throw std::runtime_error("ddc::BoundCond not handled");
    }
}

template <
        class ExecSpace,
        class MemorySpace,
        class BSplines,
        class InterpolationMesh,
        ddc::BoundCond BcXmin,
        ddc::BoundCond BcXmax,
        SplineSolver Solver,
        class... IDimX>
void SplineBuilder<
        ExecSpace,
        MemorySpace,
        BSplines,
        InterpolationMesh,
        BcXmin,
        BcXmax,
        Solver,
        IDimX...>::compute_block_sizes_non_uniform(int& lower_block_size, int& upper_block_size)
        const
{
    switch (BcXmin) {
    case ddc::BoundCond::PERIODIC:
        upper_block_size = bsplines_type::degree() - 1;
        break;
    case ddc::BoundCond::HERMITE:
        upper_block_size = s_nbc_xmin + 1;
        break;
    case ddc::BoundCond::GREVILLE:
        upper_block_size = bsplines_type::degree() - 1;
        break;
    default:
        throw std::runtime_error("ddc::BoundCond not handled");
    }
    switch (BcXmax) {
    case ddc::BoundCond::PERIODIC:
        lower_block_size = bsplines_type::degree() - 1;
        break;
    case ddc::BoundCond::HERMITE:
        lower_block_size = s_nbc_xmax + 1;
        break;
    case ddc::BoundCond::GREVILLE:
        lower_block_size = bsplines_type::degree() - 1;
        break;
    default:
        throw std::runtime_error("ddc::BoundCond not handled");
    }
}

template <
        class ExecSpace,
        class MemorySpace,
        class BSplines,
        class InterpolationMesh,
        ddc::BoundCond BcXmin,
        ddc::BoundCond BcXmax,
        SplineSolver Solver,
        class... IDimX>
void SplineBuilder<
        ExecSpace,
        MemorySpace,
        BSplines,
        InterpolationMesh,
        BcXmin,
        BcXmax,
        Solver,
        IDimX...>::
        allocate_matrix(
                [[maybe_unused]] int lower_block_size,
                [[maybe_unused]] int upper_block_size,
                std::optional<int> cols_per_chunk,
                std::optional<unsigned int> preconditionner_max_block_size)
{
    // Special case: linear spline
    // No need for matrix assembly
    // (desactivated)
    /*
    if constexpr (bsplines_type::degree() == 1)
        return;
	*/

    matrix = ddc::detail::MatrixMaker::make_new_sparse<ExecSpace>(
            ddc::discrete_space<BSplines>().nbasis(),
            cols_per_chunk,
            preconditionner_max_block_size);

    build_matrix_system();

    matrix->factorize();
}

template <
        class ExecSpace,
        class MemorySpace,
        class BSplines,
        class InterpolationMesh,
        ddc::BoundCond BcXmin,
        ddc::BoundCond BcXmax,
        SplineSolver Solver,
        class... IDimX>
void SplineBuilder<
        ExecSpace,
        MemorySpace,
        BSplines,
        InterpolationMesh,
        BcXmin,
        BcXmax,
        Solver,
        IDimX...>::build_matrix_system()
{
    // Hermite boundary conditions at xmin, if any
    if constexpr (BcXmin == ddc::BoundCond::HERMITE) {
        std::array<double, (bsplines_type::degree() / 2 + 1) * (bsplines_type::degree() + 1)>
                derivs_ptr;
        ddc::DSpan2D
                derivs(derivs_ptr.data(),
                       bsplines_type::degree() + 1,
                       bsplines_type::degree() / 2 + 1);
        ddc::discrete_space<BSplines>().eval_basis_and_n_derivs(
                derivs,
                ddc::discrete_space<BSplines>().rmin(),
                s_nbc_xmin);

        // In order to improve the condition number of the matrix, we normalize
        // all derivatives by multiplying the i-th derivative by dx^i
        for (std::size_t i = 0; i < bsplines_type::degree() + 1; ++i) {
            for (std::size_t j = 1; j < bsplines_type::degree() / 2 + 1; ++j) {
                derivs(i, j) *= ddc::detail::ipow(m_dx, j);
            }
        }

        // iterate only to deg as last bspline is 0
        for (std::size_t i = 0; i < s_nbc_xmin; ++i) {
            for (std::size_t j = 0; j < bsplines_type::degree(); ++j) {
                matrix->set_element(i, j, derivs(j, s_nbc_xmin - i - 1 + s_odd));
            }
        }
    }

    // Interpolation points
    std::array<double, bsplines_type::degree() + 1> values;
    int start = interpolation_domain().front().uid();
    ddc::for_each(interpolation_domain(), [&](auto ix) {
        auto jmin = ddc::discrete_space<BSplines>().eval_basis(
                values,
                ddc::coordinate(ddc::DiscreteElement<interpolation_mesh_type>(ix)));
        for (std::size_t s = 0; s < bsplines_type::degree() + 1; ++s) {
            int const j = ddc::detail::
                    modulo(int(jmin.uid() - m_offset + s),
                           (int)ddc::discrete_space<BSplines>().nbasis());
            matrix->set_element(ix.uid() - start + s_nbc_xmin, j, values[s]);
        }
    });

    // Hermite boundary conditions at xmax, if any
    if constexpr (BcXmax == ddc::BoundCond::HERMITE) {
        std::array<double, (bsplines_type::degree() / 2 + 1) * (bsplines_type::degree() + 1)>
                derivs_ptr;
        std::experimental::mdspan<
                double,
                std::experimental::extents<
                        std::size_t,
                        bsplines_type::degree() + 1,
                        bsplines_type::degree() / 2 + 1>> const derivs(derivs_ptr.data());

        ddc::discrete_space<BSplines>().eval_basis_and_n_derivs(
                derivs,
                ddc::discrete_space<BSplines>().rmax(),
                s_nbc_xmax);

        // In order to improve the condition number of the matrix, we normalize
        // all derivatives by multiplying the i-th derivative by dx^i
        for (std::size_t i = 0; i < bsplines_type::degree() + 1; ++i) {
            for (std::size_t j = 1; j < bsplines_type::degree() / 2 + 1; ++j) {
                derivs(i, j) *= ddc::detail::ipow(m_dx, j);
            }
        }

        int const i0 = ddc::discrete_space<BSplines>().nbasis() - s_nbc_xmax;
        int const j0 = ddc::discrete_space<BSplines>().nbasis() - bsplines_type::degree();
        for (std::size_t j = 0; j < bsplines_type::degree(); ++j) {
            for (std::size_t i = 0; i < s_nbc_xmax; ++i) {
                matrix->set_element(i0 + i, j0 + j, derivs(j + 1, i + s_odd));
            }
        }
    }
}

template <
        class ExecSpace,
        class MemorySpace,
        class BSplines,
        class InterpolationMesh,
        ddc::BoundCond BcXmin,
        ddc::BoundCond BcXmax,
        SplineSolver Solver,
        class... IDimX>
template <class Layout>
void SplineBuilder<
        ExecSpace,
        MemorySpace,
        BSplines,
        InterpolationMesh,
        BcXmin,
        BcXmax,
        Solver,
        IDimX...>::
operator()(
        ddc::ChunkSpan<double, spline_domain_type, Layout, memory_space> spline,
        ddc::ChunkSpan<double const, vals_domain_type, Layout, memory_space> vals,
        std::optional<ddc::ChunkSpan<double const, derivs_domain_type, Layout, memory_space>> const
                derivs_xmin,
        std::optional<ddc::ChunkSpan<double const, derivs_domain_type, Layout, memory_space>> const
                derivs_xmax) const
{
    assert(vals.template extent<interpolation_mesh_type>()
           == ddc::discrete_space<bsplines_type>().nbasis() - s_nbe_xmin - s_nbe_xmax);

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
        assert(derivs_xmin->template extent<deriv_type>() == s_nbc_xmin);
        auto derivs_xmin_values = *derivs_xmin;
        auto const dx_proxy = m_dx;
        ddc::parallel_for_each(
                exec_space(),
                batch_domain(),
                KOKKOS_LAMBDA(typename batch_domain_type::discrete_element_type j) {
                    for (int i = s_nbc_xmin; i > 0; --i) {
                        spline(ddc::DiscreteElement<bsplines_type>(s_nbc_xmin - i), j)
                                = derivs_xmin_values(ddc::DiscreteElement<deriv_type>(i), j)
                                  * ddc::detail::ipow(dx_proxy, i + s_odd - 1);
                    }
                });
    }

    // TODO : Consider optimizing
    // Fill spline with vals (to work in spline afterward and preserve vals)
    auto const& offset_proxy = m_offset;
    auto const& interp_size_proxy = interpolation_domain().extents();
    auto const& nbasis_proxy = ddc::discrete_space<bsplines_type>().nbasis();
    ddc::parallel_for_each(
            exec_space(),
            batch_domain(),
            KOKKOS_LAMBDA(typename batch_domain_type::discrete_element_type j) {
                for (int i = s_nbc_xmin; i < s_nbc_xmin + offset_proxy; ++i) {
                    spline(ddc::DiscreteElement<bsplines_type>(i), j) = 0.0;
                }
                for (int i = 0; i < interp_size_proxy; ++i) {
                    spline(ddc::DiscreteElement<bsplines_type>(s_nbc_xmin + i + offset_proxy), j)
                            = vals(ddc::DiscreteElement<interpolation_mesh_type>(i), j);
                }
            });

    // Hermite boundary conditions at xmax, if any
    // NOTE: For consistency with the linear system, the i-th derivative
    //       provided by the user must be multiplied by dx^i
    if constexpr (BcXmax == BoundCond::HERMITE) {
        assert(derivs_xmax->template extent<deriv_type>() == s_nbc_xmax);
        auto derivs_xmax_values = *derivs_xmax;
        auto const dx_proxy = m_dx;
        ddc::parallel_for_each(
                exec_space(),
                batch_domain(),
                KOKKOS_LAMBDA(typename batch_domain_type::discrete_element_type j) {
                    for (int i = 0; i < s_nbc_xmax; ++i) {
                        spline(ddc::DiscreteElement<bsplines_type>(nbasis_proxy - s_nbc_xmax - i),
                               j)
                                = derivs_xmax_values(ddc::DiscreteElement<deriv_type>(i + 1), j)
                                  * ddc::detail::ipow(dx_proxy, i + s_odd);
                    }
                });
    }

    // TODO : Consider optimizing
    // Allocate and fill a transposed version of spline in order to get dimension of interest as last dimension (optimal for GPU, necessary for Ginkgo)
    ddc::Chunk spline_tr_alloc(spline_tr_domain(), ddc::KokkosAllocator<double, memory_space>());
    ddc::ChunkSpan spline_tr = spline_tr_alloc.span_view();
    ddc::parallel_for_each(
            exec_space(),
            batch_domain(),
            KOKKOS_LAMBDA(typename batch_domain_type::discrete_element_type const j) {
                for (std::size_t i = 0; i < nbasis_proxy; i++) {
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
    matrix->solve_batch_inplace(bcoef_section);
    // Transpose back spline_tr in spline
    ddc::parallel_for_each(
            exec_space(),
            batch_domain(),
            KOKKOS_LAMBDA(typename batch_domain_type::discrete_element_type const j) {
                for (std::size_t i = 0; i < nbasis_proxy; i++) {
                    spline(ddc::DiscreteElement<bsplines_type>(i + offset_proxy), j)
                            = spline_tr(ddc::DiscreteElement<bsplines_type>(i), j);
                }
            });

    // Not sure yet of what this part do
    if (bsplines_type::is_periodic()) {
        ddc::parallel_for_each(
                exec_space(),
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
