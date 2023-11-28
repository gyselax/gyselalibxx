#pragma once

#include <cstdint>
#include <iostream>
#include <memory>
#include <optional>
#include <stdexcept>
#include <vector>

#include <ddc/ddc.hpp>

#include "math_tools.hpp"
#include "matrix.hpp"
#include "matrix_maker.hpp"
#include "spline_boundary_conditions.hpp"
#include "view.hpp"

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
        class interpolation_mesh_type,
        ddc::BoundCond BcXmin,
        ddc::BoundCond BcXmax,
        SplineSolver Solver = ddc::SplineSolver::GINKGO>
class SplineBuilder
{
    static_assert(
            (BSplines::is_periodic() && (BcXmin == ddc::BoundCond::PERIODIC)
             && (BcXmax == ddc::BoundCond::PERIODIC))
            || (!BSplines::is_periodic() && (BcXmin != ddc::BoundCond::PERIODIC)
                && (BcXmax != ddc::BoundCond::PERIODIC)));
    static_assert(!BSplines::is_radial());

private:
    using tag_type = typename interpolation_mesh_type::continuous_dimension_type;

public:
    using exec_space = ExecSpace;

    using memory_space = MemorySpace;

    using bsplines_type = BSplines;

    using mesh_type = interpolation_mesh_type;

    using interpolation_domain_type = ddc::DiscreteDomain<interpolation_mesh_type>;

public:
    static constexpr bool s_odd = BSplines::degree() % 2;

    static constexpr int s_nbe_xmin = n_boundary_equations(BcXmin, BSplines::degree());

    static constexpr int s_nbe_xmax = n_boundary_equations(BcXmax, BSplines::degree());

    static constexpr int s_nbc_xmin = n_user_input(BcXmin, BSplines::degree());

    static constexpr int s_nbc_xmax = n_user_input(BcXmax, BSplines::degree());

    static constexpr ddc::BoundCond s_bc_xmin = BcXmin;
    static constexpr ddc::BoundCond s_bc_xmax = BcXmax;

    // interpolator specific
    std::unique_ptr<ddc::detail::Matrix> matrix;

private:
    const int m_offset;

    interpolation_domain_type m_interpolation_domain;

    double m_dx; // average cell size for normalization of derivatives

public:
    int compute_offset(interpolation_domain_type const& interpolation_domain);

    SplineBuilder(
            interpolation_domain_type const& interpolation_domain,
            std::optional<int> cols_per_par_chunk = std::nullopt,
            std::optional<int> par_chunks_per_seq_chunk = std::nullopt,
            std::optional<unsigned int> preconditionner_max_block_size = std::nullopt)
        : matrix(nullptr)
        , m_offset(compute_offset(interpolation_domain))
        , m_interpolation_domain(interpolation_domain)
        , m_dx((ddc::discrete_space<BSplines>().rmax() - ddc::discrete_space<BSplines>().rmin())
               / ddc::discrete_space<BSplines>().ncells())
    {
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
                cols_per_par_chunk,
                par_chunks_per_seq_chunk,
                preconditionner_max_block_size);
    }

    SplineBuilder(SplineBuilder const& x) = delete;

    SplineBuilder(SplineBuilder&& x) = default;

    ~SplineBuilder() = default;

    SplineBuilder& operator=(SplineBuilder const& x) = delete;

    SplineBuilder& operator=(SplineBuilder&& x) = default;

    template <class Layout>
    void operator()(
            ddc::ChunkSpan<double, ddc::DiscreteDomain<bsplines_type>, Layout, MemorySpace> spline,
            ddc::ChunkSpan<double, interpolation_domain_type, Layout, MemorySpace> vals,
            std::optional<ddc::CDSpan1D> const derivs_xmin = std::nullopt,
            std::optional<ddc::CDSpan1D> const derivs_xmax = std::nullopt) const;

    interpolation_domain_type const& interpolation_domain() const noexcept
    {
        return m_interpolation_domain;
    }

    int offset() const noexcept
    {
        return m_offset;
    }

    ddc::DiscreteDomain<BSplines> spline_domain() const noexcept
    {
        return ddc::discrete_space<BSplines>().full_domain();
    }

    template <class Layout>
    void compute_interpolant_degree1( // Seems to need to be public for GPU compiling
            ddc::ChunkSpan<double, ddc::DiscreteDomain<bsplines_type>, Layout, MemorySpace> spline,
            ddc::ChunkSpan<double, interpolation_domain_type, Layout, MemorySpace> vals) const;


private:
    void compute_block_sizes_uniform(int& lower_block_size, int& upper_block_size) const;

    void compute_block_sizes_non_uniform(int& lower_block_size, int& upper_block_size) const;

    void allocate_matrix(
            int lower_block_size,
            int upper_block_size,
            std::optional<int> cols_per_par_chunk = std::nullopt,
            std::optional<int> par_chunks_per_seq_chunk = std::nullopt,
            std::optional<unsigned int> preconditionner_max_block_size = std::nullopt);

    void build_matrix_system();
};

template <
        class ExecSpace,
        class MemorySpace,
        class BSplines,
        class interpolation_mesh_type,
        ddc::BoundCond BcXmin,
        ddc::BoundCond BcXmax,
        SplineSolver Solver>
int SplineBuilder<
        ExecSpace,
        MemorySpace,
        BSplines,
        interpolation_mesh_type,
        BcXmin,
        BcXmax,
        Solver>::compute_offset(interpolation_domain_type const& interpolation_domain)
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

//-------------------------------------------------------------------------------------------------
/************************************************************************************
 *                         Compute interpolant functions *
 ************************************************************************************/

template <
        class ExecSpace,
        class MemorySpace,
        class BSplines,
        class interpolation_mesh_type,
        ddc::BoundCond BcXmin,
        ddc::BoundCond BcXmax,
        SplineSolver Solver>
template <class Layout>
void SplineBuilder<
        ExecSpace,
        MemorySpace,
        BSplines,
        interpolation_mesh_type,
        BcXmin,
        BcXmax,
        Solver>::
        compute_interpolant_degree1(
                ddc::ChunkSpan<double, ddc::DiscreteDomain<bsplines_type>, Layout, MemorySpace>
                        spline,
                ddc::ChunkSpan<double, interpolation_domain_type, Layout, MemorySpace> vals) const
{
    auto const& nbasis_proxy = ddc::discrete_space<bsplines_type>().nbasis();
    Kokkos::parallel_for(
            Kokkos::RangePolicy<exec_space>(0, 1),
            KOKKOS_LAMBDA(int) {
                for (std::size_t i = 0; i < nbasis_proxy; ++i) {
                    spline(ddc::DiscreteElement<bsplines_type>(i))
                            = vals(ddc::DiscreteElement<interpolation_mesh_type>(i));
                }
            });
    if constexpr (bsplines_type::is_periodic()) {
        Kokkos::parallel_for(
                Kokkos::RangePolicy<exec_space>(0, 1),
                KOKKOS_LAMBDA(int) {
                    spline(ddc::DiscreteElement<bsplines_type>(nbasis_proxy))
                            = spline(ddc::DiscreteElement<bsplines_type>(0));
                });
    }
}

//-------------------------------------------------------------------------------------------------

template <
        class ExecSpace,
        class MemorySpace,
        class BSplines,
        class interpolation_mesh_type,
        ddc::BoundCond BcXmin,
        ddc::BoundCond BcXmax,
        SplineSolver Solver>
template <class Layout>
void SplineBuilder<
        ExecSpace,
        MemorySpace,
        BSplines,
        interpolation_mesh_type,
        BcXmin,
        BcXmax,
        Solver>::
operator()(
        ddc::ChunkSpan<double, ddc::DiscreteDomain<bsplines_type>, Layout, MemorySpace> spline,
        ddc::ChunkSpan<double, interpolation_domain_type, Layout, MemorySpace> vals,
        [[maybe_unused]] std::optional<ddc::CDSpan1D> const derivs_xmin,
        [[maybe_unused]] std::optional<ddc::CDSpan1D> const derivs_xmax) const
{
    assert(vals.template extent<interpolation_mesh_type>()
           == ddc::discrete_space<BSplines>().nbasis() - s_nbe_xmin - s_nbe_xmax);
    //    assert(spline.belongs_to_space(ddc::discrete_space<BSplines>()));
    // TODO: LOG Errors
    if constexpr (bsplines_type::degree() == 1)
        return compute_interpolant_degree1(spline, vals);

    assert((BcXmin == ddc::BoundCond::HERMITE)
           != (!derivs_xmin.has_value() || derivs_xmin->extent(0) == 0));
    assert((BcXmax == ddc::BoundCond::HERMITE)
           != (!derivs_xmax.has_value() || derivs_xmax->extent(0) == 0));

    // Hermite boundary conditions at xmin, if any
    // NOTE: For consistency with the linear system, the i-th derivative
    //       provided by the user must be multiplied by dx^i
    if constexpr (BcXmin == ddc::BoundCond::HERMITE) {
        assert(derivs_xmin->extent(0) == s_nbc_xmin);
        for (int i = s_nbc_xmin; i > 0; --i) {
            spline(ddc::DiscreteElement<bsplines_type>(s_nbc_xmin - i))
                    = (*derivs_xmin)(i - 1) * ddc::detail::ipow(m_dx, i + s_odd - 1);
        }
    }
    auto const& offset_proxy = offset();
    auto const& interp_size_proxy = interpolation_domain().extents();
    auto const& nbasis_proxy = ddc::discrete_space<bsplines_type>().nbasis();
    Kokkos::parallel_for(
            Kokkos::RangePolicy<exec_space>(0, 1),
            KOKKOS_LAMBDA(int) {
                for (int i = s_nbc_xmin; i < s_nbc_xmin + offset_proxy; ++i) {
                    spline(ddc::DiscreteElement<bsplines_type>(i)) = 0.0;
                }
                for (int i = 0; i < interp_size_proxy; ++i) {
                    spline(ddc::DiscreteElement<bsplines_type>(s_nbc_xmin + i + offset_proxy))
                            = vals(ddc::DiscreteElement<interpolation_mesh_type>(i));
                }
            });

    // Hermite boundary conditions at xmax, if any
    // NOTE: For consistency with the linear system, the i-th derivative
    //       provided by the user must be multiplied by dx^i
    if constexpr (BcXmax == ddc::BoundCond::HERMITE) {
        assert(derivs_xmax->extent(0) == s_nbc_xmax);
        for (int i = 0; i < s_nbc_xmax; ++i) {
            spline(ddc::DiscreteElement<bsplines_type>(
                    ddc::discrete_space<BSplines>().nbasis() - s_nbc_xmax + i))
                    = (*derivs_xmax)(i)*ddc::detail::ipow(m_dx, i + s_odd);
        }
    }

    Kokkos::View<double**, Kokkos::LayoutRight, exec_space> bcoef_section(
            spline.data_handle() + m_offset,
            ddc::discrete_space<BSplines>().nbasis());
    matrix->solve_batch_inplace(bcoef_section);

    if constexpr (bsplines_type::is_periodic()) {
        Kokkos::parallel_for(
                Kokkos::RangePolicy<exec_space>(0, 1),
                KOKKOS_LAMBDA(int) {
                    if (offset_proxy != 0) {
                        for (int i = 0; i < offset_proxy; ++i) {
                            spline(ddc::DiscreteElement<bsplines_type>(i))
                                    = spline(ddc::DiscreteElement<bsplines_type>(nbasis_proxy + i));
                        }
                        for (std::size_t i = offset_proxy; i < bsplines_type::degree(); ++i) {
                            spline(ddc::DiscreteElement<bsplines_type>(nbasis_proxy + i))
                                    = spline(ddc::DiscreteElement<bsplines_type>(i));
                        }
                    }
                });
    }
}

//-------------------------------------------------------------------------------------------------
/************************************************************************************
 *                            Compute num diags functions *
 ************************************************************************************/

template <
        class ExecSpace,
        class MemorySpace,
        class BSplines,
        class interpolation_mesh_type,
        ddc::BoundCond BcXmin,
        ddc::BoundCond BcXmax,
        SplineSolver Solver>
void SplineBuilder<
        ExecSpace,
        MemorySpace,
        BSplines,
        interpolation_mesh_type,
        BcXmin,
        BcXmax,
        Solver>::compute_block_sizes_uniform(int& lower_block_size, int& upper_block_size) const
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

//-------------------------------------------------------------------------------------------------

template <
        class ExecSpace,
        class MemorySpace,
        class BSplines,
        class interpolation_mesh_type,
        ddc::BoundCond BcXmin,
        ddc::BoundCond BcXmax,
        SplineSolver Solver>
void SplineBuilder<
        ExecSpace,
        MemorySpace,
        BSplines,
        interpolation_mesh_type,
        BcXmin,
        BcXmax,
        Solver>::compute_block_sizes_non_uniform(int& lower_block_size, int& upper_block_size) const
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

//-------------------------------------------------------------------------------------------------
/************************************************************************************
 *                            Initialize matrix functions *
 ************************************************************************************/

template <
        class ExecSpace,
        class MemorySpace,
        class BSplines,
        class interpolation_mesh_type,
        ddc::BoundCond BcXmin,
        ddc::BoundCond BcXmax,
        SplineSolver Solver>
void SplineBuilder<
        ExecSpace,
        MemorySpace,
        BSplines,
        interpolation_mesh_type,
        BcXmin,
        BcXmax,
        Solver>::
        allocate_matrix(
                [[maybe_unused]] int lower_block_size,
                [[maybe_unused]] int upper_block_size,
                std::optional<int> cols_per_par_chunk,
                std::optional<int> par_chunks_per_seq_chunk,
                std::optional<unsigned int> preconditionner_max_block_size)
{
    // Special case: linear spline
    // No need for matrix assembly
    if constexpr (bsplines_type::degree() == 1)
        return;

    if constexpr (bsplines_type::is_periodic()) {
        if (Solver == SplineSolver::GINKGO) {
            matrix = ddc::detail::MatrixMaker::make_new_sparse<ExecSpace>(
                    ddc::discrete_space<BSplines>().nbasis(),
                    cols_per_par_chunk,
                    par_chunks_per_seq_chunk,
                    preconditionner_max_block_size);
        }
    } else {
        if (Solver == SplineSolver::GINKGO) {
            matrix = ddc::detail::MatrixMaker::make_new_sparse<ExecSpace>(
                    ddc::discrete_space<BSplines>().nbasis(),
                    cols_per_par_chunk,
                    par_chunks_per_seq_chunk,
                    preconditionner_max_block_size);
        }
    }

    build_matrix_system();

    matrix->factorize();
}

//-------------------------------------------------------------------------------------------------

template <
        class ExecSpace,
        class MemorySpace,
        class BSplines,
        class interpolation_mesh_type,
        ddc::BoundCond BcXmin,
        ddc::BoundCond BcXmax,
        SplineSolver Solver>
void SplineBuilder<
        ExecSpace,
        MemorySpace,
        BSplines,
        interpolation_mesh_type,
        BcXmin,
        BcXmax,
        Solver>::build_matrix_system()
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
    int start = m_interpolation_domain.front().uid();
    ddc::for_each(m_interpolation_domain, [&](auto ix) {
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
} // namespace ddc
