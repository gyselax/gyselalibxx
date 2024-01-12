#pragma once

#include <cstdint>
#include <iostream>
#include <memory>
#include <optional>
#include <stdexcept>
#include <vector>

#include <ddc/ddc.hpp>

#include <Kokkos_Core.hpp>

#include "sll/math_tools.hpp"
#include "sll/matrix.hpp"
#include "sll/spline_boundary_conditions.hpp"
#include "sll/view.hpp"

constexpr bool is_spline_interpolation_mesh_uniform(
        bool const is_uniform,
        BoundCond const BcXmin,
        BoundCond const BcXmax,
        int degree)
{
    int N_BE_MIN = n_boundary_equations(BcXmin, degree);
    int N_BE_MAX = n_boundary_equations(BcXmax, degree);
    bool is_periodic = (BcXmin == BoundCond::PERIODIC) && (BcXmax == BoundCond::PERIODIC);
    return is_uniform && ((N_BE_MIN != 0 && N_BE_MAX != 0) || is_periodic);
}

/**
 * @brief A class for creating a spline approximation of a function.
 *
 * A class which contains an operator () which can be used to build a spline approximation
 * of a function. A spline approximation is represented by coefficients stored in a Chunk
 * of BSplines. The spline is constructed such that it respects the boundary conditions
 * BcXmin and BcXmax, and it interpolates the function at the points on the interpolation_mesh
 * associated with interpolation_mesh_type.
 */
template <class BSplines, class interpolation_mesh_type, BoundCond BcXmin, BoundCond BcXmax>
class SplineBuilder
{
    static_assert(
            (BSplines::is_periodic() && (BcXmin == BoundCond::PERIODIC)
             && (BcXmax == BoundCond::PERIODIC))
            || (!BSplines::is_periodic() && (BcXmin != BoundCond::PERIODIC)
                && (BcXmax != BoundCond::PERIODIC)));
    static_assert(!BSplines::is_radial());

private:
    using tag_type = typename interpolation_mesh_type::continuous_dimension_type;

public:
    /**
     * @brief The type of the BSplines which are compatible with this class.
     */
    using bsplines_type = BSplines;

    /**
     * @brief The type of the interpolation mesh used by this class.
     */
    using mesh_type = interpolation_mesh_type;

    /**
     * @brief The type of the domain for the interpolation mesh used by this class.
     */
    using interpolation_domain_type = ddc::DiscreteDomain<interpolation_mesh_type>;

public:
    /**
     * @brief Indicates if the degree of the splines is odd or even.
     */
    static constexpr bool s_odd = BSplines::degree() % 2;

    /**
     * @brief The number of equations which define the boundary conditions at the lower bound.
     */
    static constexpr int s_nbe_xmin = n_boundary_equations(BcXmin, BSplines::degree());

    /**
     * @brief The number of equations which define the boundary conditions at the upper bound.
     */
    static constexpr int s_nbe_xmax = n_boundary_equations(BcXmax, BSplines::degree());

    /**
     * @brief The number of boundary conditions which must be provided by the user at the lower bound.
     *
     * This value is usually equal to s_nbe_xmin, but it may be difference if the chosen boundary
     * conditions impose a specific value (e.g. no values need to be provided for Dirichlet boundary
     * conditions).
     */
    static constexpr int s_nbc_xmin = n_user_input(BcXmin, BSplines::degree());

    /**
     * @brief The number of boundary conditions which must be provided by the user at the upper bound.
     *
     * This value is usually equal to s_nbe_xmin, but it may be difference if the chosen boundary
     * conditions impose a specific value (e.g. no values need to be provided for Dirichlet boundary
     * conditions).
     */
    static constexpr int s_nbc_xmax = n_user_input(BcXmax, BSplines::degree());

    /**
     * @brief The boundary condition implemented at the lower bound.
     */
    static constexpr BoundCond s_bc_xmin = BcXmin;
    /**
     * @brief The boundary condition implemented at the upper bound.
     */
    static constexpr BoundCond s_bc_xmax = BcXmax;

private:
    interpolation_domain_type m_interpolation_domain;

    double m_dx; // average cell size for normalization of derivatives

    // interpolator specific
    std::unique_ptr<Matrix> matrix;

    int m_offset;

public:
    /**
     * @brief Create a new SplineBuilder.
     *
     * @param interpolation_domain The domain on which points will be provided in order to
     * 		create the spline approximation.
     */
    SplineBuilder(interpolation_domain_type const& interpolation_domain);

    /**
     * @brief Create a new SplineBuilder by copy
     *
     * @param x The SplineBuilder being copied.
     */
    SplineBuilder(SplineBuilder const& x) = delete;

    /**
     * @brief Create a new SplineBuilder by copy
     *
     * @param x The SplineBuilder being copied.
     */
    SplineBuilder(SplineBuilder&& x) = default;

    ~SplineBuilder() = default;

    SplineBuilder& operator=(SplineBuilder const& x) = delete;

    /**
     * @brief Copy a SplineBuilder.
     *
     * @param x The SplineBuilder being copied.
     * @returns A reference to this object.
     */
    SplineBuilder& operator=(SplineBuilder&& x) = default;

    /**
     * @brief Build a spline approximation of a function.
     *
     * Use the values of a function at known grid points (as specified by
     * SplineBuilder::interpolation_domain) and the derivatives of the
     * function at the boundaries (if necessary for the chosen boundary
     * conditions) to calculate a spline approximation of a function.
     *
     * The spline approximation is stored as a ChunkSpan of coefficients
     * associated with basis-splines.
     *
     * @param[out] spline The coefficients of the spline calculated by the function.
     * @param[in] vals The values of the function at the grid points.
     * @param[in] derivs_xmin The values of the derivatives at the lower boundary.
     * @param[in] derivs_xmax The values of the derivatives at the upper boundary.
     */
    void operator()(
            ddc::ChunkSpan<double, ddc::DiscreteDomain<bsplines_type>> spline,
            ddc::ChunkSpan<double const, interpolation_domain_type> vals,
            std::optional<CDSpan1D> const derivs_xmin = std::nullopt,
            std::optional<CDSpan1D> const derivs_xmax = std::nullopt) const;

    /**
     * @brief Get the domain from which the approximation is defined.
     *
     * Get the domain on which values of the function must be provided in order
     * to build a spline approximation of the function.
     *
     * @return The domain for the grid points.
     */
    interpolation_domain_type const& interpolation_domain() const noexcept
    {
        return m_interpolation_domain;
    }

    /**
     * @brief Get the domain on which the approximation is defined.
     *
     * Get the domain of the basis-splines for which the coefficients of the spline
     * approximation must be calculated.
     *
     * @return The domain for the splines.
     */
    ddc::DiscreteDomain<BSplines> spline_domain() const noexcept
    {
        return ddc::discrete_space<BSplines>().full_domain();
    }

    /**
     * @brief Get the interpolation matrix.
     *
     * Get the interpolation matrix. This can be useful for debugging (as it allows
     * one to print the matrix) or for more complex quadrature schemes.
     *
     * @return A reference to the interpolation matrix.
     */
    const Matrix& get_interpolation_matrix() const noexcept
    {
        return *matrix;
    }

    /**
     * @brief Get the spline quadrature coefficients.
     *
     * To integrate a function with a spline quadrature, we use:
     *
     * @f$ \int_a^b f(x)dx
     * \simeq \sum_{i = 0}^{N_{\text{basis}} -1 } c_i  \int_a^b b_{i,d}()x dx @f$,
     *
     * which rewritten gives
     *
     * @f$ \int_a^b f(x)dx
     * \simeq \sum_{i = 0}^{N_{\text{basis}} - 1} q_i f_i @f$,
     *
     * with
     *  - @f$\{ f_i\}_i @f$ the values of the function at the interpolation points;
     *  - @f$ q = \{ q_i\}_i @f$ the quadrature coefficients we compute thanks to
     *  @f$ q B^T = I_b @f$,
     *      - with @f$ B @f$ the matrix of B-splines @f$ B_{ij} = b_{j,d}(x_i)@f$,
     *      - and @f$ I_b = \int_a^b b_{i,d}(x)dx @f$ the integrated B-splines.
     *
     * More details are given in Emily Bourne's thesis
     * "Non-Uniform Numerical Schemes for the Modelling of Turbulence
     * in the 5D GYSELA Code". December 2022.
     *
     *
     * @param[in] domain
     *      The domain where the functions we want to integrate
     *      are defined.
     *
     * @return A chunk with the quadrature coefficients @f$ q @f$.
     */
    template <class IDim>
    ddc::Chunk<double, ddc::DiscreteDomain<IDim>> quadrature_coefficients(
            ddc::DiscreteDomain<IDim> const& domain) const noexcept
    {
        static_assert(s_nbe_xmin == 0 && s_nbe_xmax == 0);
        ddc::Chunk<double, ddc::DiscreteDomain<IDim>> coefficients(domain);

        // Vector of integrals of B-splines
        ddc::Chunk<double, ddc::DiscreteDomain<bsplines_type>> integral_bsplines(spline_domain());
        ddc::discrete_space<bsplines_type>().integrals(integral_bsplines);

        // Coefficients of quadrature in integral_bsplines
        ddc::DiscreteDomain<bsplines_type> slice = spline_domain().take_first(
                ddc::DiscreteVector<bsplines_type> {ddc::discrete_space<BSplines>().nbasis()});

        Kokkos::deep_copy(
                coefficients.allocation_kokkos_view(),
                integral_bsplines[slice].allocation_kokkos_view());

        matrix->solve_transpose_inplace(coefficients.allocation_mdspan());

        return coefficients;
    }


private:
    void compute_block_sizes_uniform(int& lower_block_size, int& upper_block_size) const;

    void compute_block_sizes_non_uniform(int& lower_block_size, int& upper_block_size) const;

    void allocate_matrix(int lower_block_size, int upper_block_size);

    void compute_interpolant_degree1(
            ddc::ChunkSpan<double, ddc::DiscreteDomain<bsplines_type>> spline,
            ddc::ChunkSpan<double const, interpolation_domain_type> vals) const;

    void build_matrix_system();
};

template <class BSplines, class interpolation_mesh_type, BoundCond BcXmin, BoundCond BcXmax>
SplineBuilder<BSplines, interpolation_mesh_type, BcXmin, BcXmax>::SplineBuilder(
        interpolation_domain_type const& interpolation_domain)
    : m_interpolation_domain(interpolation_domain)
    , m_dx((ddc::discrete_space<BSplines>().rmax() - ddc::discrete_space<BSplines>().rmin())
           / ddc::discrete_space<BSplines>().ncells())
    , matrix(nullptr)
    , m_offset(0)
{
    if constexpr (bsplines_type::is_periodic()) {
        // Calculate offset so that the matrix is diagonally dominant
        std::array<double, bsplines_type::degree() + 1> values_ptr;
        DSpan1D values(values_ptr.data(), bsplines_type::degree() + 1);
        ddc::DiscreteElement<interpolation_mesh_type> start(interpolation_domain.front());
        auto jmin = ddc::discrete_space<BSplines>()
                            .eval_basis(values, ddc::coordinate(start + BSplines::degree()));
        if constexpr (bsplines_type::degree() % 2 == 0) {
            m_offset = jmin.uid() - start.uid() + bsplines_type::degree() / 2 - BSplines::degree();
        } else {
            int const mid = bsplines_type::degree() / 2;
            m_offset = jmin.uid() - start.uid() + (values(mid) > values(mid + 1) ? mid : mid + 1)
                       - BSplines::degree();
        }
    }

    // Calculate block sizes
    int lower_block_size, upper_block_size;
    if constexpr (bsplines_type::is_uniform()) {
        compute_block_sizes_uniform(lower_block_size, upper_block_size);
    } else {
        compute_block_sizes_non_uniform(lower_block_size, upper_block_size);
    }
    allocate_matrix(lower_block_size, upper_block_size);
}

//-------------------------------------------------------------------------------------------------
/************************************************************************************
 *                         Compute interpolant functions *
 ************************************************************************************/

template <class BSplines, class interpolation_mesh_type, BoundCond BcXmin, BoundCond BcXmax>
void SplineBuilder<BSplines, interpolation_mesh_type, BcXmin, BcXmax>::compute_interpolant_degree1(
        ddc::ChunkSpan<double, ddc::DiscreteDomain<bsplines_type>> const spline,
        ddc::ChunkSpan<double const, interpolation_domain_type> const vals) const
{
    for (std::size_t i = 0; i < ddc::discrete_space<BSplines>().nbasis(); ++i) {
        spline(ddc::DiscreteElement<bsplines_type>(i))
                = vals(ddc::DiscreteElement<interpolation_mesh_type>(i));
    }
    if constexpr (bsplines_type::is_periodic()) {
        spline(ddc::DiscreteElement<bsplines_type>(ddc::discrete_space<BSplines>().nbasis()))
                = spline(ddc::DiscreteElement<bsplines_type>(0));
    }
}

//-------------------------------------------------------------------------------------------------

template <class BSplines, class interpolation_mesh_type, BoundCond BcXmin, BoundCond BcXmax>
void SplineBuilder<BSplines, interpolation_mesh_type, BcXmin, BcXmax>::operator()(
        ddc::ChunkSpan<double, ddc::DiscreteDomain<bsplines_type>> const spline,
        ddc::ChunkSpan<double const, interpolation_domain_type> const vals,
        std::optional<CDSpan1D> const derivs_xmin,
        std::optional<CDSpan1D> const derivs_xmax) const
{
    assert(vals.template extent<interpolation_mesh_type>()
           == ddc::discrete_space<BSplines>().nbasis() - s_nbe_xmin - s_nbe_xmax);
    // assert(spline.belongs_to_space(ddc::discrete_space<BSplines>()));
    // TODO: LOG Errors
    if constexpr (bsplines_type::degree() == 1)
        return compute_interpolant_degree1(spline, vals);

    assert((BcXmin == BoundCond::HERMITE)
           != (!derivs_xmin.has_value() || derivs_xmin->extent(0) == 0));
    assert((BcXmax == BoundCond::HERMITE)
           != (!derivs_xmax.has_value() || derivs_xmax->extent(0) == 0));

    // Hermite boundary conditions at xmin, if any
    // NOTE: For consistency with the linear system, the i-th derivative
    //       provided by the user must be multiplied by dx^i
    if constexpr (BcXmin == BoundCond::HERMITE) {
        assert(derivs_xmin->extent(0) == s_nbc_xmin);
        for (int i = s_nbc_xmin; i > 0; --i) {
            spline(ddc::DiscreteElement<bsplines_type>(s_nbc_xmin - i))
                    = (*derivs_xmin)(i - 1) * ipow(m_dx, i + s_odd - 1);
        }
    }
    for (int i = s_nbc_xmin; i < s_nbc_xmin + m_offset; ++i) {
        spline(ddc::DiscreteElement<bsplines_type>(i)) = 0.0;
    }

    for (int i = 0; i < m_interpolation_domain.extents(); ++i) {
        spline(ddc::DiscreteElement<bsplines_type>(s_nbc_xmin + i + m_offset))
                = vals(ddc::DiscreteElement<interpolation_mesh_type>(i));
    }

    // Hermite boundary conditions at xmax, if any
    // NOTE: For consistency with the linear system, the i-th derivative
    //       provided by the user must be multiplied by dx^i
    if constexpr (BcXmax == BoundCond::HERMITE) {
        assert(derivs_xmax->extent(0) == s_nbc_xmax);
        for (int i = 0; i < s_nbc_xmax; ++i) {
            spline(ddc::DiscreteElement<bsplines_type>(
                    ddc::discrete_space<BSplines>().nbasis() - s_nbc_xmax + i))
                    = (*derivs_xmax)(i)*ipow(m_dx, i + s_odd);
        }
    }

    DSpan1D const bcoef_section(
            spline.data_handle() + m_offset,
            ddc::discrete_space<BSplines>().nbasis());
    matrix->solve_inplace(bcoef_section);

    if constexpr (bsplines_type::is_periodic()) {
        if (m_offset != 0) {
            for (int i = 0; i < m_offset; ++i) {
                spline(ddc::DiscreteElement<bsplines_type>(i))
                        = spline(ddc::DiscreteElement<bsplines_type>(
                                ddc::discrete_space<BSplines>().nbasis() + i));
            }
            for (std::size_t i = m_offset; i < bsplines_type::degree(); ++i) {
                spline(ddc::DiscreteElement<bsplines_type>(
                        ddc::discrete_space<BSplines>().nbasis() + i))
                        = spline(ddc::DiscreteElement<bsplines_type>(i));
            }
        }
    }
}

//-------------------------------------------------------------------------------------------------
/************************************************************************************
 *                            Compute num diags functions *
 ************************************************************************************/

template <class BSplines, class interpolation_mesh_type, BoundCond BcXmin, BoundCond BcXmax>
void SplineBuilder<BSplines, interpolation_mesh_type, BcXmin, BcXmax>::compute_block_sizes_uniform(
        int& lower_block_size,
        int& upper_block_size) const
{
    switch (BcXmin) {
    case BoundCond::PERIODIC:
        upper_block_size = (bsplines_type::degree()) / 2;
        break;
    case BoundCond::NATURAL:
    case BoundCond::HERMITE:
        upper_block_size = s_nbc_xmin;
        break;
    case BoundCond::GREVILLE:
        upper_block_size = bsplines_type::degree() - 1;
        break;
    default:
        throw std::runtime_error("BoundCond not handled");
    }
    switch (BcXmax) {
    case BoundCond::PERIODIC:
        lower_block_size = (bsplines_type::degree()) / 2;
        break;
    case BoundCond::NATURAL:
    case BoundCond::HERMITE:
        lower_block_size = s_nbc_xmax;
        break;
    case BoundCond::GREVILLE:
        lower_block_size = bsplines_type::degree() - 1;
        break;
    default:
        throw std::runtime_error("BoundCond not handled");
    }
}

//-------------------------------------------------------------------------------------------------

template <class BSplines, class interpolation_mesh_type, BoundCond BcXmin, BoundCond BcXmax>
void SplineBuilder<BSplines, interpolation_mesh_type, BcXmin, BcXmax>::
        compute_block_sizes_non_uniform(int& lower_block_size, int& upper_block_size) const
{
    switch (BcXmin) {
    case BoundCond::PERIODIC:
        upper_block_size = bsplines_type::degree() - 1;
        break;
    case BoundCond::HERMITE:
        upper_block_size = s_nbc_xmin + 1;
        break;
    case BoundCond::GREVILLE:
        upper_block_size = bsplines_type::degree() - 1;
        break;
    default:
        throw std::runtime_error("BoundCond not handled");
    }
    switch (BcXmax) {
    case BoundCond::PERIODIC:
        lower_block_size = bsplines_type::degree() - 1;
        break;
    case BoundCond::HERMITE:
        lower_block_size = s_nbc_xmax + 1;
        break;
    case BoundCond::GREVILLE:
        lower_block_size = bsplines_type::degree() - 1;
        break;
    default:
        throw std::runtime_error("BoundCond not handled");
    }
}

//-------------------------------------------------------------------------------------------------
/************************************************************************************
 *                            Initialize matrix functions *
 ************************************************************************************/

template <class BSplines, class interpolation_mesh_type, BoundCond BcXmin, BoundCond BcXmax>
void SplineBuilder<BSplines, interpolation_mesh_type, BcXmin, BcXmax>::allocate_matrix(
        int lower_block_size,
        int upper_block_size)
{
    // Special case: linear spline
    // No need for matrix assembly
    if constexpr (bsplines_type::degree() == 1)
        return;

    int upper_band_width;
    if (bsplines_type::is_uniform()) {
        upper_band_width = bsplines_type::degree() / 2;
    } else {
        upper_band_width = bsplines_type::degree() - 1;
    }

    if constexpr (bsplines_type::is_periodic()) {
        matrix = Matrix::make_new_periodic_banded(
                ddc::discrete_space<BSplines>().nbasis(),
                upper_band_width,
                upper_band_width,
                bsplines_type::is_uniform());
    } else {
        matrix = Matrix::make_new_block_with_banded_region(
                ddc::discrete_space<BSplines>().nbasis(),
                upper_band_width,
                upper_band_width,
                bsplines_type::is_uniform(),
                upper_block_size,
                lower_block_size);
    }

    build_matrix_system();

    matrix->factorize();
}

//-------------------------------------------------------------------------------------------------

template <class BSplines, class interpolation_mesh_type, BoundCond BcXmin, BoundCond BcXmax>
void SplineBuilder<BSplines, interpolation_mesh_type, BcXmin, BcXmax>::build_matrix_system()
{
    // Hermite boundary conditions at xmin, if any
    if constexpr (BcXmin == BoundCond::HERMITE) {
        double derivs_ptr[(bsplines_type::degree() / 2 + 1) * (bsplines_type::degree() + 1)];
        DSpan2D derivs(derivs_ptr, bsplines_type::degree() + 1, bsplines_type::degree() / 2 + 1);
        ddc::discrete_space<BSplines>().eval_basis_and_n_derivs(
                derivs,
                ddc::discrete_space<BSplines>().rmin(),
                s_nbc_xmin);

        // In order to improve the condition number of the matrix, we normalize
        // all derivatives by multiplying the i-th derivative by dx^i
        for (std::size_t i = 0; i < bsplines_type::degree() + 1; ++i) {
            for (std::size_t j = 1; j < bsplines_type::degree() / 2 + 1; ++j) {
                derivs(i, j) *= ipow(m_dx, j);
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
    std::array<double, bsplines_type::degree() + 1> values_ptr;
    std::experimental::mdspan<
            double,
            std::experimental::extents<std::size_t, bsplines_type::degree() + 1>> const
            values(values_ptr.data());
    int start = m_interpolation_domain.front().uid();
    ddc::for_each(m_interpolation_domain, [&](auto ix) {
        auto jmin = ddc::discrete_space<BSplines>().eval_basis(
                values,
                ddc::coordinate(ddc::DiscreteElement<interpolation_mesh_type>(ix)));
        for (std::size_t s = 0; s < bsplines_type::degree() + 1; ++s) {
            int const j = modulo(
                    int(jmin.uid() - m_offset + s),
                    (int)ddc::discrete_space<BSplines>().nbasis());
            matrix->set_element(ix.uid() - start + s_nbc_xmin, j, values(s));
        }
    });

    // Hermite boundary conditions at xmax, if any
    if constexpr (BcXmax == BoundCond::HERMITE) {
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
                derivs(i, j) *= ipow(m_dx, j);
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
