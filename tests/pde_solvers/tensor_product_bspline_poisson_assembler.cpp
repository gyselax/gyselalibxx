// SPDX-License-Identifier: MIT

#include <cstddef>
#include <memory>

#include <ddc/ddc.hpp>
#include <ddc/kernels/splines.hpp>

#include <ginkgo/ginkgo.hpp>
#include <gtest/gtest.h>

#include "ddc_alias_inline_functions.hpp"
#include "ddc_aliases.hpp"
#include "gauss_legendre_integration.hpp"
#include "i_interpolation_builder.hpp"
#include "identity_coordinate_change.hpp"
#include "matrix_batch_csr.hpp"
#include "matrix_utils.hpp"
#include "mesh_builder.hpp"
#include "metric_tensor_evaluator.hpp"
#include "tensor_product_bspline_poisson_assembler.hpp"
#include "vector_index_tools.hpp"

namespace {

using ExecSpace = Kokkos::DefaultExecutionSpace;
using MatrixType = MatrixBatchCsr<ExecSpace, MatrixBatchCsrSolver::CG>;

constexpr std::size_t s_degree = 3;
// The number of cells differs between dimensions to reveal errors in the index ordering.
constexpr std::size_t s_ncells_1 = 32;
constexpr std::size_t s_ncells_2 = 40;

struct X
{
    static bool constexpr IS_COVARIANT = true;
    static bool constexpr IS_CONTRAVARIANT = true;
    using Dual = X;
};

struct Y
{
    static bool constexpr IS_COVARIANT = true;
    static bool constexpr IS_CONTRAVARIANT = true;
    using Dual = Y;
};

template <bool Periodic>
struct Dim1Cov;

template <bool Periodic>
struct Dim1
{
    static bool constexpr PERIODIC = Periodic;
    static bool constexpr IS_COVARIANT = false;
    static bool constexpr IS_CONTRAVARIANT = true;
    using Dual = Dim1Cov<Periodic>;
};

template <bool Periodic>
struct Dim1Cov
{
    static bool constexpr PERIODIC = Periodic;
    static bool constexpr IS_COVARIANT = true;
    static bool constexpr IS_CONTRAVARIANT = false;
    using Dual = Dim1<Periodic>;
};

template <bool Periodic>
struct Dim2Cov;

template <bool Periodic>
struct Dim2
{
    static bool constexpr PERIODIC = Periodic;
    static bool constexpr IS_COVARIANT = false;
    static bool constexpr IS_CONTRAVARIANT = true;
    using Dual = Dim2Cov<Periodic>;
};

template <bool Periodic>
struct Dim2Cov
{
    static bool constexpr PERIODIC = Periodic;
    static bool constexpr IS_COVARIANT = true;
    static bool constexpr IS_CONTRAVARIANT = false;
    using Dual = Dim2<Periodic>;
};

// Non-uniform B-splines are used because their knots are repeated at a non-periodic boundary.
// The end B-spline is therefore the only B-spline which is non-zero on that boundary.
template <bool Periodic>
struct BSplines1 : ddc::NonUniformBSplines<Dim1<Periodic>, s_degree, Periodic>
{
};

template <bool Periodic>
struct BSplines2 : ddc::NonUniformBSplines<Dim2<Periodic>, s_degree, Periodic>
{
};

template <bool Periodic>
inline constexpr ddc::SplineBuilderClosure s_closure
        = Periodic ? ddc::SplineBuilderClosure::PERIODIC : ddc::SplineBuilderClosure::GREVILLE;

template <bool Periodic>
using InterpPoints1 = ddc::
        GrevilleInterpolationPoints<BSplines1<Periodic>, s_closure<Periodic>, s_closure<Periodic>>;

template <bool Periodic>
using InterpPoints2 = ddc::
        GrevilleInterpolationPoints<BSplines2<Periodic>, s_closure<Periodic>, s_closure<Periodic>>;

template <bool Periodic>
struct Grid1 : NonUniformGridBase<Dim1<Periodic>>
{
};

template <bool Periodic>
struct Grid2 : NonUniformGridBase<Dim2<Periodic>>
{
};

template <bool Periodic>
struct QDim1 : NonUniformGridBase<Dim1<Periodic>>
{
};

template <bool Periodic>
struct QDim2 : NonUniformGridBase<Dim2<Periodic>>
{
};

/// A function which is non-zero on the boundary of the domain and periodic in both dimensions.
struct CosCos
{
    KOKKOS_FUNCTION double operator()(double x, double y) const
    {
        return Kokkos::cos(x) * Kokkos::cos(2 * y);
    }
};

/// A function which is zero on the boundary of the domain and periodic in both dimensions.
struct SinSin
{
    KOKKOS_FUNCTION double operator()(double x, double y) const
    {
        return Kokkos::sin(x) * Kokkos::sin(2 * y);
    }
};

/// Initialise the discrete spaces of the B-splines on uniform break points.
template <bool Periodic1, bool Periodic2>
void init_bsplines()
{
    ddc::init_discrete_space<BSplines1<Periodic1>>(build_uniform_break_points(
            Coord<Dim1<Periodic1>>(0.0),
            Coord<Dim1<Periodic1>>(2 * M_PI),
            IdxStep<Grid1<Periodic1>>(s_ncells_1)));
    ddc::init_discrete_space<BSplines2<Periodic2>>(build_uniform_break_points(
            Coord<Dim2<Periodic2>>(0.0),
            Coord<Dim2<Periodic2>>(2 * M_PI),
            IdxStep<Grid2<Periodic2>>(s_ncells_2)));
}

/// Get the index range of all the distinct B-splines of a basis.
template <class BSplines>
IdxRange<BSplines> get_full_basis()
{
    return ddc::discrete_space<BSplines>().full_domain().take_first(
            IdxStep<BSplines>(ddc::discrete_space<BSplines>().nbasis()));
}

/// Get the index range of the B-splines of a basis which are zero on the boundary of the domain.
template <class BSplines>
IdxRange<BSplines> get_basis_without_boundary()
{
    IdxRange<BSplines> full_basis = get_full_basis<BSplines>();
    if constexpr (BSplines::is_periodic()) {
        return full_basis;
    } else {
        return full_basis.remove(IdxStep<BSplines>(1), IdxStep<BSplines>(1));
    }
}

/**
 * @brief Calculate c^T A c with c the coefficients and A the assembled matrix.
 *
 * @param[in] matrix The assembled matrix.
 * @param[in] coeffs The contiguous spline coefficients on the B-splines of the FEM basis.
 *
 * @return The squared norm induced by the matrix.
 */
template <class BSplinesA, class BSplinesB>
double squared_norm_from_matrix(
        std::unique_ptr<MatrixType> const& matrix,
        DField<IdxRange<BSplinesA, BSplinesB>> coeffs)
{
    auto [values, col_idx, nnz_per_row] = matrix->get_batch_csr();
    std::shared_ptr const gko_exec = gko::ext::kokkos::create_executor(ExecSpace());
    int const n_basis = coeffs.size();

    auto gko_matrix = gko::batch::matrix::Csr<double, int>::
            create(gko_exec,
                   gko::batch_dim<2>(1, gko::dim<2>(n_basis, n_basis)),
                   gko::array<double>::view(gko_exec, values.span(), values.data()),
                   gko::array<int>::view(gko_exec, col_idx.span(), col_idx.data()),
                   gko::array<int>::view(gko_exec, nnz_per_row.span(), nnz_per_row.data()));

    Kokkos::View<double**, Kokkos::LayoutRight, ExecSpace>
            coeffs_view(coeffs.data_handle(), 1, n_basis);
    Kokkos::View<double**, Kokkos::LayoutRight, ExecSpace>
            matrix_times_coeffs_view("matrix_times_coeffs", 1, n_basis);
    auto gko_coeffs = to_gko_multivector(gko_exec, coeffs_view);
    auto gko_matrix_times_coeffs = to_gko_multivector(gko_exec, matrix_times_coeffs_view);

    gko_matrix->apply(gko_coeffs, gko_matrix_times_coeffs);

    auto result = gko::batch::MultiVector<
            double>::create(gko_exec->get_master(), gko::batch_dim<2>(1, gko::dim<2>(1, 1)));
    gko_coeffs->compute_dot(gko_matrix_times_coeffs, result);
    return result->get_const_values()[0];
}

/**
 * @brief Check that the matrix assembled on a basis gives the expected H1-like norm of a
 * spline interpolating the function f(x, y) = g(x) h(2y), with g and h either both sin or
 * both cos, on [0, 2 pi]^2.
 *
 * @f$ \int f^2 dx dy = \pi^2 @f$ and @f$ \int |\nabla f|^2 dx dy = 5 \pi^2 @f$.
 *
 * The B-splines must be initialised before calling this function (see init_bsplines).
 *
 * @param[in] function The function to be interpolated. It must be expressible on the basis.
 * @param[in] idx_range_basis The index range of the tensor product B-splines in the FEM basis.
 */
template <bool Periodic1, bool Periodic2, class Function>
void check_h1_norm(
        Function function,
        IdxRange<BSplines1<Periodic1>, BSplines2<Periodic2>> idx_range_basis)
{
    using D1 = Dim1<Periodic1>;
    using D2 = Dim2<Periodic2>;
    using BS1 = BSplines1<Periodic1>;
    using BS2 = BSplines2<Periodic2>;
    using G1 = Grid1<Periodic1>;
    using G2 = Grid2<Periodic2>;
    using Q1 = QDim1<Periodic1>;
    using Q2 = QDim2<Periodic2>;
    using IdxRangeBS2D = IdxRange<BS1, BS2>;
    using IdxRangeGrid2D = IdxRange<G1, G2>;
    using Mapping = IdentityCoordinateChange<VectorIndexSet<D1, D2>, VectorIndexSet<X, Y>>;
    using Metric = MetricTensorEvaluator<Mapping, Coord<D1, D2>>;
    using Assembler = TensorProductBSplinePoissonAssembler<G1, G2, BS1, BS2, Q1, Q2, Metric>;
    using Builder = ddc::SplineBuilder2D<
            ExecSpace,
            typename ExecSpace::memory_space,
            BS1,
            BS2,
            G1,
            G2,
            s_closure<Periodic1>,
            s_closure<Periodic1>,
            s_closure<Periodic2>,
            s_closure<Periodic2>,
            ddc::SplineSolver::LAPACK>;

    ddc::init_discrete_space<G1>(InterpPoints1<Periodic1>::template get_sampling<G1>());
    ddc::init_discrete_space<G2>(InterpPoints2<Periodic2>::template get_sampling<G2>());

    IdxRangeGrid2D idx_range_grid(
            InterpPoints1<Periodic1>::template get_domain<G1>(),
            InterpPoints2<Periodic2>::template get_domain<G2>());

    // Interpolate the function
    DFieldMem<IdxRangeGrid2D> vals_alloc(idx_range_grid);
    DField<IdxRangeGrid2D> vals = get_field(vals_alloc);
    ddc::parallel_for_each(
            ExecSpace(),
            idx_range_grid,
            KOKKOS_LAMBDA(Idx<G1, G2> const idx) {
                double const x = ddc::coordinate(Idx<G1>(idx));
                double const y = ddc::coordinate(Idx<G2>(idx));
                vals(idx) = function(x, y);
            });

    Builder const builder(idx_range_grid);
    DFieldMem<IdxRangeBS2D> all_coeffs(batched_basis_idx_range(builder, idx_range_grid));
    builder(get_field(all_coeffs), get_const_field(vals_alloc));

    // Restrict the coefficients to the FEM basis while avoiding strided layouts
    DFieldMem<IdxRangeBS2D> coeffs(idx_range_basis);
    ddc::parallel_deepcopy(get_field(coeffs), get_const_field(all_coeffs)[idx_range_basis]);

    // Define the integral volume on the physical domain (the mapping is the identity)
    GaussLegendre<Q1, s_degree + 1> gl_coeffs_1(ddc::discrete_space<BS1>().break_point_domain());
    GaussLegendre<Q2, s_degree + 1> gl_coeffs_2(ddc::discrete_space<BS2>().break_point_domain());
    DFieldMem<IdxRange<Q1, Q2>> int_volume
            = gauss_legendre_quadrature_coefficients<ExecSpace>(gl_coeffs_1, gl_coeffs_2);

    // Assemble the matrix
    double const alpha = 2.0;
    double const beta = 0.5;
    Metric const metric {Mapping {}};
    Assembler const assembler(get_field(int_volume), metric, idx_range_grid, idx_range_basis);
    std::unique_ptr<MatrixType> matrix;
    assembler.setup_sparse_matrix(matrix);
    assembler(
            matrix,
            KOKKOS_LAMBDA(Idx<>, Coord<D1, D2>) { return alpha; },
            KOKKOS_LAMBDA(Idx<>, Coord<D1, D2>) { return beta; });

    double const norm = squared_norm_from_matrix(matrix, get_field(coeffs));
    double const expected_norm = (5 * alpha + beta) * M_PI * M_PI;
    // The interpolation error dominates. It is O(h^4) for cubic splines.
    EXPECT_NEAR(norm / expected_norm, 1.0, 2e-4);
}

template <bool Periodic1, bool Periodic2>
void check_h1_norm_full_basis()
{
    init_bsplines<Periodic1, Periodic2>();
    check_h1_norm<Periodic1, Periodic2>(
            CosCos {},
            IdxRange<BSplines1<Periodic1>, BSplines2<Periodic2>>(
                    get_full_basis<BSplines1<Periodic1>>(),
                    get_full_basis<BSplines2<Periodic2>>()));
}

template <bool Periodic1, bool Periodic2>
void check_h1_norm_without_boundary_basis()
{
    init_bsplines<Periodic1, Periodic2>();
    check_h1_norm<Periodic1, Periodic2>(
            SinSin {},
            IdxRange<BSplines1<Periodic1>, BSplines2<Periodic2>>(
                    get_basis_without_boundary<BSplines1<Periodic1>>(),
                    get_basis_without_boundary<BSplines2<Periodic2>>()));
}

TEST(TensorProductBSplinePoissonAssembler, FullBasisNonPeriodicNonPeriodic)
{
    check_h1_norm_full_basis<false, false>();
}

TEST(TensorProductBSplinePoissonAssembler, FullBasisNonPeriodicPeriodic)
{
    check_h1_norm_full_basis<false, true>();
}

TEST(TensorProductBSplinePoissonAssembler, FullBasisPeriodicNonPeriodic)
{
    check_h1_norm_full_basis<true, false>();
}

TEST(TensorProductBSplinePoissonAssembler, FullBasisPeriodicPeriodic)
{
    check_h1_norm_full_basis<true, true>();
}

TEST(TensorProductBSplinePoissonAssembler, NoBoundaryBasisNonPeriodicNonPeriodic)
{
    check_h1_norm_without_boundary_basis<false, false>();
}

TEST(TensorProductBSplinePoissonAssembler, NoBoundaryBasisNonPeriodicPeriodic)
{
    check_h1_norm_without_boundary_basis<false, true>();
}

TEST(TensorProductBSplinePoissonAssembler, NoBoundaryBasisPeriodicNonPeriodic)
{
    check_h1_norm_without_boundary_basis<true, false>();
}

} // namespace
