// SPDX-License-Identifier: MIT
#include <optional>
#include <type_traits>

#include <ddc/ddc.hpp>
#include <ddc/kernels/splines.hpp>

#include <gtest/gtest.h>

#include "composite_spline_builder.hpp"
#include "ddc_alias_inline_functions.hpp"
#include "ddc_aliases.hpp"
#include "i_interpolation_builder.hpp"
#include "type_seq_tools.hpp"

namespace {

/*
 * The test function is a product of 1D functions sin(k * x + phase). The wavenumber k
 * and the phase are defined for each dimension. In periodic dimensions the wavenumber
 * is chosen so the function is periodic over the domain.
 */
struct X
{
    static bool constexpr PERIODIC = false;
    static double constexpr wavenumber = 1.3;
    static double constexpr phase = 0.2;
};
struct Y
{
    static bool constexpr PERIODIC = true;
    // The domain is [0, 2]
    static double constexpr wavenumber = Kokkos::numbers::pi;
    static double constexpr phase = 0.7;
};
struct Z
{
    static bool constexpr PERIODIC = false;
    static double constexpr wavenumber = 0.9;
    static double constexpr phase = 0.4;
};
struct Batch
{
};

int constexpr Degree = 3;

struct BSplinesX : ddc::UniformBSplines<X, Degree, X::PERIODIC>
{
};
struct BSplinesY : ddc::UniformBSplines<Y, Degree, Y::PERIODIC>
{
};
struct BSplinesZ : ddc::UniformBSplines<Z, Degree, Z::PERIODIC>
{
};

auto constexpr ClosureX = ddc::SplineBuilderClosure::GREVILLE;
auto constexpr ClosureY = ddc::SplineBuilderClosure::PERIODIC;
auto constexpr ClosureZ = ddc::SplineBuilderClosure::GREVILLE;

using SplineInterpPointsX = ddc::GrevilleInterpolationPoints<BSplinesX, ClosureX, ClosureX>;
using SplineInterpPointsY = ddc::GrevilleInterpolationPoints<BSplinesY, ClosureY, ClosureY>;
using SplineInterpPointsZ = ddc::GrevilleInterpolationPoints<BSplinesZ, ClosureZ, ClosureZ>;

struct GridX : SplineInterpPointsX::interpolation_discrete_dimension_type
{
};
struct GridY : SplineInterpPointsY::interpolation_discrete_dimension_type
{
};
struct GridZ : SplineInterpPointsZ::interpolation_discrete_dimension_type
{
};
struct GridBatch : UniformGridBase<Batch>
{
};

using IdxRangeXY = IdxRange<GridX, GridY>;
using IdxRangeXBatchY = IdxRange<GridX, GridBatch, GridY>;
using IdxRangeXYZ = IdxRange<GridX, GridY, GridZ>;
using IdxRangeXBatchYZ = IdxRange<GridX, GridBatch, GridY, GridZ>;

using ExecSpace = Kokkos::DefaultExecutionSpace;
using MemSpace = ExecSpace::memory_space;

using BuilderX = ddc::SplineBuilder<
        ExecSpace,
        MemSpace,
        BSplinesX,
        GridX,
        ClosureX,
        ClosureX,
        ddc::SplineSolver::LAPACK>;
using BuilderY = ddc::SplineBuilder<
        ExecSpace,
        MemSpace,
        BSplinesY,
        GridY,
        ClosureY,
        ClosureY,
        ddc::SplineSolver::LAPACK>;

using BuilderXY = CompositeSplineBuilder2D<BuilderX, BuilderY>;
using RefBuilderXY = ddc::SplineBuilder2D<
        ExecSpace,
        MemSpace,
        BSplinesX,
        BSplinesY,
        GridX,
        GridY,
        ClosureX,
        ClosureX,
        ClosureY,
        ClosureY,
        ddc::SplineSolver::LAPACK>;

using BuilderZ = ddc::SplineBuilder<
        ExecSpace,
        MemSpace,
        BSplinesZ,
        GridZ,
        ClosureZ,
        ClosureZ,
        ddc::SplineSolver::LAPACK>;
using BuilderYZ = ddc::SplineBuilder2D<
        ExecSpace,
        MemSpace,
        BSplinesY,
        BSplinesZ,
        GridY,
        GridZ,
        ClosureY,
        ClosureY,
        ClosureZ,
        ClosureZ,
        ddc::SplineSolver::LAPACK>;

using BuilderXYZ = CompositeSplineBuilder3D<BuilderX, BuilderYZ>;
using RefBuilderXYZ = ddc::SplineBuilder3D<
        ExecSpace,
        MemSpace,
        BSplinesX,
        BSplinesY,
        BSplinesZ,
        GridX,
        GridY,
        GridZ,
        ClosureX,
        ClosureX,
        ClosureY,
        ClosureY,
        ClosureZ,
        ClosureZ,
        ddc::SplineSolver::LAPACK>;

void init_discrete_spaces()
{
    ddc::init_discrete_space<BSplinesX>(Coord<X>(0.0), Coord<X>(1.0), IdxStep<GridX>(10));
    ddc::init_discrete_space<BSplinesY>(Coord<Y>(0.0), Coord<Y>(2.0), IdxStep<GridY>(12));
    ddc::init_discrete_space<BSplinesZ>(Coord<Z>(0.0), Coord<Z>(1.5), IdxStep<GridZ>(8));
    ddc::init_discrete_space<GridX>(SplineInterpPointsX::get_sampling<GridX>());
    ddc::init_discrete_space<GridY>(SplineInterpPointsY::get_sampling<GridY>());
    ddc::init_discrete_space<GridZ>(SplineInterpPointsZ::get_sampling<GridZ>());
    ddc::init_discrete_space<GridBatch>(
            GridBatch::init(Coord<Batch>(0.0), Coord<Batch>(1.0), IdxStep<GridBatch>(4)));
}

IdxRangeXY get_interpolation_idx_range()
{
    return IdxRangeXY(
            SplineInterpPointsX::get_domain<GridX>(),
            SplineInterpPointsY::get_domain<GridY>());
}

IdxRangeXYZ get_interpolation_idx_range_3d()
{
    return IdxRangeXYZ(
            SplineInterpPointsX::get_domain<GridX>(),
            SplineInterpPointsY::get_domain<GridY>(),
            SplineInterpPointsZ::get_domain<GridZ>());
}

/**
 * Describes the 1D spline builder and interpolation grid used along one dimension.
 * The nested Grid is unique to each instantiation so each configuration has its own grid.
 */
template <class BSplines, ddc::SplineBuilderClosure Closure>
struct SplineConfig
{
    using bsplines_type = BSplines;
    static constexpr ddc::SplineBuilderClosure closure = Closure;
    using InterpPoints = ddc::GrevilleInterpolationPoints<BSplines, Closure, Closure>;
    struct Grid : InterpPoints::interpolation_discrete_dimension_type
    {
    };
    using Builder = ddc::SplineBuilder<
            ExecSpace,
            MemSpace,
            BSplines,
            Grid,
            Closure,
            Closure,
            ddc::SplineSolver::LAPACK>;

    static IdxRange<Grid> init_grid()
    {
        ddc::init_discrete_space<Grid>(InterpPoints::template get_sampling<Grid>());
        return InterpPoints::template get_domain<Grid>();
    }
};

using ConfigXGreville = SplineConfig<BSplinesX, ddc::SplineBuilderClosure::GREVILLE>;
using ConfigXHermite = SplineConfig<BSplinesX, ddc::SplineBuilderClosure::HERMITE>;
using ConfigYPeriodic = SplineConfig<BSplinesY, ddc::SplineBuilderClosure::PERIODIC>;
using ConfigZGreville = SplineConfig<BSplinesZ, ddc::SplineBuilderClosure::GREVILLE>;
using ConfigZHermite = SplineConfig<BSplinesZ, ddc::SplineBuilderClosure::HERMITE>;

void init_bsplines()
{
    ddc::init_discrete_space<BSplinesX>(Coord<X>(0.0), Coord<X>(1.0), IdxStep<GridX>(10));
    ddc::init_discrete_space<BSplinesY>(Coord<Y>(0.0), Coord<Y>(2.0), IdxStep<GridY>(12));
    ddc::init_discrete_space<BSplinesZ>(Coord<Z>(0.0), Coord<Z>(1.5), IdxStep<GridZ>(8));
}

/**
 * Evaluate the 1D test function sin(k * x + phase), or one of its derivatives.
 *
 * @param coord The coordinate at which the function is evaluated.
 * @param deriv_order The order of the derivative.
 * @returns The value of the derivative of the test function.
 */
template <class Dim>
KOKKOS_FUNCTION double test_function_1d(Coord<Dim> coord, int deriv_order = 0)
{
    double const k = Dim::wavenumber;
    return Kokkos::pow(k, deriv_order)
           * Kokkos::sin(
                   k * static_cast<double>(coord) + Dim::phase
                   + deriv_order * Kokkos::numbers::pi / 2);
}

/**
 * Fill a field with the values of the test function at the interpolation points.
 * The test function is constant along the batch dimensions.
 *
 * @tparam Dims The interpolation dimensions.
 */
template <class... Dims, class... Grids>
void fill_with_test_values(DField<IdxRange<Grids...>> vals)
{
    using GridTypeSeq = ddc::detail::TypeSeq<Grids...>;
    ddc::parallel_for_each(
            ExecSpace(),
            get_idx_range(vals),
            KOKKOS_LAMBDA(Idx<Grids...> idx) {
                vals(idx)
                        = (test_function_1d(
                                   ddc::coordinate(Idx<find_grid_t<Dims, GridTypeSeq>>(idx)))
                           * ...);
            });
}

/**
 * Evaluate the 1D test function along Dim at an index. If the index is defined along the
 * derivative dimension ddc::Deriv<Dim> then the derivative of the order given by the index
 * is evaluated at the boundary. Otherwise the test function is evaluated at the coordinate
 * of the index and the boundary is unused.
 */
template <class Dim, class... Grids>
KOKKOS_FUNCTION double test_function_on_dim(Idx<Grids...> idx, Coord<Dim> boundary)
{
    using GridTypeSeq = ddc::detail::TypeSeq<Grids...>;
    if constexpr (ddc::in_tags_v<ddc::Deriv<Dim>, GridTypeSeq>) {
        int const deriv_order = (Idx<ddc::Deriv<Dim>>(idx) - Idx<ddc::Deriv<Dim>>(0)).value();
        return test_function_1d(boundary, deriv_order);
    } else {
        return test_function_1d(ddc::coordinate(Idx<find_grid_t<Dim, GridTypeSeq>>(idx)));
    }
}

/**
 * Fill a field with the values of the test function where the field is defined along
 * derivative dimensions. Along each interpolation dimension the field is defined either on
 * an interpolation grid or on the derivative dimension ddc::Deriv<Dim>. Along a derivative
 * dimension the derivative of the order given by the index is evaluated at the boundary.
 * The test function is constant along any other (batch) dimensions.
 *
 * @param[out] derivs The field to be filled.
 * @param[in] boundaries The boundary of each interpolation dimension at which derivatives
 *          are evaluated. A boundary is unused if the field is not defined along the
 *          derivative dimension of that interpolation dimension.
 */
template <class... Grids, class... Dims>
void fill_with_test_derivs(DField<IdxRange<Grids...>> derivs, Coord<Dims>... boundaries)
{
    ddc::parallel_for_each(
            ExecSpace(),
            get_idx_range(derivs),
            KOKKOS_LAMBDA(Idx<Grids...> idx) {
                derivs(idx) = (test_function_on_dim(idx, boundaries) * ...);
            });
}

/// Get the maximum absolute difference between two fields.
template <class IdxRangeType>
double max_abs_difference(DConstField<IdxRangeType> field1, DConstField<IdxRangeType> field2)
{
    return ddc::parallel_transform_reduce(
            ExecSpace(),
            get_idx_range(field1),
            0.0,
            ddc::reducer::max<double>(),
            KOKKOS_LAMBDA(typename IdxRangeType::discrete_element_type idx) {
                return Kokkos::abs(field1(idx) - field2(idx));
            });
}

/**
 * Check that the CompositeSplineBuilder2D gives the same coefficients as ddc::SplineBuilder2D.
 *
 * @param idx_range_batch The index range of the batch dimensions (if any). These grids must
 *          already be initialised.
 */
template <class Config1, class Config2, class... BatchGrids>
void check_composite_2d_against_ddc(IdxRange<BatchGrids...> idx_range_batch = IdxRange<>())
{
    using Grid1 = typename Config1::Grid;
    using Grid2 = typename Config2::Grid;
    using Builder1 = typename Config1::Builder;
    using Builder2 = typename Config2::Builder;
    using Builder = CompositeSplineBuilder2D<Builder1, Builder2>;
    using RefBuilder = ddc::SplineBuilder2D<
            ExecSpace,
            MemSpace,
            typename Config1::bsplines_type,
            typename Config2::bsplines_type,
            Grid1,
            Grid2,
            Config1::closure,
            Config1::closure,
            Config2::closure,
            Config2::closure,
            ddc::SplineSolver::LAPACK>;

    static_assert(concepts::InterpolationBuilder<Builder>);

    using IdxRangeBatched = IdxRange<BatchGrids..., Grid1, Grid2>;

    init_bsplines();
    IdxRange<Grid1> idx_range_1 = Config1::init_grid();
    IdxRange<Grid2> idx_range_2 = Config2::init_grid();
    IdxRange<Grid1, Grid2> idx_range(idx_range_1, idx_range_2);
    IdxRangeBatched idx_range_batched(idx_range_batch, idx_range);

    Builder1 builder1(idx_range);
    Builder2 builder2(idx_range);
    Builder builder(builder1, builder2);
    RefBuilder ref_builder(idx_range);

    using Dim1 = typename Grid1::continuous_dimension_type;
    using Dim2 = typename Grid2::continuous_dimension_type;
    Coord<Dim1> const min1 = ddc::coordinate(idx_range_1.front());
    Coord<Dim1> const max1 = ddc::coordinate(idx_range_1.back());
    Coord<Dim2> const min2 = ddc::coordinate(idx_range_2.front());
    Coord<Dim2> const max2 = ddc::coordinate(idx_range_2.back());

    DFieldMem<IdxRangeBatched> vals_alloc(idx_range_batched);
    fill_with_test_values<Dim1, Dim2>(get_field(vals_alloc));

    auto idx_range_derivs1 = builder1.batched_derivs_xmin_domain(idx_range_batched);
    DFieldMem<decltype(idx_range_derivs1)> derivs_min1_alloc(idx_range_derivs1);
    DFieldMem<decltype(idx_range_derivs1)> derivs_max1_alloc(idx_range_derivs1);
    fill_with_test_derivs(get_field(derivs_min1_alloc), min1, min2);
    fill_with_test_derivs(get_field(derivs_max1_alloc), max1, min2);

    auto idx_range_derivs2 = builder2.batched_derivs_xmin_domain(idx_range_batched);
    DFieldMem<decltype(idx_range_derivs2)> derivs_min2_alloc(idx_range_derivs2);
    DFieldMem<decltype(idx_range_derivs2)> derivs_max2_alloc(idx_range_derivs2);
    fill_with_test_derivs(get_field(derivs_min2_alloc), min1, min2);
    fill_with_test_derivs(get_field(derivs_max2_alloc), min1, max2);

    auto idx_range_mixed_derivs = builder1.batched_derivs_xmin_domain(idx_range_derivs2);
    DFieldMem<decltype(idx_range_mixed_derivs)> mixed_derivs_min1_min2_alloc(
            idx_range_mixed_derivs);
    DFieldMem<decltype(idx_range_mixed_derivs)> mixed_derivs_max1_min2_alloc(
            idx_range_mixed_derivs);
    DFieldMem<decltype(idx_range_mixed_derivs)> mixed_derivs_min1_max2_alloc(
            idx_range_mixed_derivs);
    DFieldMem<decltype(idx_range_mixed_derivs)> mixed_derivs_max1_max2_alloc(
            idx_range_mixed_derivs);
    fill_with_test_derivs(get_field(mixed_derivs_min1_min2_alloc), min1, min2);
    fill_with_test_derivs(get_field(mixed_derivs_max1_min2_alloc), max1, min2);
    fill_with_test_derivs(get_field(mixed_derivs_min1_max2_alloc), min1, max2);
    fill_with_test_derivs(get_field(mixed_derivs_max1_max2_alloc), max1, max2);

    using IdxRangeCoeffs = typename Builder::template batched_basis_idx_range_type<IdxRangeBatched>;
    DFieldMem<IdxRangeCoeffs> coeffs_alloc(builder.batched_basis_idx_range(idx_range_batched));
    DFieldMem<IdxRangeCoeffs> ref_coeffs_alloc(
            ref_builder.batched_spline_domain(idx_range_batched));

    builder(get_field(coeffs_alloc),
            get_const_field(vals_alloc),
            std::optional(get_const_field(derivs_min1_alloc)),
            std::optional(get_const_field(derivs_max1_alloc)),
            std::optional(get_const_field(derivs_min2_alloc)),
            std::optional(get_const_field(derivs_max2_alloc)),
            std::optional(get_const_field(mixed_derivs_min1_min2_alloc)),
            std::optional(get_const_field(mixed_derivs_max1_min2_alloc)),
            std::optional(get_const_field(mixed_derivs_min1_max2_alloc)),
            std::optional(get_const_field(mixed_derivs_max1_max2_alloc)));
    ref_builder(
            get_field(ref_coeffs_alloc),
            get_const_field(vals_alloc),
            std::optional(get_const_field(derivs_min1_alloc)),
            std::optional(get_const_field(derivs_max1_alloc)),
            std::optional(get_const_field(derivs_min2_alloc)),
            std::optional(get_const_field(derivs_max2_alloc)),
            std::optional(get_const_field(mixed_derivs_min1_min2_alloc)),
            std::optional(get_const_field(mixed_derivs_max1_min2_alloc)),
            std::optional(get_const_field(mixed_derivs_min1_max2_alloc)),
            std::optional(get_const_field(mixed_derivs_max1_max2_alloc)));

    EXPECT_LE(
            max_abs_difference(get_const_field(coeffs_alloc), get_const_field(ref_coeffs_alloc)),
            1e-13);
}

/**
 * Check that the CompositeSplineBuilder3D gives the same coefficients as ddc::SplineBuilder3D.
 *
 * @param idx_range_batch The index range of the batch dimensions (if any). These grids must
 *          already be initialised.
 */
template <class Config1, class Config2, class Config3, class... BatchGrids>
void check_composite_3d_against_ddc(IdxRange<BatchGrids...> idx_range_batch = IdxRange<>())
{
    using Grid1 = typename Config1::Grid;
    using Grid2 = typename Config2::Grid;
    using Grid3 = typename Config3::Grid;
    using Builder1 = typename Config1::Builder;
    using Builder23 = ddc::SplineBuilder2D<
            ExecSpace,
            MemSpace,
            typename Config2::bsplines_type,
            typename Config3::bsplines_type,
            Grid2,
            Grid3,
            Config2::closure,
            Config2::closure,
            Config3::closure,
            Config3::closure,
            ddc::SplineSolver::LAPACK>;
    using Builder = CompositeSplineBuilder3D<Builder1, Builder23>;
    using RefBuilder = ddc::SplineBuilder3D<
            ExecSpace,
            MemSpace,
            typename Config1::bsplines_type,
            typename Config2::bsplines_type,
            typename Config3::bsplines_type,
            Grid1,
            Grid2,
            Grid3,
            Config1::closure,
            Config1::closure,
            Config2::closure,
            Config2::closure,
            Config3::closure,
            Config3::closure,
            ddc::SplineSolver::LAPACK>;

    static_assert(concepts::InterpolationBuilder<Builder>);

    using IdxRangeBatched = IdxRange<BatchGrids..., Grid1, Grid2, Grid3>;

    init_bsplines();
    IdxRange<Grid1> idx_range_1 = Config1::init_grid();
    IdxRange<Grid2> idx_range_2 = Config2::init_grid();
    IdxRange<Grid3> idx_range_3 = Config3::init_grid();
    IdxRange<Grid1, Grid2, Grid3> idx_range(idx_range_1, idx_range_2, idx_range_3);
    IdxRangeBatched idx_range_batched(idx_range_batch, idx_range);

    Builder1 builder1(idx_range);
    Builder23 builder23(idx_range);
    Builder builder(builder1, builder23);
    RefBuilder ref_builder(idx_range);

    using Dim1 = typename Grid1::continuous_dimension_type;
    using Dim2 = typename Grid2::continuous_dimension_type;
    using Dim3 = typename Grid3::continuous_dimension_type;
    Coord<Dim1> const min1 = ddc::coordinate(idx_range_1.front());
    Coord<Dim1> const max1 = ddc::coordinate(idx_range_1.back());
    Coord<Dim2> const min2 = ddc::coordinate(idx_range_2.front());
    Coord<Dim3> const min3 = ddc::coordinate(idx_range_3.front());

    DFieldMem<IdxRangeBatched> vals_alloc(idx_range_batched);
    fill_with_test_values<Dim1, Dim2, Dim3>(get_field(vals_alloc));

    auto idx_range_derivs1 = builder1.batched_derivs_xmin_domain(idx_range_batched);
    DFieldMem<decltype(idx_range_derivs1)> derivs_min1_alloc(idx_range_derivs1);
    DFieldMem<decltype(idx_range_derivs1)> derivs_max1_alloc(idx_range_derivs1);
    fill_with_test_derivs(get_field(derivs_min1_alloc), min1, min2, min3);
    fill_with_test_derivs(get_field(derivs_max1_alloc), max1, min2, min3);

    using IdxRangeCoeffs = typename Builder::template batched_basis_idx_range_type<IdxRangeBatched>;
    DFieldMem<IdxRangeCoeffs> coeffs_alloc(builder.batched_basis_idx_range(idx_range_batched));
    DFieldMem<IdxRangeCoeffs> ref_coeffs_alloc(
            ref_builder.batched_spline_domain(idx_range_batched));

    builder(get_field(coeffs_alloc),
            get_const_field(vals_alloc),
            std::optional(get_const_field(derivs_min1_alloc)),
            std::optional(get_const_field(derivs_max1_alloc)));
    ref_builder(
            get_field(ref_coeffs_alloc),
            get_const_field(vals_alloc),
            std::optional(get_const_field(derivs_min1_alloc)),
            std::optional(get_const_field(derivs_max1_alloc)));

    EXPECT_LE(
            max_abs_difference(get_const_field(coeffs_alloc), get_const_field(ref_coeffs_alloc)),
            1e-13);
}

} // namespace

TEST(CompositeSplineBuilder2D, Types)
{
    using RefTraits = InterpolationBuilderTraits<RefBuilderXY>;
    static_assert(std::is_same_v<BuilderXY::exec_space, RefBuilderXY::exec_space>);
    static_assert(std::is_same_v<BuilderXY::memory_space, RefBuilderXY::memory_space>);
    static_assert(std::is_same_v<BuilderXY::data_type, RefTraits::data_type>);
    static_assert(std::is_same_v<
                  BuilderXY::interpolation_idx_range_type,
                  RefTraits::interpolation_idx_range_type>);
    static_assert(std::is_same_v<BuilderXY::coeff_idx_range_type, RefTraits::coeff_idx_range_type>);
    static_assert(std::is_same_v<
                  BuilderXY::batched_basis_idx_range_type<IdxRangeXBatchY>,
                  RefBuilderXY::batched_spline_domain_type<IdxRangeXBatchY>>);
    static_assert(std::is_same_v<
                  BuilderXY::batched_derivs_idx_range_type1<IdxRangeXBatchY>,
                  RefBuilderXY::batched_derivs_domain_type1<IdxRangeXBatchY>>);
    static_assert(std::is_same_v<
                  BuilderXY::batched_derivs_idx_range_type2<IdxRangeXBatchY>,
                  RefBuilderXY::batched_derivs_domain_type2<IdxRangeXBatchY>>);
    static_assert(std::is_same_v<
                  BuilderXY::batched_derivs_idx_range_type<IdxRangeXBatchY>,
                  RefBuilderXY::batched_derivs_domain_type<IdxRangeXBatchY>>);
}

TEST(CompositeSplineBuilder2D, BatchedBasisIdxRange)
{
    init_discrete_spaces();
    IdxRangeXY idx_range_xy = get_interpolation_idx_range();
    IdxRangeXBatchY idx_range_xby(
            ddc::select<GridX>(idx_range_xy),
            IdxRange<GridBatch>(Idx<GridBatch>(0), IdxStep<GridBatch>(4)),
            ddc::select<GridY>(idx_range_xy));

    BuilderX builder_x(idx_range_xy);
    BuilderY builder_y(idx_range_xy);
    BuilderXY builder(builder_x, builder_y);
    RefBuilderXY ref_builder(idx_range_xy);

    EXPECT_EQ(
            builder.batched_basis_idx_range(idx_range_xy),
            ref_builder.batched_spline_domain(idx_range_xy));
    EXPECT_EQ(
            builder.batched_basis_idx_range(idx_range_xby),
            ref_builder.batched_spline_domain(idx_range_xby));
}

TEST(CompositeSplineBuilder2D, GrevillePeriodic)
{
    check_composite_2d_against_ddc<ConfigXGreville, ConfigYPeriodic>();
}

TEST(CompositeSplineBuilder2D, HermitePeriodic)
{
    check_composite_2d_against_ddc<ConfigXHermite, ConfigYPeriodic>();
}

TEST(CompositeSplineBuilder2D, GrevilleHermite)
{
    check_composite_2d_against_ddc<ConfigXGreville, ConfigZHermite>();
}

TEST(CompositeSplineBuilder2D, HermiteHermite)
{
    check_composite_2d_against_ddc<ConfigXHermite, ConfigZHermite>();
}

TEST(CompositeSplineBuilder2D, HermiteHermiteBatched)
{
    ddc::init_discrete_space<GridBatch>(
            GridBatch::init(Coord<Batch>(0.0), Coord<Batch>(1.0), IdxStep<GridBatch>(4)));
    IdxRange<GridBatch> idx_range_batch(Idx<GridBatch>(0), IdxStep<GridBatch>(4));
    check_composite_2d_against_ddc<ConfigXHermite, ConfigZHermite>(idx_range_batch);
}

TEST(CompositeSplineBuilder3D, Types)
{
    using RefTraits = InterpolationBuilderTraits<RefBuilderXYZ>;
    static_assert(std::is_same_v<BuilderXYZ::exec_space, RefBuilderXYZ::exec_space>);
    static_assert(std::is_same_v<BuilderXYZ::memory_space, RefBuilderXYZ::memory_space>);
    static_assert(std::is_same_v<BuilderXYZ::data_type, RefTraits::data_type>);
    static_assert(std::is_same_v<
                  BuilderXYZ::interpolation_idx_range_type,
                  RefTraits::interpolation_idx_range_type>);
    static_assert(
            std::is_same_v<BuilderXYZ::coeff_idx_range_type, RefTraits::coeff_idx_range_type>);
    static_assert(std::is_same_v<
                  BuilderXYZ::batched_basis_idx_range_type<IdxRangeXBatchYZ>,
                  RefBuilderXYZ::batched_spline_domain_type<IdxRangeXBatchYZ>>);
    static_assert(std::is_same_v<
                  BuilderXYZ::batched_derivs_idx_range_type1<IdxRangeXBatchYZ>,
                  RefBuilderXYZ::batched_derivs_domain_type1<IdxRangeXBatchYZ>>);
    static_assert(std::is_same_v<
                  BuilderXYZ::batched_derivs_idx_range_type2<IdxRangeXBatchYZ>,
                  RefBuilderXYZ::batched_derivs_domain_type2<IdxRangeXBatchYZ>>);
    static_assert(std::is_same_v<
                  BuilderXYZ::batched_derivs_idx_range_type3<IdxRangeXBatchYZ>,
                  RefBuilderXYZ::batched_derivs_domain_type3<IdxRangeXBatchYZ>>);
    static_assert(std::is_same_v<
                  BuilderXYZ::batched_derivs_idx_range_type1_2<IdxRangeXBatchYZ>,
                  RefBuilderXYZ::batched_derivs_domain_type1_2<IdxRangeXBatchYZ>>);
    static_assert(std::is_same_v<
                  BuilderXYZ::batched_derivs_idx_range_type2_3<IdxRangeXBatchYZ>,
                  RefBuilderXYZ::batched_derivs_domain_type2_3<IdxRangeXBatchYZ>>);
    static_assert(std::is_same_v<
                  BuilderXYZ::batched_derivs_idx_range_type1_3<IdxRangeXBatchYZ>,
                  RefBuilderXYZ::batched_derivs_domain_type1_3<IdxRangeXBatchYZ>>);
    static_assert(std::is_same_v<
                  BuilderXYZ::batched_derivs_idx_range_type<IdxRangeXBatchYZ>,
                  RefBuilderXYZ::batched_derivs_domain_type<IdxRangeXBatchYZ>>);
}

TEST(CompositeSplineBuilder3D, BatchedBasisIdxRange)
{
    init_discrete_spaces();
    IdxRangeXYZ idx_range_xyz = get_interpolation_idx_range_3d();
    IdxRangeXBatchYZ idx_range_xbyz(
            ddc::select<GridX>(idx_range_xyz),
            IdxRange<GridBatch>(Idx<GridBatch>(0), IdxStep<GridBatch>(4)),
            ddc::select<GridY, GridZ>(idx_range_xyz));

    BuilderX builder_x(idx_range_xyz);
    BuilderYZ builder_yz(idx_range_xyz);
    BuilderXYZ builder(builder_x, builder_yz);
    RefBuilderXYZ ref_builder(idx_range_xyz);

    EXPECT_EQ(
            builder.batched_basis_idx_range(idx_range_xyz),
            ref_builder.batched_spline_domain(idx_range_xyz));
    EXPECT_EQ(
            builder.batched_basis_idx_range(idx_range_xbyz),
            ref_builder.batched_spline_domain(idx_range_xbyz));
}

TEST(CompositeSplineBuilder3D, GrevillePeriodicGreville)
{
    check_composite_3d_against_ddc<ConfigXGreville, ConfigYPeriodic, ConfigZGreville>();
}

TEST(CompositeSplineBuilder3D, HermitePeriodicGreville)
{
    check_composite_3d_against_ddc<ConfigXHermite, ConfigYPeriodic, ConfigZGreville>();
}
