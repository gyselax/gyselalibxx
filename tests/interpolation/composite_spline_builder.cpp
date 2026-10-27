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

auto constexpr ClosureX = ddc::SplineBuilderClosure::GREVILLE;
auto constexpr ClosureY = ddc::SplineBuilderClosure::PERIODIC;

using SplineInterpPointsX = ddc::GrevilleInterpolationPoints<BSplinesX, ClosureX, ClosureX>;
using SplineInterpPointsY = ddc::GrevilleInterpolationPoints<BSplinesY, ClosureY, ClosureY>;

struct GridX : SplineInterpPointsX::interpolation_discrete_dimension_type
{
};
struct GridY : SplineInterpPointsY::interpolation_discrete_dimension_type
{
};
struct GridBatch : UniformGridBase<Batch>
{
};

using IdxRangeXY = IdxRange<GridX, GridY>;
using IdxRangeXBatchY = IdxRange<GridX, GridBatch, GridY>;

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

void init_discrete_spaces()
{
    ddc::init_discrete_space<BSplinesX>(Coord<X>(0.0), Coord<X>(1.0), IdxStep<GridX>(10));
    ddc::init_discrete_space<BSplinesY>(Coord<Y>(0.0), Coord<Y>(2.0), IdxStep<GridY>(12));
    ddc::init_discrete_space<GridX>(SplineInterpPointsX::get_sampling<GridX>());
    ddc::init_discrete_space<GridY>(SplineInterpPointsY::get_sampling<GridY>());
    ddc::init_discrete_space<GridBatch>(
            GridBatch::init(Coord<Batch>(0.0), Coord<Batch>(1.0), IdxStep<GridBatch>(4)));
}

IdxRangeXY get_interpolation_idx_range()
{
    return IdxRangeXY(
            SplineInterpPointsX::get_domain<GridX>(),
            SplineInterpPointsY::get_domain<GridY>());
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

void init_bsplines()
{
    ddc::init_discrete_space<BSplinesX>(Coord<X>(0.0), Coord<X>(1.0), IdxStep<GridX>(10));
    ddc::init_discrete_space<BSplinesY>(Coord<Y>(0.0), Coord<Y>(2.0), IdxStep<GridY>(12));
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
 */
template <class... Grids>
void fill_with_test_values(DField<IdxRange<Grids...>> vals)
{
    using Grid1 = ddc::type_seq_element_t<sizeof...(Grids)-2, ddc::detail::TypeSeq<Grids...>>;
    using Grid2 = ddc::type_seq_element_t<sizeof...(Grids)-1, ddc::detail::TypeSeq<Grids...>>;
    ddc::parallel_for_each(
            ExecSpace(),
            get_idx_range(vals),
            KOKKOS_LAMBDA(Idx<Grids...> idx) {
                vals(idx) = test_function_1d(ddc::coordinate(Idx<Grid1>(idx)))
                            * test_function_1d(ddc::coordinate(Idx<Grid2>(idx)));
            });
}

/**
 * Fill a field with the derivatives of the test function along DimDeriv on the boundary
 * of that dimension. The test function is constant along the batch dimensions.
 */
template <class DimDeriv, class... Grids>
void fill_with_test_derivs(DField<IdxRange<Grids...>> derivs, Coord<DimDeriv> boundary)
{
    using GridDeriv = ddc::type_seq_element_t<sizeof...(Grids)-2, ddc::detail::TypeSeq<Grids...>>;
    using Grid2 = ddc::type_seq_element_t<sizeof...(Grids)-1, ddc::detail::TypeSeq<Grids...>>;
    static_assert(std::is_same_v<GridDeriv, ddc::Deriv<DimDeriv>>);
    ddc::parallel_for_each(
            ExecSpace(),
            get_idx_range(derivs),
            KOKKOS_LAMBDA(Idx<Grids...> idx) {
                int const deriv_order = (Idx<ddc::Deriv<DimDeriv>>(idx)
                                         - Idx<ddc::Deriv<DimDeriv>>(0))
                                                .value();
                derivs(idx) = test_function_1d(boundary, deriv_order)
                              * test_function_1d(ddc::coordinate(Idx<Grid2>(idx)));
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

/// Check that the CompositeSplineBuilder2D gives the same coefficients as ddc::SplineBuilder2D.
template <class Config1, class Config2>
void check_composite_2d_against_ddc()
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

    init_bsplines();
    IdxRange<Grid1, Grid2> idx_range(Config1::init_grid(), Config2::init_grid());

    Builder1 builder1(idx_range);
    Builder2 builder2(idx_range);
    Builder builder(builder1, builder2);
    RefBuilder ref_builder(idx_range);

    DFieldMem<IdxRange<Grid1, Grid2>> vals_alloc(idx_range);
    fill_with_test_values(get_field(vals_alloc));

    auto idx_range_derivs1 = builder1.batched_derivs_xmin_domain(idx_range);
    DFieldMem<decltype(idx_range_derivs1)> derivs_min1_alloc(idx_range_derivs1);
    DFieldMem<decltype(idx_range_derivs1)> derivs_max1_alloc(idx_range_derivs1);
    IdxRange<Grid1> idx_range_1(idx_range);
    fill_with_test_derivs(get_field(derivs_min1_alloc), ddc::coordinate(idx_range_1.front()));
    fill_with_test_derivs(get_field(derivs_max1_alloc), ddc::coordinate(idx_range_1.back()));

    using IdxRangeCoeffs =
            typename Builder::template batched_basis_idx_range_type<IdxRange<Grid1, Grid2>>;
    DFieldMem<IdxRangeCoeffs> coeffs_alloc(builder.batched_basis_idx_range(idx_range));
    DFieldMem<IdxRangeCoeffs> ref_coeffs_alloc(ref_builder.batched_spline_domain(idx_range));

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
