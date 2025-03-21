// SPDX-License-Identifier: MIT
#include <ddc/ddc.hpp>
#include <ddc/kernels/splines.hpp>

#include <gtest/gtest.h>

#include "central_fdm_partial_derivatives.hpp"
#include "ddc_aliases.hpp"
#include "mesh_builder.hpp"
#include "spline_1d_partial_derivative.hpp"

namespace {

struct X
{
    static bool constexpr PERIODIC = false;
};
using CoordX = Coord<X>;
struct BSplinesX : ddc::UniformBSplines<X, 3>
{
};
auto constexpr SplineXBoundary = ddc::BoundCond::GREVILLE;
using SplineInterpPointsX
        = ddc::GrevilleInterpolationPoints<BSplinesX, SplineXBoundary, SplineXBoundary>;

struct GridX : NonUniformGridBase<X>
{
};
using IdxX = Idx<GridX>;
using IdxStepX = IdxStep<GridX>;
using IdxRangeX = IdxRange<GridX>;

struct Y
{
    static bool constexpr PERIODIC = false;
};
using CoordY = Coord<Y>;
struct BSplinesY : ddc::UniformBSplines<Y, 3>
{
};
auto constexpr SplineYBoundary = ddc::BoundCond::GREVILLE;
using SplineInterpPointsY
        = ddc::GrevilleInterpolationPoints<BSplinesY, SplineYBoundary, SplineYBoundary>;
struct GridY : SplineInterpPointsY::interpolation_discrete_dimension_type
{
};
using IdxY = Idx<GridY>;
using IdxStepY = IdxStep<GridY>;
using IdxRangeY = IdxRange<GridY>;

using IdxXY = Idx<GridX, GridY>;
using IdxRangeXY = IdxRange<GridX, GridY>;
using DFieldMemXY = DFieldMem<IdxRangeXY>;
using DFieldXY = DField<IdxRangeXY>;
using DConstFieldXY = DConstField<IdxRangeXY>;

// --- Operators ---
using SplineXBuilder = ddc::SplineBuilder<
        Kokkos::DefaultExecutionSpace,
        Kokkos::DefaultExecutionSpace::memory_space,
        BSplinesX,
        GridX,
        SplineXBoundary,
        SplineXBoundary,
        ddc::SplineSolver::LAPACK,
        GridX,
        GridY>;
using SplineXEvaluator = ddc::SplineEvaluator<
        Kokkos::DefaultExecutionSpace,
        Kokkos::DefaultExecutionSpace::memory_space,
        BSplinesX,
        GridX,
        ddc::ConstantExtrapolationRule<X>,
        ddc::ConstantExtrapolationRule<X>,
        GridX,
        GridY>;

template <class Dimension>
void test_partial_derivative_dx(
        IPartialDerivativeCreator<IdxRangeXY, Dimension> const& partial_dx_creator,
        IdxRangeXY const& idxrange_xy)
{
    DFieldMemXY field_xy_alloc(idxrange_xy);
    DFieldXY field_xy = get_field(field_xy_alloc);

    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            idxrange_xy,
            KOKKOS_LAMBDA(IdxXY const idx_xy) {
                IdxX idx_x(idx_xy);
                IdxY idx_y(idx_xy);
                field_xy(idx_xy)
                        = ddc::coordinate(idx_x) * ddc::coordinate(idx_x) * ddc::coordinate(idx_y);
            });


    std::unique_ptr<IPartialDerivative<IdxRangeXY, X>> const partial_dx_pointer
            = partial_dx_creator.create_instance(get_const_field(field_xy));
    IPartialDerivative<IdxRangeXY, X> const& partial_dx = *partial_dx_pointer;

    DFieldMemXY dfield_dx_xy_alloc(idxrange_xy);
    DFieldXY dfield_dx_xy = get_field(dfield_dx_xy_alloc);
    partial_dx(dfield_dx_xy);
    double max_error = ddc::parallel_transform_reduce(
            Kokkos::DefaultExecutionSpace(),
            idxrange_xy,
            0.,
            ddc::reducer::max<double>(),
            KOKKOS_LAMBDA(IdxXY const idx_xy) {
                IdxX idx_x(idx_xy);
                IdxY idx_y(idx_xy);
                double const dfield_dx_anal = 2. * ddc::coordinate(idx_x) * ddc::coordinate(idx_y);
                return Kokkos::abs(dfield_dx_xy(idx_xy) - dfield_dx_anal);
            });
    EXPECT_LE(max_error, 1e-12);
}


TEST(PartialDerivative, Spline1DPartialDerivative)
{
    int n_elems_x(10);
    int n_elems_y(20);

    Coord<X> const x_min(0.0);
    Coord<X> const x_max(1.0);
    IdxStepX x_ncells(n_elems_x);

    Coord<Y> const y_min(0.0);
    Coord<Y> const y_max(2.0);
    IdxStepY y_ncells(n_elems_y);

    ddc::init_discrete_space<BSplinesX>(x_min, x_max, x_ncells);
    ddc::init_discrete_space<GridX>(SplineInterpPointsX::get_sampling<GridX>());
    IdxRangeX idxrange_x(SplineInterpPointsX::get_domain<GridX>());

    ddc::init_discrete_space<BSplinesY>(y_min, y_max, y_ncells);
    ddc::init_discrete_space<GridY>(SplineInterpPointsY::get_sampling<GridY>());
    IdxRangeY idxrange_y(SplineInterpPointsY::get_domain<GridY>());

    IdxRangeXY idxrange_xy(idxrange_x, idxrange_y);

    SplineXBuilder const builder_x(idxrange_xy);
    ddc::ConstantExtrapolationRule<X> bv_x_min(x_min);
    ddc::ConstantExtrapolationRule<X> bv_x_max(x_max);
    SplineXEvaluator const spline_evaluator_x(bv_x_min, bv_x_max);

    Spline1DPartialDerivativeCreator<SplineXBuilder, SplineXEvaluator> const
            partial_dx_creator(builder_x, spline_evaluator_x);

    test_partial_derivative_dx<X>(partial_dx_creator, idxrange_xy);
}


TEST(PartialDerivative, CentralFDMPartialDerivativeDx)
{
    int n_elems_x(10);
    int n_elems_y(20);

    Coord<X> const x_min(0.0);
    Coord<X> const x_max(1.0);
    IdxStepX x_ncells(n_elems_x);

    Coord<Y> const y_min(0.0);
    Coord<Y> const y_max(2.0);
    IdxStepY y_ncells(n_elems_y);

    ddc::init_discrete_space<GridX>(build_random_non_uniform_break_points(x_min, x_max, x_ncells));
    IdxRangeX idxrange_x(IdxX {0}, x_ncells);

    ddc::init_discrete_space<GridY>(build_random_non_uniform_break_points(y_min, y_max, y_ncells));
    IdxRangeY idxrange_y(IdxY {0}, y_ncells);

    IdxRangeXY idxrange_xy(idxrange_x, idxrange_y);

    CentralFDMPartialDerivativeCreator<IdxRangeXY, X> const partial_dx_creator;
    test_partial_derivative_dx(partial_dx_creator, idxrange_xy);
}
} // namespace
