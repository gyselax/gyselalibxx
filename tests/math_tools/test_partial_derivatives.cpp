// SPDX-License-Identifier: MIT
#include <ddc/ddc.hpp>
#include <ddc/kernels/splines.hpp>

#include <gtest/gtest.h>

#include "ddc_aliases.hpp"
#include "math_tools.hpp"
#include "mesh_builder.hpp"
#include "spline_1d_partial_derivative.hpp"
#include "spline_2d_partial_derivative.hpp"

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

using CoordXY = Coord<X, Y>;

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

// --- Operators ---
using SplineYBuilder = ddc::SplineBuilder<
        Kokkos::DefaultExecutionSpace,
        Kokkos::DefaultExecutionSpace::memory_space,
        BSplinesY,
        GridY,
        SplineYBoundary,
        SplineYBoundary,
        ddc::SplineSolver::LAPACK,
        GridX,
        GridY>;
using SplineYEvaluator = ddc::SplineEvaluator<
        Kokkos::DefaultExecutionSpace,
        Kokkos::DefaultExecutionSpace::memory_space,
        BSplinesY,
        GridY,
        ddc::ConstantExtrapolationRule<Y>,
        ddc::ConstantExtrapolationRule<Y>,
        GridX,
        GridY>;

using SplineXYBuilder = ddc::SplineBuilder2D<
        Kokkos::DefaultExecutionSpace,
        Kokkos::DefaultExecutionSpace::memory_space,
        BSplinesX,
        BSplinesY,
        GridX,
        GridY,
        SplineXBoundary,
        SplineXBoundary,
        SplineYBoundary,
        SplineYBoundary,
        ddc::SplineSolver::LAPACK,
        GridX,
        GridY>;

using SplineXYEvaluator = ddc::SplineEvaluator2D<
        Kokkos::DefaultExecutionSpace,
        Kokkos::DefaultExecutionSpace::memory_space,
        BSplinesX,
        BSplinesY,
        GridX,
        GridY,
        ddc::ConstantExtrapolationRule<X, Y>,
        ddc::ConstantExtrapolationRule<X, Y>,
        ddc::ConstantExtrapolationRule<Y, X>,
        ddc::ConstantExtrapolationRule<Y, X>,
        GridX,
        GridY>;


using IdxRangeBSXY = SplineXYBuilder::batched_spline_domain_type;

/**
 * @brief A class that represents a polynomial test function for computing partial derivatives.
 */
class FunctionToDifferentiatePolynomial
{
public:
    /**
     * @brief Get the value of the function at given coordinate.
     *
     * @param[in] coord_xy The coordinate where we want to evaluate
     * the function.
     *
     * @return The value of the function at the coordinate.
     */
    KOKKOS_FUNCTION double operator()(CoordXY const coord_xy) const
    {
        double const x = ddc::get<X>(coord_xy);
        double const y = ddc::get<Y>(coord_xy);
        return (ipow(x, 3) + ipow(y, 2)) * y;
    }

    /**
     * @brief Get the value of the partial derivative of the function
     * in the DerivativeDimension direction at a given coordinate.
     *
     * @param[in] coord_xy The coordinate where we want to evaluate 
     * the partial derivative.
     *
     * @return The value of the partial derivative of the function
     * at the coordinate.
     */
    template <class DerivativeDimension>
    KOKKOS_FUNCTION double differentiate(CoordXY const coord_xy) const
    {
        static_assert(
                std::is_same_v<DerivativeDimension, X> || std::is_same_v<DerivativeDimension, Y>);
        double const x = ddc::get<X>(coord_xy);
        double const y = ddc::get<Y>(coord_xy);

        if constexpr (std::is_same_v<DerivativeDimension, X>) {
            return 3. * ipow(x, 2) * y;
        } else {
            return ipow(x, 3) + 3. * ipow(y, 2);
        }
    }
};


template <class FunctionToDifferentiate, class DerivativeDimension>
void test_partial_derivative(
        FunctionToDifferentiate const function_to_differentiate,
        IPartialDerivativeCreator<IdxRangeXY, DerivativeDimension> const&
                partial_derivative_creator,
        IdxRangeXY const& idxrange_xy)
{
    DFieldMemXY field_xy_alloc(idxrange_xy);
    DFieldXY field_xy = get_field(field_xy_alloc);
    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            idxrange_xy,
            KOKKOS_LAMBDA(IdxXY const idx_xy) {
                field_xy(idx_xy) = function_to_differentiate(ddc::coordinate(idx_xy));
            });


    std::unique_ptr<IPartialDerivative<IdxRangeXY, DerivativeDimension>> const
            partial_derivative_creator_pointer
            = partial_derivative_creator.create_instance(get_const_field(field_xy));
    IPartialDerivative<IdxRangeXY, DerivativeDimension> const& partial_derivative
            = *partial_derivative_creator_pointer;

    DFieldMemXY field_differentiated_alloc(idxrange_xy);
    DFieldXY field_differentiated = get_field(field_differentiated_alloc);
    partial_derivative(field_differentiated);

    double max_error = ddc::parallel_transform_reduce(
            Kokkos::DefaultExecutionSpace(),
            idxrange_xy,
            0.,
            ddc::reducer::max<double>(),
            KOKKOS_LAMBDA(IdxXY const idx_xy) {
                double const field_differentiated_anal
                        = function_to_differentiate.template differentiate<DerivativeDimension>(
                                ddc::coordinate(idx_xy));
                ;
                return Kokkos::abs(field_differentiated(idx_xy) - field_differentiated_anal);
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

    FunctionToDifferentiatePolynomial function_to_differentiate;

    // partial derivatives evaluated with 1d splines in X direction
    SplineXBuilder const builder_x(idxrange_xy);
    ddc::ConstantExtrapolationRule<X> bv_x_min(x_min);
    ddc::ConstantExtrapolationRule<X> bv_x_max(x_max);
    SplineXEvaluator const spline_evaluator_x(bv_x_min, bv_x_max);

    Spline1DPartialDerivativeCreator<SplineXBuilder, SplineXEvaluator> const
            partial_dx_creator(builder_x, spline_evaluator_x);
    test_partial_derivative<
            FunctionToDifferentiatePolynomial,
            X>(function_to_differentiate, partial_dx_creator, idxrange_xy);

    // partial derivatives evaluated with 1d splines in Y direction
    SplineYBuilder const builder_y(idxrange_xy);
    ddc::ConstantExtrapolationRule<Y> bv_y_min(y_min);
    ddc::ConstantExtrapolationRule<Y> bv_y_max(y_max);
    SplineYEvaluator const spline_evaluator_y(bv_y_min, bv_y_max);

    Spline1DPartialDerivativeCreator<SplineYBuilder, SplineYEvaluator> const
            partial_dy_creator(builder_y, spline_evaluator_y);
    test_partial_derivative<
            FunctionToDifferentiatePolynomial,
            Y>(function_to_differentiate, partial_dy_creator, idxrange_xy);


    // partial derivatives evaluated with 2d splines
    SplineXYBuilder builder_xy(idxrange_xy);
    SplineBuilder2DCache<SplineXYBuilder> builder_cache(builder_xy);
    ddc::ConstantExtrapolationRule<X, Y> bv_xy_x_min(x_min, y_min, y_max);
    ddc::ConstantExtrapolationRule<X, Y> bv_xy_x_max(x_max, y_min, y_max);
    ddc::ConstantExtrapolationRule<Y, X> bv_xy_y_min(y_min, x_min, x_max);
    ddc::ConstantExtrapolationRule<Y, X> bv_xy_y_max(y_max, x_min, x_max);
    SplineXYEvaluator evaluator_xy(bv_xy_x_min, bv_xy_x_max, bv_xy_y_min, bv_xy_y_max);

    // derivatives in X direction
    Spline2DPartialDerivativeCreator<
            SplineBuilder2DCache<SplineXYBuilder>,
            SplineXYEvaluator,
            X> const partial2d_dx_creator(builder_cache, evaluator_xy);
    test_partial_derivative<
            FunctionToDifferentiatePolynomial,
            X>(function_to_differentiate, partial2d_dx_creator, idxrange_xy);

    // derivatives in Y direction
    Spline2DPartialDerivativeCreator<
            SplineBuilder2DCache<SplineXYBuilder>,
            SplineXYEvaluator,
            Y> const partial2d_dy_creator(builder_cache, evaluator_xy);
    test_partial_derivative<
            FunctionToDifferentiatePolynomial,
            Y>(function_to_differentiate, partial2d_dy_creator, idxrange_xy);}
} // namespace
