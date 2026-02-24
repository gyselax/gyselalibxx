/*
    Advection along X on (X). 
*/

#include <ddc/ddc.hpp>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "bsl_advection_1d.hpp"
#include "ddc_helper.hpp"
#include "identity_interpolation_builder.hpp"
#include "itimestepper.hpp"
#include "lagrange_basis_uniform.hpp"
#include "lagrange_evaluator.hpp"
#include "rk2.hpp"
#include "vector_field_common.hpp"

namespace {
// Continuous dimension
/// @brief A class which describes the real space in the first spatial direction X.
struct X
{
    /// @brief A boolean indicating if the dimension is periodic.
    static bool constexpr PERIODIC = true;
};

using CoordX = Coord<X>;

// Spline
struct BSplinesX : ddc::UniformBSplines<X, 3>
{
};

ddc::BoundCond constexpr SplineXBoundary = ddc::BoundCond::PERIODIC;


// Discrete dimension
struct GridX : UniformGridBase<X>
{
};

using SplineInterpPointsX
        = ddc::GrevilleInterpolationPoints<BSplinesX, SplineXBoundary, SplineXBoundary>;

using IdxRangeX = IdxRange<GridX>;
using IdxX = Idx<GridX>;
using IdxStepX = IdxStep<GridX>;


// Field types
template <class ElementType>
using FieldMemX = FieldMem<ElementType, IdxRangeX>;
using DFieldMemX = FieldMemX<double>;

template <class ElementType>
using FieldX = Field<ElementType, IdxRangeX>;
using DFieldX = FieldX<double>;


// Operators
using SplineXBuilder = ddc::SplineBuilder<
        Kokkos::DefaultExecutionSpace,
        Kokkos::DefaultExecutionSpace::memory_space,
        BSplinesX,
        GridX,
        SplineXBoundary,
        SplineXBoundary,
        ddc::SplineSolver::LAPACK>;

using SplineXEvaluator = ddc::SplineEvaluator<
        Kokkos::DefaultExecutionSpace,
        Kokkos::DefaultExecutionSpace::memory_space,
        BSplinesX,
        GridX,
        ddc::PeriodicExtrapolationRule<X>,
        ddc::PeriodicExtrapolationRule<X>>;

// Lagrange basis for the advection field interpolation
struct LagBasisX : UniformLagrangeBasis<X, 3, double>
{
};

using LagBuilderX = IdentityInterpolationBuilder<
        Kokkos::DefaultExecutionSpace,
        Kokkos::DefaultExecutionSpace::memory_space,
        double,
        GridX,
        LagBasisX>;

using LagEvaluatorX = LagrangeEvaluator<
        Kokkos::DefaultExecutionSpace,
        Kokkos::DefaultExecutionSpace::memory_space,
        double,
        LagBasisX,
        GridX,
        ddc::PeriodicExtrapolationRule<X>,
        ddc::PeriodicExtrapolationRule<X>>;


class XAdvection1DTest : public ::testing::Test
{
protected:
    static constexpr CoordX x_min = CoordX(-M_PI);
    static constexpr CoordX x_max = CoordX(M_PI);
    static constexpr IdxStepX x_size = IdxStepX(32);

    IdxRangeX const interpolation_idx_range;

public:
    XAdvection1DTest() : interpolation_idx_range(SplineInterpPointsX::get_domain<GridX>()) {};

    ~XAdvection1DTest() = default;


    static void SetUpTestSuite()
    {
        ddc::init_discrete_space<BSplinesX>(x_min, x_max, x_size);
        ddc::init_discrete_space<GridX>(SplineInterpPointsX::get_sampling<GridX>());
        IdxRangeX interpolation_idx_range(SplineInterpPointsX::get_domain<GridX>());
        IdxRangeX lagrange_break_point_idx_range(
                interpolation_idx_range.front(),
                interpolation_idx_range.extents() + 1);
        ddc::init_discrete_space<LagBasisX>(lagrange_break_point_idx_range);
    }

    template <class AdvectionOperator>
    double AdvectionX(AdvectionOperator const& advection)
    {
        // TIME PARAMETERS ---------------------------------------------------------------------------
        double const dt = 0.05;
        double const final_t = 0.4;
        int const time_iter = int(final_t / dt);

        // INITIALISATION ----------------------------------------------------------------------------
        DFieldMemX function_alloc(interpolation_idx_range);
        DFieldX function = get_field(function_alloc);

        DFieldMemX advection_field_alloc(interpolation_idx_range);
        DFieldX advection_field = get_field(advection_field_alloc);

        ddc::parallel_for_each(
                Kokkos::DefaultExecutionSpace(),
                interpolation_idx_range,
                KOKKOS_LAMBDA(IdxX const idx) {
                    double const x = ddc::coordinate(idx);
                    function(idx) = Kokkos::sin(2 * x);
                    advection_field(idx) = Kokkos::sin(x);
                });


        // EXACT ADVECTED FUNCTION -------------------------------------------------------------------
        host_t<DFieldMemX> exact_function(interpolation_idx_range);
        ddc::host_for_each(interpolation_idx_range, [&](IdxX const idx) {
            double const x0 = ddc::coordinate(idx);
            double x = 2 * std::atan(std::tan(x0 / 2.) * std::exp(-final_t));
            // Replace the feet inside the domain if the dimension is periodic
            if (X::PERIODIC) {
                x = fmod(x - double(x_min), double(x_max - x_min)) + double(x_min);
                x = x > double(x_min) ? x : x + double(x_max - x_min);
            }
            exact_function(idx) = std::sin(2 * x);
        });


        // SIMULATION --------------------------------------------------------------------------------
        for (int i(0); i < time_iter; i++) {
            advection(function, advection_field, dt);
        };

        // CHECK ERRORS ------------------------------------------------------------------------------
        /*
            Simulation launched on GPU but error checking on CPU. 
        */
        auto function_host = ddc::create_mirror_view_and_copy(function);
        double max_relative_error = 0;
        ddc::host_for_each(interpolation_idx_range, [&](IdxX const idx) {
            double const relative_error = std::abs(function_host(idx) - exact_function(idx));
            max_relative_error
                    = max_relative_error > relative_error ? max_relative_error : relative_error;
        });
        return max_relative_error;
    }
};


} // end namespace



TEST_F(XAdvection1DTest, AdvectionX)
{
    // CREATING OPERATORS ------------------------------------------------------------------------
    SplineXBuilder const builder(interpolation_idx_range);

    ddc::PeriodicExtrapolationRule<X> bv_x_min;
    ddc::PeriodicExtrapolationRule<X> bv_x_max;
    SplineXEvaluator const spline_evaluator(bv_x_min, bv_x_max);

    RK2Builder time_stepper;
    BslAdvection1D<
            GridX,
            IdxRangeX,
            IdxRangeX,
            SplineXBuilder,
            SplineXEvaluator,
            SplineXBuilder,
            SplineXEvaluator,
            RK2Builder> const advection(builder, spline_evaluator, builder, spline_evaluator, time_stepper);

    double const max_relative_error = AdvectionX(advection);
    EXPECT_LE(max_relative_error, 5.e-3);
    std::cout << "Test on " << x_size << " grid: max relative error = " << max_relative_error
              << std::endl;
}

TEST_F(XAdvection1DTest, AdvectionXLagrange)
{
    // CREATING OPERATORS ------------------------------------------------------------------------
    SplineXBuilder const function_builder(interpolation_idx_range);

    ddc::PeriodicExtrapolationRule<X> bv_x_min;
    ddc::PeriodicExtrapolationRule<X> bv_x_max;
    SplineXEvaluator const spline_evaluator(bv_x_min, bv_x_max);

    LagBuilderX const lag_builder;
    LagEvaluatorX const lag_evaluator(bv_x_min, bv_x_max);

    RK2Builder time_stepper;
    BslAdvection1D<
            GridX,
            IdxRangeX,
            IdxRangeX,
            SplineXBuilder,
            SplineXEvaluator,
            LagBuilderX,
            LagEvaluatorX,
            RK2Builder> const advection(function_builder, spline_evaluator, lag_builder, lag_evaluator, time_stepper);

    double const max_relative_error = AdvectionX(advection);
    EXPECT_LE(max_relative_error, 5.e-2);
    std::cout << "Lagrange test on " << x_size
              << " grid: max relative error = " << max_relative_error << std::endl;
}
