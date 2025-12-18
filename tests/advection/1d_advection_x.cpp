/*
    Advection along X on (X). 
*/

#include <ddc/ddc.hpp>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "bsl_advection_1d.hpp"
#include "ddc_helper.hpp"
#include "itimestepper.hpp"
#include "rk2.hpp"
#include "spline_interpolator.hpp"
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

template <class ElementType>
using FieldX = Field<ElementType, IdxRangeX>;


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


template <class DataType>
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
    }

    template <class AdvectionOperator>
    double AdvectionX(AdvectionOperator const& advection)
    {
        // TIME PARAMETERS ---------------------------------------------------------------------------
        DataType const dt = 0.05;
        DataType const final_t = 0.4;
        int const time_iter = int(final_t / dt);

        // INITIALISATION ----------------------------------------------------------------------------
        FieldMemX<DataType> function_alloc(interpolation_idx_range);
        FieldX<DataType> function = get_field(function_alloc);

        FieldMemX<DataType> advection_field_alloc(interpolation_idx_range);
        FieldX<DataType> advection_field = get_field(advection_field_alloc);

        ddc::parallel_for_each(
                Kokkos::DefaultExecutionSpace(),
                interpolation_idx_range,
                KOKKOS_LAMBDA(IdxX const idx) {
                    DataType const x = ddc::coordinate(idx);
                    function(idx) = Kokkos::sin(2 * x);
                    advection_field(idx) = Kokkos::sin(x);
                });


        // EXACT ADVECTED FUNCTION -------------------------------------------------------------------
        host_t<FieldMemX<DataType>> exact_function(interpolation_idx_range);
        ddc::host_for_each(interpolation_idx_range, [&](IdxX const idx) {
            DataType const x0 = ddc::coordinate(idx);
            DataType x = 2 * std::atan(std::tan(x0 / 2.) * std::exp(-final_t));
            // Replace the feet inside the domain if the dimension is periodic
            if (X::PERIODIC) {
                x = std::fmod(x - DataType(x_min), DataType(x_max - x_min)) + DataType(x_min);
                x = x > DataType(x_min) ? x : x + DataType(x_max - x_min);
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

class XAdvection1DTestDouble : public XAdvection1DTest<double>
{
};

class XAdvection1DTestFloat : public XAdvection1DTest<float>
{
};

} // end namespace



TEST_F(XAdvection1DTestDouble, AdvectionX)
{
    // CREATING OPERATORS ------------------------------------------------------------------------
    SplineXBuilder const builder(interpolation_idx_range);

    ddc::PeriodicExtrapolationRule<X> bv_x_min;
    ddc::PeriodicExtrapolationRule<X> bv_x_max;
    SplineXEvaluator const spline_evaluator(bv_x_min, bv_x_max);

    PreallocatableSplineInterpolator const
            spline_interpolator(builder, spline_evaluator, interpolation_idx_range);


    RK2Builder time_stepper;
    BslAdvection1D<GridX, IdxRangeX, IdxRangeX, SplineXBuilder, SplineXEvaluator, RK2Builder> const
            advection(spline_interpolator, builder, spline_evaluator, time_stepper);

    double const max_relative_error = AdvectionX(advection);
    EXPECT_LE(max_relative_error, 5.e-3);
    std::cout << "Test on " << x_size << " grid: max relative error = " << max_relative_error
              << std::endl;
}

TEST_F(XAdvection1DTestFloat, AdvectionX)
{
    // CREATING OPERATORS ------------------------------------------------------------------------
    // TODO: Use a multi-precision interpolator (E.g. Lagrange)
    SplineXBuilder const builder(interpolation_idx_range);

    ddc::PeriodicExtrapolationRule<X> bv_x_min;
    ddc::PeriodicExtrapolationRule<X> bv_x_max;
    SplineXEvaluator const spline_evaluator(bv_x_min, bv_x_max);

    PreallocatableSplineInterpolator const
            spline_interpolator(builder, spline_evaluator, interpolation_idx_range);


    RK2Builder time_stepper;
    BslAdvection1D<
            GridX,
            IdxRangeX,
            IdxRangeX,
            SplineXBuilder,
            SplineXEvaluator,
            RK2Builder,
            float> const advection(spline_interpolator, builder, spline_evaluator, time_stepper);

    double const max_relative_error = AdvectionX(advection);
    EXPECT_LE(max_relative_error, 5.e-3);
    std::cout << "Test on " << x_size << " grid: max relative error = " << max_relative_error
              << std::endl;
}
