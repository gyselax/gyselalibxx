/*
    Advection along X on (X). 
*/

#include <ddc/ddc.hpp>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "bsl_advection_1d.hpp"
#include "ddc_helper.hpp"
#include "itimestepper.hpp"
#include "lagrange_basis_uniform.hpp"
#include "lagrange_interpolation.hpp"
#include "rk2.hpp"
#include "spline_interpolation.hpp"
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
using SplineInterpolatorX = SplineInterpolator<
        Kokkos::DefaultExecutionSpace,
        BSplinesX,
        GridX,
        PERIODIC,
        PERIODIC,
        SplineXBoundary,
        SplineXBoundary>;

// Lagrange basis for the advection field interpolation
struct LagBasisX : UniformLagrangeBasis<X, 3, double>
{
};

struct LagBasisFloatX : UniformLagrangeBasis<X, 3, float>
{
};

using LagrangeInterpolatorX = LagrangeInterpolator<
        Kokkos::DefaultExecutionSpace,
        LagBasisX,
        GridX,
        PERIODIC,
        PERIODIC,
        ddc::BoundCond::PERIODIC,
        ddc::BoundCond::PERIODIC>;

using LagrangeInterpolatorFloatX = LagrangeInterpolator<
        Kokkos::DefaultExecutionSpace,
        LagBasisFloatX,
        GridX,
        PERIODIC,
        PERIODIC,
        ddc::BoundCond::PERIODIC,
        ddc::BoundCond::PERIODIC,
        float>;


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
        IdxRangeX interpolation_idx_range(SplineInterpPointsX::get_domain<GridX>());
        IdxRangeX lagrange_break_point_idx_range(
                interpolation_idx_range.front(),
                interpolation_idx_range.extents() + 1);
        ddc::init_discrete_space<LagBasisX>(lagrange_break_point_idx_range);
        ddc::init_discrete_space<LagBasisFloatX>(lagrange_break_point_idx_range);
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
    SplineInterpolatorX spline_interpolation(interpolation_idx_range);

    RK2Builder time_stepper;
    BslAdvection1D<
            GridX,
            IdxRangeX,
            IdxRangeX,
            SplineInterpolatorX,
            SplineInterpolatorX,
            RK2Builder> const advection(spline_interpolation, spline_interpolation, time_stepper);

    double const max_relative_error = AdvectionX(advection);
    EXPECT_LE(max_relative_error, 5.e-3);
    std::cout << "Test on " << x_size << " grid: max relative error = " << max_relative_error
              << std::endl;
}

TEST_F(XAdvection1DTestDouble, AdvectionXLagrange)
{
    // CREATING OPERATORS ------------------------------------------------------------------------
    SplineInterpolatorX spline_interpolation(interpolation_idx_range);
    LagrangeInterpolatorX lag_interpolation;

    RK2Builder time_stepper;
    BslAdvection1D<
            GridX,
            IdxRangeX,
            IdxRangeX,
            SplineInterpolatorX,
            LagrangeInterpolatorX,
            RK2Builder> const advection(spline_interpolation, lag_interpolation, time_stepper);

    double const max_relative_error = AdvectionX(advection);
    EXPECT_LE(max_relative_error, 5.e-2);
    std::cout << "Lagrange test on " << x_size
              << " grid: max relative error = " << max_relative_error << std::endl;
}

TEST_F(XAdvection1DTestFloat, AdvectionX)
{
    // CREATING OPERATORS ------------------------------------------------------------------------
    LagrangeInterpolatorFloatX lag_interpolation;

    RK2Builder time_stepper;
    BslAdvection1D<
            GridX,
            IdxRangeX,
            IdxRangeX,
            LagrangeInterpolatorFloatX,
            LagrangeInterpolatorFloatX,
            RK2Builder,
            float> const advection(lag_interpolation, lag_interpolation, time_stepper);

    double const max_relative_error = AdvectionX(advection);
    EXPECT_LE(max_relative_error, 5.e-3);
    std::cout << "Test on " << x_size << " grid: max relative error = " << max_relative_error
              << std::endl;
}
