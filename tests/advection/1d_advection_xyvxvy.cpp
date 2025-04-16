
#include <cmath>

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
// Continuous dimensions
/// @brief A class which describes the real space in the first spatial direction X.
struct X
{
    /// @brief A boolean indicating if the dimension is periodic.
    static bool constexpr PERIODIC = true;
};

/// @brief A class which describes the real space in the second spatial direction Y.
struct Y
{
    /// @brief A boolean indicating if the dimension is periodic.
    static bool constexpr PERIODIC = true;
};

/// @brief A class which describes the real space in the second velocity direction X.
struct Vx
{
};

/// @brief A class which describes the real space in the second velocity direction Y.
struct Vy
{
};


using CoordX = Coord<X>;
using CoordY = Coord<Y>;
using CoordXY = Coord<X, Y>;

using CoordVx = Coord<Vx>;
using CoordVy = Coord<Vy>;

// Splines
struct BSplinesX : ddc::UniformBSplines<X, 3>
{
};
struct BSplinesY : ddc::UniformBSplines<Y, 3>
{
};

ddc::BoundCond constexpr SplineXBoundary = ddc::BoundCond::PERIODIC;
ddc::BoundCond constexpr SplineYBoundary = ddc::BoundCond::PERIODIC;

// Discrete dimensions
struct GridX : UniformGridBase<X>
{
};
struct GridY : UniformGridBase<Y>
{
};
struct GridVx : UniformGridBase<Vx>
{
};
struct GridVy : UniformGridBase<Vy>
{
};


using SplineInterpPointsX
        = ddc::GrevilleInterpolationPoints<BSplinesX, SplineXBoundary, SplineXBoundary>;
using SplineInterpPointsY
        = ddc::GrevilleInterpolationPoints<BSplinesY, SplineYBoundary, SplineYBoundary>;


using IdxXY = Idx<GridX, GridY>;
using IdxVx = Idx<GridVx>;
using IdxVy = Idx<GridVy>;
using IdxXYVxVy = Idx<GridX, GridY, GridVx, GridVy>;


using IdxStepX = IdxStep<GridX>;
using IdxStepY = IdxStep<GridY>;
using IdxStepVx = IdxStep<GridVx>;
using IdxStepVy = IdxStep<GridVy>;


using IdxRangeX = IdxRange<GridX>;
using IdxRangeY = IdxRange<GridY>;
using IdxRangeXY = IdxRange<GridX, GridY>;
using IdxRangeVx = IdxRange<GridVx>;
using IdxRangeVy = IdxRange<GridVy>;
using IdxRangeXYVxVy = IdxRange<GridX, GridY, GridVx, GridVy>;



// Field types
template <class ElementType>
using FieldMemXY = FieldMem<ElementType, IdxRangeXY>;
using DFieldMemXY = FieldMemXY<double>;

template <class ElementType>
using FieldMemXYVxVy = FieldMem<ElementType, IdxRangeXYVxVy>;
using DFieldMemXYVxVy = FieldMemXYVxVy<double>;


template <class ElementType>
using FieldXY = Field<ElementType, IdxRangeXY>;
using DFieldXY = FieldXY<double>;

template <class ElementType>
using FieldXYVxVy = Field<ElementType, IdxRangeXYVxVy>;
using DFieldXYVxVy = FieldXYVxVy<double>;


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
        ddc::PeriodicExtrapolationRule<X>,
        GridX,
        GridY>;

using SplineYBuilder = ddc::SplineBuilder<
        Kokkos::DefaultExecutionSpace,
        Kokkos::DefaultExecutionSpace::memory_space,
        BSplinesY,
        GridY,
        SplineYBoundary,
        SplineYBoundary,
        ddc::SplineSolver::LAPACK>;
using SplineYEvaluator = ddc::SplineEvaluator<
        Kokkos::DefaultExecutionSpace,
        Kokkos::DefaultExecutionSpace::memory_space,
        BSplinesY,
        GridY,
        ddc::PeriodicExtrapolationRule<Y>,
        ddc::PeriodicExtrapolationRule<Y>,
        GridX,
        GridY>;



class XYVxVyAdvection1DTest : public ::testing::Test
{
protected:
    static constexpr CoordX x_min = CoordX(-.5);
    static constexpr CoordX x_max = CoordX(.5);
    static constexpr IdxStepX x_size = IdxStepX(60);

    static constexpr CoordY y_min = CoordY(-.5);
    static constexpr CoordY y_max = CoordY(.5);
    static constexpr IdxStepY y_size = IdxStepY(60);

    static constexpr IdxVx idx0_vx = IdxVx(0);
    static constexpr IdxStepVx vx_size = IdxStepVx(2);

    static constexpr IdxVy idx0_vy = IdxVy(0);
    static constexpr IdxStepVy vy_size = IdxStepVy(2);


    IdxRangeX const interpolation_idx_range_x;
    IdxRangeY const interpolation_idx_range_y;
    IdxRangeXY const xy_grid;

    IdxRangeVx const idx_range_vx;
    IdxRangeVy const idx_range_vy;
    IdxRangeXYVxVy const xyvxvy_grid;

public:
    XYVxVyAdvection1DTest()
        : interpolation_idx_range_x(SplineInterpPointsX::get_domain<GridX>())
        , interpolation_idx_range_y(SplineInterpPointsY::get_domain<GridY>())
        , xy_grid(interpolation_idx_range_x, interpolation_idx_range_y)
        , idx_range_vx(idx0_vx, vx_size)
        , idx_range_vy(idx0_vy, vy_size)
        , xyvxvy_grid(
                  interpolation_idx_range_x,
                  interpolation_idx_range_y,
                  idx_range_vx,
                  idx_range_vy) {};

    ~XYVxVyAdvection1DTest() = default;


    static void SetUpTestSuite()
    {
        ddc::init_discrete_space<BSplinesX>(x_min, x_max, x_size);
        ddc::init_discrete_space<BSplinesY>(y_min, y_max, y_size);

        ddc::init_discrete_space<GridX>(SplineInterpPointsX::get_sampling<GridX>());
        ddc::init_discrete_space<GridY>(SplineInterpPointsY::get_sampling<GridY>());
    }

    template <class AdvectionOperatorX, class AdvectionOperatorY>
    double AdvectionXY(AdvectionOperatorX const& advection_x, AdvectionOperatorY const& advection_y)
    {
        // TIME PARAMETERS ---------------------------------------------------------------------------
        double const dt = 0.05;
        double const final_t = 0.2;
        int const time_iter = int(final_t / dt);

        // INITIALISATION ----------------------------------------------------------------------------
        // Parameters
        double const xc = 0.1;
        double const yc = 0.2;
        double const a = 0.5;

        DFieldMemXYVxVy function_alloc(xyvxvy_grid);
        DFieldXYVxVy function = get_field(function_alloc);
        ddc::parallel_for_each(
                Kokkos::DefaultExecutionSpace(),
                xyvxvy_grid,
                KOKKOS_LAMBDA(IdxXYVxVy const idx) {
                    CoordXY coord_xy = CoordXY(ddc::coordinate(idx));
                    double const x = CoordX(coord_xy);
                    double const y = CoordY(coord_xy);

                    double const r1 = Kokkos::sqrt((x - xc) * (x - xc) + 8 * (y - yc) * (y - yc));
                    double const r2 = Kokkos::sqrt(8 * (x - xc) * (x - xc) + (y - yc) * (y - yc));
                    double const G1
                            = Kokkos::pow(Kokkos::cos(M_PI * r1 / 2. / a), 4) * (abs(r1) < a);
                    double const G2
                            = Kokkos::pow(Kokkos::cos(M_PI * r2 / 2. / a), 4) * (abs(r2) < a);
                    function(idx) = 0.5 * (G1 + G2);
                });

        /*
            The advection field is a Field of doubles instead of a Field 
            of coordinates, because inside the advection operator
            we use an interpolator which needs the DFieldMemXY format.
        */
        DFieldMemXY advection_field_x_alloc(xy_grid);
        DFieldXY advection_field_x = get_field(advection_field_x_alloc);

        DFieldMemXY advection_field_y_alloc(xy_grid);
        DFieldXY advection_field_y = get_field(advection_field_y_alloc);

        ddc::parallel_for_each(
                Kokkos::DefaultExecutionSpace(),
                xy_grid,
                KOKKOS_LAMBDA(IdxXY const idx) {
                    CoordXY coord_xy(ddc::coordinate(idx));
                    double const x = CoordX(coord_xy);
                    advection_field_x(idx) = Kokkos::sin(x * 2 * M_PI) / M_PI / 2.;

                    double const y = CoordY(coord_xy);
                    advection_field_y(idx) = Kokkos::sin(y * 2 * M_PI) / M_PI / 2.;
                });


        // EXACT ADVECTED FUNCTION -------------------------------------------------------------------
        host_t<DFieldMemXYVxVy> exact_function(xyvxvy_grid);
        ddc::for_each(xyvxvy_grid, [&](IdxXYVxVy const idx) {
            CoordXY coord_xy = CoordXY(ddc::coordinate(idx));
            double const x0 = CoordX(coord_xy);
            double const y0 = CoordY(coord_xy);

            double x = 2 * std::atan(std::tan(x0 * M_PI) * std::exp(-final_t)) / M_PI / 2.;
            double y = 2 * std::atan(std::tan(y0 * M_PI) * std::exp(-final_t)) / M_PI / 2.;

            // Replace the feet inside the domain if the dimension is periodic
            if (X::PERIODIC) {
                x = fmod(x - double(x_min), double(x_max - x_min)) + double(x_min);
                x = x > double(x_min) ? x : x + double(x_max - x_min);
            }
            if (Y::PERIODIC) {
                y = fmod(y - double(y_min), double(y_max - y_min)) + double(y_min);
                y = y > double(y_min) ? y : y + double(y_max - y_min);
            }

            double const r1 = std::sqrt((x - xc) * (x - xc) + 8 * (y - yc) * (y - yc));
            double const r2 = std::sqrt(8 * (x - xc) * (x - xc) + (y - yc) * (y - yc));
            double const G1 = std::pow(std::cos(M_PI * r1 / 2. / a), 4) * (abs(r1) < a);
            double const G2 = std::pow(std::cos(M_PI * r2 / 2. / a), 4) * (abs(r2) < a);
            exact_function(idx) = 0.5 * (G1 + G2);
        });


        // SIMULATION --------------------------------------------------------------------------------
        for (int i(0); i < time_iter; i++) {
            advection_x(function, advection_field_x, dt / 2);
            advection_y(function, advection_field_y, dt);
            advection_x(function, advection_field_x, dt / 2);
        };


        // CHECK ERRORS ------------------------------------------------------------------------------
        /*
            Simulation launched on GPU but error checking on CPU. 
        */
        auto function_host = ddc::create_mirror_view_and_copy(function);
        double max_relative_error = 0;
        ddc::for_each(xyvxvy_grid, [&](IdxXYVxVy const idx) {
            double const relative_error = abs(function_host(idx) - exact_function(idx));
            max_relative_error
                    = max_relative_error > relative_error ? max_relative_error : relative_error;
        });
        return max_relative_error;
    };
};

} // end namespace



TEST_F(XYVxVyAdvection1DTest, AdvectionXY)
{
    // CREATING OPERATORS ------------------------------------------------------------------------
    SplineXBuilder_2d const adv_field_builder_x(xy_grid);
    SplineYBuilder_2d const adv_field_builder_y(xy_grid);
    SplineXBuilder_4d const function_builder_x(xyvxvy_grid);
    SplineYBuilder_4d const function_builder_y(xyvxvy_grid);


    ddc::PeriodicExtrapolationRule<X> bv_x_min;
    ddc::PeriodicExtrapolationRule<X> bv_x_max;
    SplineXEvaluator_2d const adv_field_spline_evaluator_x(bv_x_min, bv_x_max);
    SplineXEvaluator_4d const function_spline_evaluator_x(bv_x_min, bv_x_max);

    ddc::PeriodicExtrapolationRule<Y> bv_y_min;
    ddc::PeriodicExtrapolationRule<Y> bv_y_max;
    SplineYEvaluator_2d const adv_field_spline_evaluator_y(bv_y_min, bv_y_max);
    SplineYEvaluator_4d const function_spline_evaluator_y(bv_y_min, bv_y_max);


    PreallocatableSplineInterpolator const
            function_spline_x_interpolator(function_builder_x, function_spline_evaluator_x);
    PreallocatableSplineInterpolator const
            function_spline_y_interpolator(function_builder_y, function_spline_evaluator_y);


    RK2<FieldMemXY<CoordX>, DFieldMemXY> time_stepper_x(xy_grid);
    RK2<FieldMemXY<CoordY>, DFieldMemXY> time_stepper_y(xy_grid);

    BslAdvection1D<
            GridX,
            IdxRangeXY,
            IdxRangeXYVxVy,
            SplineXBuilder_2d,
            SplineXEvaluator_2d,
            RK2<FieldMemXY<CoordX>, DFieldMemXY>> const
            advection_x(
                    function_spline_x_interpolator,
                    adv_field_builder_x,
                    adv_field_spline_evaluator_x,
                    time_stepper_x);
    BslAdvection1D<
            GridY,
            IdxRangeXY,
            IdxRangeXYVxVy,
            SplineYBuilder_2d,
            SplineYEvaluator_2d,
            RK2<FieldMemXY<CoordY>, DFieldMemXY>> const
            advection_y(
                    function_spline_y_interpolator,
                    adv_field_builder_y,
                    adv_field_spline_evaluator_y,
                    time_stepper_y);


    double const max_relative_error = AdvectionXY(advection_x, advection_y);
    EXPECT_LE(max_relative_error, 7.e-2);
    std::cout << "Test on " << x_size << "x" << y_size << "x" << vx_size << "x" << vy_size
              << " grid: max relative error = " << max_relative_error << std::endl;
}
