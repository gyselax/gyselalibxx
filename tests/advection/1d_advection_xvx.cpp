/*
    Advection along X on (X, Vx).  The advection field is given by Vx.
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
/// @brief A class which describes the real space in the first velocity direction Vx.
struct Vx
{
    /// @brief A boolean indicating if the dimension is periodic.
    static bool constexpr PERIODIC = false;
};


using CoordX = Coord<X>;
using CoordVx = Coord<Vx>;

// Spline
struct BSplinesX : ddc::UniformBSplines<X, 3>
{
};
struct BSplinesVx : ddc::UniformBSplines<Vx, 3>
{
};

ddc::BoundCond constexpr SplineXBoundary = ddc::BoundCond::PERIODIC;
ddc::BoundCond constexpr SplineVxBoundary = ddc::BoundCond::HERMITE;


// Discrete dimension
struct GridX : UniformGridBase<X>
{
};
struct GridVx : UniformGridBase<Vx>
{
};


using SplineInterpPointsX
        = ddc::GrevilleInterpolationPoints<BSplinesX, SplineXBoundary, SplineXBoundary>;
using SplineInterpPointsVx
        = ddc::GrevilleInterpolationPoints<BSplinesVx, SplineVxBoundary, SplineVxBoundary>;


using IdxRangeX = IdxRange<GridX>;
using IdxX = Idx<GridX>;
using IdxStepX = IdxStep<GridX>;

using IdxRangeVx = IdxRange<GridVx>;
using IdxVx = Idx<GridVx>;
using IdxStepVx = IdxStep<GridVx>;

using IdxRangeXVx = IdxRange<GridX, GridVx>;
using IdxXVx = Idx<GridX, GridVx>;


// Field types
template <class ElementType>
using FieldMemX = FieldMem<ElementType, IdxRangeX>;
using DFieldMemX = FieldMemX<double>;

template <class ElementType>
using FieldMemVx = FieldMem<ElementType, IdxRangeVx>;

template <class ElementType>
using FieldMemXVx = FieldMem<ElementType, IdxRangeXVx>;
using DFieldMemXVx = FieldMemXVx<double>;

template <class ElementType>
using FieldX = Field<ElementType, IdxRangeX>;
using DFieldX = FieldX<double>;

template <class ElementType>
using FieldVx = Field<ElementType, IdxRangeVx>;

template <class ElementType>
using FieldXVx = Field<ElementType, IdxRangeXVx>;
using DFieldXVx = FieldXVx<double>;



// Operators
using SplineXBuilder = ddc::SplineBuilder<
        Kokkos::DefaultExecutionSpace,
        Kokkos::DefaultExecutionSpace::memory_space,
        BSplinesX,
        GridX,
        SplineXBoundary,
        SplineXBoundary,
        ddc::SplineSolver::LAPACK,
        GridX,
        GridVx>;
using SplineXEvaluator = ddc::SplineEvaluator<
        Kokkos::DefaultExecutionSpace,
        Kokkos::DefaultExecutionSpace::memory_space,
        BSplinesX,
        GridX,
        ddc::PeriodicExtrapolationRule<X>,
        ddc::PeriodicExtrapolationRule<X>,
        GridX,
        GridVx>;


class XVxAdvection1DTest : public ::testing::Test
{
protected:
    static constexpr CoordX x_min = CoordX(-M_PI);
    static constexpr CoordX x_max = CoordX(M_PI);
    static constexpr IdxStepX x_size = IdxStepX(100);

    static constexpr CoordVx vx_min = CoordVx(-6);
    static constexpr CoordVx vx_max = CoordVx(6);
    static constexpr IdxStepVx vx_size = IdxStepVx(50);

    IdxRangeX const x_dom;
    IdxRangeVx const vx_dom;
    IdxRangeXVx const xvx_dom;

public:
    XVxAdvection1DTest()
        : x_dom(SplineInterpPointsX::get_domain<GridX>())
        , vx_dom(SplineInterpPointsVx::get_domain<GridVx>())
        , xvx_dom(x_dom, vx_dom) {};

    ~XVxAdvection1DTest() = default;


    static void SetUpTestSuite()
    {
        ddc::init_discrete_space<BSplinesX>(x_min, x_max, x_size);
        ddc::init_discrete_space<BSplinesVx>(vx_min, vx_max, vx_size);

        ddc::init_discrete_space<GridX>(SplineInterpPointsX::get_sampling<GridX>());
        ddc::init_discrete_space<GridVx>(SplineInterpPointsVx::get_sampling<GridVx>());
    }

    template <class AdvectionOperator>
    double AdvectionXVx(AdvectionOperator const& advection)
    {
        // TIME PARAMETERS ---------------------------------------------------------------------------
        double const dt = 0.1;
        double const final_t = 0.4;
        int const time_iter = int(final_t / dt);

        // INITIALISATION ----------------------------------------------------------------------------
        DFieldMemXVx function_alloc(xvx_dom);
        DFieldXVx function = get_field(function_alloc);

        DFieldMemXVx advection_field_alloc(xvx_dom);
        DFieldXVx advection_field = get_field(advection_field_alloc);

        ddc::parallel_for_each(
                Kokkos::DefaultExecutionSpace(),
                xvx_dom,
                KOKKOS_LAMBDA(IdxXVx const xv_idx) {
                    double const x = ddc::coordinate(IdxX(xv_idx));
                    function(xv_idx) = Kokkos::cos(x);

                    double const v = ddc::coordinate(IdxVx(xv_idx));
                    advection_field(xv_idx) = v;
                });


        // EXACT ADVECTED FUNCTION -------------------------------------------------------------------
        host_t<DFieldMemXVx> exact_function(xvx_dom);
        ddc::for_each(xvx_dom, [&](IdxXVx const xv_idx) {
            double const x0 = ddc::coordinate(IdxX(xv_idx));
            double const v = ddc::coordinate(IdxVx(xv_idx));
            double x = x0 - final_t * v;
            // Replace inside the index range the feet if the dimension if periodic
            if (X::PERIODIC) {
                x = fmod(x - double(x_min), double(x_max - x_min)) + double(x_min);
                x = x > double(x_min) ? x : x + double(x_max - x_min);
            }
            exact_function(xv_idx) = Kokkos::cos(x);
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
        ddc::for_each(xvx_dom, [&](IdxXVx const xv_idx) {
            double const relative_error = abs(function_host(xv_idx) - exact_function(xv_idx));
            EXPECT_LE(relative_error, 5.e-7);
            max_relative_error
                    = max_relative_error > relative_error ? max_relative_error : relative_error;
        });
        return max_relative_error;
    };
};


} // end namespace



TEST_F(XVxAdvection1DTest, AdvectionXVx)
{
    // CREATING OPERATORS ------------------------------------------------------------------------
    SplineXBuilder const builder_x(xvx_dom);

    ddc::PeriodicExtrapolationRule<X> bv_x_min;
    ddc::PeriodicExtrapolationRule<X> bv_x_max;
    SplineXEvaluator const spline_evaluator_x(bv_x_min, bv_x_max);

    PreallocatableSplineInterpolator const spline_interpolator_x(builder_x, spline_evaluator_x);

    RK2<FieldMemXVx<CoordX>, DFieldMemXVx> time_stepper(xvx_dom);
    BslAdvection1D<
            GridX,
            IdxRangeXVx,
            IdxRangeXVx,
            SplineXBuilder,
            SplineXEvaluator,
            RK2<FieldMemXVx<CoordX>, DFieldMemXVx>> const
            advection(spline_interpolator_x, builder_x, spline_evaluator_x, time_stepper);

    double const max_relative_error = AdvectionXVx(advection);
    EXPECT_LE(max_relative_error, 5.e-7);
    std::cout << "Test on " << x_size << "x" << vx_size
              << " grid: max relative error = " << max_relative_error << std::endl;
}
