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
struct RDimX
{
    /// @brief A boolean indicating if the dimension is periodic.
    static bool constexpr PERIODIC = true;
};
/// @brief A class which describes the real space in the first velocity direction Vx.
struct RDimVx
{
    /// @brief A boolean indicating if the dimension is periodic.
    static bool constexpr PERIODIC = false;
};


using CoordX = ddc::Coordinate<RDimX>;
using CoordVx = ddc::Coordinate<RDimVx>;

// Spline
struct BSplinesX : ddc::UniformBSplines<RDimX, 3>
{
};
struct BSplinesVx : ddc::UniformBSplines<RDimVx, 3>
{
};

ddc::BoundCond constexpr SplineXBoundary = ddc::BoundCond::PERIODIC;
ddc::BoundCond constexpr SplineVxBoundary = ddc::BoundCond::HERMITE;


// Discrete dimension
struct IDimX : ddc::UniformPointSampling<RDimX>
{
};
struct IDimVx : ddc::UniformPointSampling<RDimVx>
{
};


using SplineInterpPointsX
        = ddc::GrevilleInterpolationPoints<BSplinesX, SplineXBoundary, SplineXBoundary>;
using SplineInterpPointsVx
        = ddc::GrevilleInterpolationPoints<BSplinesVx, SplineVxBoundary, SplineVxBoundary>;


using IDomainX = ddc::DiscreteDomain<IDimX>;
using IndexX = ddc::DiscreteElement<IDimX>;
using IVectX = ddc::DiscreteVector<IDimX>;

using IDomainVx = ddc::DiscreteDomain<IDimVx>;
using IndexVx = ddc::DiscreteElement<IDimVx>;
using IVectVx = ddc::DiscreteVector<IDimVx>;

using IDomainXVx = ddc::DiscreteDomain<IDimX, IDimVx>;
using IndexXVx = ddc::DiscreteElement<IDimX, IDimVx>;


// Chunks, Spans and Views
template <class ElementType>
using FieldX = device_t<ddc::Chunk<ElementType, IDomainX>>;
using DFieldX = FieldX<double>;

template <class ElementType>
using FieldVx = device_t<ddc::Chunk<ElementType, IDomainVx>>;

template <class ElementType>
using FieldXVx = device_t<ddc::Chunk<ElementType, IDomainXVx>>;
using DFieldXVx = FieldXVx<double>;

template <class ElementType>
using SpanX = device_t<ddc::ChunkSpan<ElementType, IDomainX>>;
using DSpanX = SpanX<double>;

template <class ElementType>
using SpanVx = device_t<ddc::ChunkSpan<ElementType, IDomainVx>>;

template <class ElementType>
using SpanXVx = device_t<ddc::ChunkSpan<ElementType, IDomainXVx>>;
using DSpanXVx = SpanXVx<double>;



// Operators
using SplineXBuilder = ddc::SplineBuilder<
        Kokkos::DefaultExecutionSpace,
        Kokkos::DefaultExecutionSpace::memory_space,
        BSplinesX,
        IDimX,
        SplineXBoundary,
        SplineXBoundary,
        ddc::SplineSolver::LAPACK,
        IDimX,
        IDimVx>;
using SplineXEvaluator = ddc::SplineEvaluator<
        Kokkos::DefaultExecutionSpace,
        Kokkos::DefaultExecutionSpace::memory_space,
        BSplinesX,
        IDimX,
        ddc::PeriodicExtrapolationRule<RDimX>,
        ddc::PeriodicExtrapolationRule<RDimX>,
        IDimX,
        IDimVx>;


class XVxAdvection1DTest : public ::testing::Test
{
protected:
    static constexpr CoordX x_min = CoordX(-M_PI);
    static constexpr CoordX x_max = CoordX(M_PI);
    static constexpr IVectX x_size = IVectX(100);

    static constexpr CoordVx vx_min = CoordVx(-6);
    static constexpr CoordVx vx_max = CoordVx(6);
    static constexpr IVectVx vx_size = IVectVx(50);

    IDomainX const x_dom;
    IDomainVx const vx_dom;
    IDomainXVx const xvx_dom;

public:
    XVxAdvection1DTest()
        : x_dom(SplineInterpPointsX::get_domain<IDimX>())
        , vx_dom(SplineInterpPointsVx::get_domain<IDimVx>())
        , xvx_dom(x_dom, vx_dom) {};

    ~XVxAdvection1DTest() = default;


    static void SetUpTestSuite()
    {
        ddc::init_discrete_space<BSplinesX>(x_min, x_max, x_size);
        ddc::init_discrete_space<BSplinesVx>(vx_min, vx_max, vx_size);

        ddc::init_discrete_space<IDimX>(SplineInterpPointsX::get_sampling<IDimX>());
        ddc::init_discrete_space<IDimVx>(SplineInterpPointsVx::get_sampling<IDimVx>());
    }

    template <class AdvectionOperator>
    double AdvectionXVx(AdvectionOperator const& advection)
    {
        // TIME PARAMETERS ---------------------------------------------------------------------------
        double const dt = 0.1;
        double const final_t = 0.4;
        int const time_iter = int(final_t / dt);

        // INITIALISATION ----------------------------------------------------------------------------
        DFieldXVx function_alloc(xvx_dom);
        DSpanXVx function = function_alloc.span_view();

        DFieldXVx advection_field_alloc(xvx_dom);
        DSpanXVx advection_field = advection_field_alloc.span_view();

        ddc::parallel_for_each(
                Kokkos::DefaultExecutionSpace(),
                xvx_dom,
                KOKKOS_LAMBDA(IndexXVx const xv_idx) {
                    double const x = ddc::coordinate(IndexX(xv_idx));
                    function(xv_idx) = Kokkos::cos(x);

                    double const v = ddc::coordinate(IndexVx(xv_idx));
                    advection_field(xv_idx) = v;
                });


        // EXACT ADVECTED FUNCTION -------------------------------------------------------------------
        host_t<DFieldXVx> exact_function(xvx_dom);
        ddc::for_each(xvx_dom, [&](IndexXVx const xv_idx) {
            double const x0 = ddc::coordinate(IndexX(xv_idx));
            double const v = ddc::coordinate(IndexVx(xv_idx));
            double x = x0 - final_t * v;
            // Replace inside the domain the feet if the dimension if periodic
            if (RDimX::PERIODIC) {
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
        ddc::for_each(xvx_dom, [&](IndexXVx const xv_idx) {
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

    ddc::PeriodicExtrapolationRule<RDimX> bv_x_min;
    ddc::PeriodicExtrapolationRule<RDimX> bv_x_max;
    SplineXEvaluator const spline_evaluator_x(bv_x_min, bv_x_max);

    PreallocatableSplineInterpolator const spline_interpolator_x(builder_x, spline_evaluator_x);

    RK2<FieldXVx<CoordX>, DFieldXVx> time_stepper(xvx_dom);
    BslAdvection1D<
            IDimX,
            IDomainXVx,
            IDomainXVx,
            SplineXBuilder,
            SplineXEvaluator,
            RK2<FieldXVx<CoordX>, DFieldXVx>> const
            advection(spline_interpolator_x, builder_x, spline_evaluator_x, time_stepper);

    double const max_relative_error = AdvectionXVx(advection);
    EXPECT_LE(max_relative_error, 5.e-7);
    std::cout << "Test on " << x_size << "x" << vx_size
              << " grid: max relative error = " << max_relative_error << std::endl;
}
