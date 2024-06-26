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
struct RDimX
{
    /// @brief A boolean indicating if the dimension is periodic.
    static bool constexpr PERIODIC = true;
};

using CoordX = ddc::Coordinate<RDimX>;

// Spline
struct BSplinesX : ddc::UniformBSplines<RDimX, 3>
{
};

ddc::BoundCond constexpr SplineXBoundary = ddc::BoundCond::PERIODIC;


// Discrete dimension
struct IDimX : ddc::UniformPointSampling<RDimX>
{
};

using SplineInterpPointsX
        = ddc::GrevilleInterpolationPoints<BSplinesX, SplineXBoundary, SplineXBoundary>;

using IDomainX = ddc::DiscreteDomain<IDimX>;
using IndexX = ddc::DiscreteElement<IDimX>;
using IVectX = ddc::DiscreteVector<IDimX>;


// Chunks, Spans and Views
template <class ElementType>
using FieldX = device_t<ddc::Chunk<ElementType, IDomainX>>;
using DFieldX = FieldX<double>;

template <class ElementType>
using SpanX = device_t<ddc::ChunkSpan<ElementType, IDomainX>>;
using DSpanX = SpanX<double>;


// Operators
using SplineXBuilder = ddc::SplineBuilder<
        Kokkos::DefaultExecutionSpace,
        Kokkos::DefaultExecutionSpace::memory_space,
        BSplinesX,
        IDimX,
        SplineXBoundary,
        SplineXBoundary,
        ddc::SplineSolver::GINKGO,
        IDimX>;

using SplineXEvaluator = ddc::SplineEvaluator<
        Kokkos::DefaultExecutionSpace,
        Kokkos::DefaultExecutionSpace::memory_space,
        BSplinesX,
        IDimX,
        ddc::PeriodicExtrapolationRule<RDimX>,
        ddc::PeriodicExtrapolationRule<RDimX>,
        IDimX>;


class XAdvection1DTest : public ::testing::Test
{
protected:
    static constexpr CoordX x_min = CoordX(-M_PI);
    static constexpr CoordX x_max = CoordX(M_PI);
    static constexpr IVectX x_size = IVectX(32);

    IDomainX const interpolation_domain;

public:
    XAdvection1DTest() : interpolation_domain(SplineInterpPointsX::get_domain<IDimX>()) {};

    ~XAdvection1DTest() = default;


    static void SetUpTestSuite()
    {
        ddc::init_discrete_space<BSplinesX>(x_min, x_max, x_size);
        ddc::init_discrete_space<IDimX>(SplineInterpPointsX::get_sampling<IDimX>());
    }

    template <class AdvectionOperator>
    double AdvectionX(AdvectionOperator const& advection)
    {
        // TIME PARAMETERS ---------------------------------------------------------------------------
        double const dt = 0.05;
        double const final_t = 0.4;
        int const time_iter = int(final_t / dt);

        // INITIALISATION ----------------------------------------------------------------------------
        DFieldX function_alloc(interpolation_domain);
        DSpanX function = function_alloc.span_view();

        DFieldX advection_field_alloc(interpolation_domain);
        DSpanX advection_field = advection_field_alloc.span_view();

        ddc::parallel_for_each(
                Kokkos::DefaultExecutionSpace(),
                interpolation_domain,
                KOKKOS_LAMBDA(IndexX const idx) {
                    double const x = ddc::coordinate(idx);
                    function(idx) = Kokkos::sin(2 * x);
                    advection_field(idx) = Kokkos::sin(x);
                });


        // EXACT ADVECTED FUNCTION -------------------------------------------------------------------
        host_t<DFieldX> exact_function(interpolation_domain);
        ddc::for_each(interpolation_domain, [&](IndexX const idx) {
            double const x0 = ddc::coordinate(idx);
            double x = 2 * std::atan(std::tan(x0 / 2.) * std::exp(-final_t));
            // Replace inside the domain the feet if the dimension if periodic
            if (RDimX::PERIODIC) {
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
        auto function_host
                = ddc::create_mirror_view_and_copy(Kokkos::DefaultHostExecutionSpace(), function);
        double max_relative_error = 0;
        ddc::for_each(interpolation_domain, [&](IndexX const idx) {
            double const relative_error = abs(function_host(idx) - exact_function(idx));
            max_relative_error
                    = max_relative_error > relative_error ? max_relative_error : relative_error;
        });
        return max_relative_error;
    };
};


} // end namespace



TEST_F(XAdvection1DTest, AdvectionX)
{
    // CREATING OPERATORS ------------------------------------------------------------------------
    SplineXBuilder const builder(interpolation_domain);

    ddc::PeriodicExtrapolationRule<RDimX> bv_x_min;
    ddc::PeriodicExtrapolationRule<RDimX> bv_x_max;
    SplineXEvaluator const spline_evaluator(bv_x_min, bv_x_max);

    PreallocatableSplineInterpolator const spline_interpolator(builder, spline_evaluator);


    RK2<FieldX<CoordX>, DFieldX> time_stepper(interpolation_domain);
    BslAdvection1D<
            IDimX,
            IDomainX,
            IDomainX,
            SplineXBuilder,
            SplineXEvaluator,
            RK2<FieldX<CoordX>, DFieldX>> const
            advection(spline_interpolator, builder, spline_evaluator, time_stepper);

    double const max_relative_error = AdvectionX(advection);
    EXPECT_LE(max_relative_error, 5.e-3);
    std::cout << "Test on " << x_size << " grid: max relative error = " << max_relative_error
              << std::endl;
}
