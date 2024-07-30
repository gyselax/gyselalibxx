/*
    Advection along X on (Sp, X, Vx). The advection field is given by Vx. 
    The test is the same as the one applied to BslAdvectionSpatial operator.
*/

#include <ddc/ddc.hpp>

#include <gtest/gtest.h>

#include "bsl_advection_1d.hpp"
#include "euler.hpp"
#include "species_info.hpp"
#include "spline_interpolator.hpp"


namespace {
// Continuous dimensions
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

// Splines
struct BSplinesX : ddc::UniformBSplines<RDimX, 3>
{
};
struct BSplinesVx : ddc::UniformBSplines<RDimVx, 3>
{
};

ddc::BoundCond constexpr SplineXBoundary = ddc::BoundCond::PERIODIC;
ddc::BoundCond constexpr SplineVxBoundary = ddc::BoundCond::HERMITE;


// Discrete dimensions
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

using IDomainSp = ddc::DiscreteDomain<IDimSp>;
using IndexSp = ddc::DiscreteElement<IDimSp>;
using IVectSp = ddc::DiscreteVector<IDimSp>;

using IDomainSpXVx = ddc::DiscreteDomain<IDimSp, IDimX, IDimVx>;
using IndexSpXVx = ddc::DiscreteElement<IDimSp, IDimX, IDimVx>;

// Chunks, Spans and Views
template <class ElementType>
using FieldSp = device_t<ddc::Chunk<ElementType, IDomainSp>>;
using DFieldSp = DFieldSp;

template <class ElementType>
using FieldSpXVx = device_t<ddc::Chunk<ElementType, IDomainSpXVx>>;
using DFieldSpXVx = FieldSpXVx<double>;

template <class ElementType>
using SpanSp = device_t<ddc::ChunkSpan<ElementType, IDomainSp>>;
using DSpanSp = SpanSp<double>;

template <class ElementType>
using SpanSpXVx = device_t<ddc::ChunkSpan<ElementType, IDomainSpXVx>>;
using DSpanSpXVx = SpanSpXVx<double>;


template <class ElementType>
using FieldXVx = device_t<ddc::Chunk<ElementType, IDomainXVx>>;

template <class ElementType>
using SpanXVx = device_t<ddc::ChunkSpan<ElementType, IDomainXVx>>;



// Operators
using SplineXBuilder = ddc::SplineBuilder<
        Kokkos::DefaultExecutionSpace,
        Kokkos::DefaultExecutionSpace::memory_space,
        BSplinesX,
        IDimX,
        SplineXBoundary,
        SplineXBoundary,
        ddc::SplineSolver::LAPACK,
        IDimSp,
        IDimX,
        IDimVx>;
using SplineXEvaluator = ddc::SplineEvaluator<
        Kokkos::DefaultExecutionSpace,
        Kokkos::DefaultExecutionSpace::memory_space,
        BSplinesX,
        IDimX,
        ddc::PeriodicExtrapolationRule<RDimX>,
        ddc::PeriodicExtrapolationRule<RDimX>,
        IDimSp,
        IDimX,
        IDimVx>;


class Spatial1DAdvectionTest : public ::testing::Test
{
protected:
    IDomainX const x_dom;
    IDomainVx const vx_dom;
    IDomainSp const dom_allsp;

    static constexpr IVectSp nb_species = IVectSp(2);

public:
    Spatial1DAdvectionTest()
        : x_dom(SplineInterpPointsX::get_domain<IDimX>())
        , vx_dom(SplineInterpPointsVx::get_domain<IDimVx>())
        , dom_allsp(IndexSp(0), nb_species) {};

    ~Spatial1DAdvectionTest() = default;


    static void SetUpTestSuite()
    {
        CoordX const x_min(-M_PI);
        CoordX const x_max(M_PI);
        IVectX const x_size(100);

        CoordVx const vx_min(-6);
        CoordVx const vx_max(6);
        IVectVx const vx_size(50);

        ddc::init_discrete_space<BSplinesX>(x_min, x_max, x_size);
        ddc::init_discrete_space<BSplinesVx>(vx_min, vx_max, vx_size);

        ddc::init_discrete_space<IDimX>(SplineInterpPointsX::get_sampling<IDimX>());
        ddc::init_discrete_space<IDimVx>(SplineInterpPointsVx::get_sampling<IDimVx>());
    }

    template <class AdvectionOperator>
    double SpatialAdvection(AdvectionOperator const& advection_x)
    {
        // Mesh ----------------------------------------------------------------------------------
        IDomainSpXVx const meshSpXVx(dom_allsp, x_dom, vx_dom);
        IndexSp const i_elec = dom_allsp.front();
        IndexSp const i_ion = dom_allsp.back();


        // INITIALISATION ------------------------------------------------------------------------
        // Initialization of the masses
        host_t<FieldSp<int>> masses_host(dom_allsp);
        masses_host(i_elec) = 1;
        masses_host(i_ion) = 1;
        auto masses_alloc = ddc::
                create_mirror_and_copy(Kokkos::DefaultExecutionSpace(), masses_host.span_view());
        ddc::ChunkSpan masses = masses_alloc.span_view();


        // Initialization of the distribution function and advection field
        DFieldSpXVx allfdistribu_alloc(meshSpXVx);
        DSpanSpXVx allfdistribu = allfdistribu_alloc.span_view();

        DFieldSpXVx advection_field_alloc(meshSpXVx);
        DSpanSpXVx advection_field = advection_field_alloc.span_view();

        ddc::parallel_for_each(
                Kokkos::DefaultExecutionSpace(),
                meshSpXVx,
                KOKKOS_LAMBDA(IndexSpXVx const ispxvx) {
                    IndexX const ix(ispxvx);
                    allfdistribu(ispxvx) = Kokkos::cos(ddc::coordinate(ix));

                    IndexVx const ivx(ispxvx);
                    IndexSp const isp(ispxvx);
                    advection_field(ispxvx)
                            = -Kokkos::sqrt(masses(i_elec) / masses(isp)) * ddc::coordinate(ivx);
                });


        // SIMULATION ----------------------------------------------------------------------------
        double const timestep = .1;

        advection_x(allfdistribu, advection_field, timestep);


        double const max_advection_error = ddc::parallel_transform_reduce(
                Kokkos::DefaultExecutionSpace(),
                meshSpXVx,
                0.0,
                ddc::reducer::max<double>(),
                KOKKOS_LAMBDA(IndexSpXVx const ispxvx) {
                    IndexX const ix(ispxvx);
                    IndexVx const ivx(ispxvx);
                    return Kokkos::abs(
                            allfdistribu(ispxvx)
                            - Kokkos::cos(
                                    ddc::coordinate(ix) - advection_field(ispxvx) * timestep));
                });
        return max_advection_error;
    }
};

} // end namespace


TEST_F(Spatial1DAdvectionTest, SpatialAdvection)
{
    IDomainSpXVx meshSpXVx(dom_allsp, x_dom, vx_dom);

    SplineXBuilder const builder_x(meshSpXVx);

    ddc::PeriodicExtrapolationRule<RDimX> bv_x_min;
    ddc::PeriodicExtrapolationRule<RDimX> bv_x_max;
    SplineXEvaluator const spline_x_evaluator(bv_x_min, bv_x_max);

    PreallocatableSplineInterpolator const spline_x_interpolator(builder_x, spline_x_evaluator);

    Euler<FieldSpXVx<CoordX>, DFieldSpXVx> euler(meshSpXVx);
    BslAdvection1D<
            IDimX,
            IDomainSpXVx,
            IDomainSpXVx,
            SplineXBuilder,
            SplineXEvaluator,
            Euler<FieldSpXVx<CoordX>, DFieldSpXVx>> const
            spline_advection_x(spline_x_interpolator, builder_x, spline_x_evaluator, euler);

    double const err = SpatialAdvection(spline_advection_x);
    EXPECT_LE(err, 1.e-6);
    std::cout << "Max absolute difference to the exact function: " << err << std::endl;
}