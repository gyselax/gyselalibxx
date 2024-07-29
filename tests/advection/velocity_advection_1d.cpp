/*
    Advection along Vx on (Sp, X, Vx). The advection field is given by function on X.
    The test is the same as the one applied to BslAdvectionVelocity operator.
*/

#include <ddc/ddc.hpp>

#include <gtest/gtest.h>

#include "Lagrange_interpolator.hpp"
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

using IdxRangeSp = ddc::DiscreteDomain<Species>;
using IdxSp = ddc::DiscreteElement<Species>;
using IdxStepSp = ddc::DiscreteVector<Species>;

using IDomainSpX = ddc::DiscreteDomain<Species, IDimX>;
using IndexSpX = ddc::DiscreteElement<Species, IDimX>;
using IVectSpX = ddc::DiscreteVector<Species, IDimX>;

using IDomainSpXVx = ddc::DiscreteDomain<Species, IDimX, IDimVx>;
using IndexSpXVx = ddc::DiscreteElement<Species, IDimX, IDimVx>;

// Chunks, Spans and Views
template <class ElementType>
using FieldMemSp = device_t<ddc::Chunk<ElementType, IdxRangeSp>>;
using DFieldMemSp = DFieldMemSp;

template <class ElementType>
using FieldSpXVx = device_t<ddc::Chunk<ElementType, IDomainSpXVx>>;
using DFieldSpXVx = FieldSpXVx<double>;

template <class ElementType>
using FieldSp = device_t<ddc::ChunkSpan<ElementType, IdxRangeSp>>;
using DFieldSp = FieldSp<double>;

template <class ElementType>
using SpanSpXVx = device_t<ddc::ChunkSpan<ElementType, IDomainSpXVx>>;
using DSpanSpXVx = SpanSpXVx<double>;


// For the derivatives
using IDomainSpVDerivVx = ddc::DiscreteDomain<Species, IDimX, ddc::Deriv<RDimVx>>;

template <class ElementType>
using DerivFieldSpX = device_t<ddc::Chunk<ElementType, IDomainSpVDerivVx>>;
using DDerivFieldSpX = DerivFieldSpX<double>;

template <class ElementType>
using DerivSpanSpX = device_t<ddc::ChunkSpan<ElementType, IDomainSpVDerivVx>>;
using DDerivSpanSpX = DerivSpanSpX<double>;



// Operators
using SplineVxBuilder = ddc::SplineBuilder<
        Kokkos::DefaultExecutionSpace,
        Kokkos::DefaultExecutionSpace::memory_space,
        BSplinesVx,
        IDimVx,
        SplineVxBoundary,
        SplineVxBoundary,
        ddc::SplineSolver::LAPACK,
        Species,
        IDimX,
        IDimVx>;
using SplineVxEvaluator = ddc::SplineEvaluator<
        Kokkos::DefaultExecutionSpace,
        Kokkos::DefaultExecutionSpace::memory_space,
        BSplinesVx,
        IDimVx,
        ddc::ConstantExtrapolationRule<RDimVx>,
        ddc::ConstantExtrapolationRule<RDimVx>,
        Species,
        IDimX,
        IDimVx>;


class Velocity1DAdvectionTest : public ::testing::Test
{
protected:
    IDomainX const x_dom;
    IDomainVx const vx_dom;
    IdxRangeSp const dom_allsp;

    static constexpr IdxStepSp nb_species = IdxStepSp(2);

public:
    Velocity1DAdvectionTest()
        : x_dom(SplineInterpPointsX::get_domain<IDimX>())
        , vx_dom(SplineInterpPointsVx::get_domain<IDimVx>())
        , dom_allsp(IdxSp(0), nb_species) {};

    ~Velocity1DAdvectionTest() = default;


    static void SetUpTestSuite()
    {
        CoordX const x_min(0);
        CoordX const x_max(2 * M_PI);
        IVectX const x_size(50);

        CoordVx const vx_min(-10);
        CoordVx const vx_max(10);
        IVectVx const vx_size(100);

        ddc::init_discrete_space<BSplinesX>(x_min, x_max, x_size);
        ddc::init_discrete_space<BSplinesVx>(vx_min, vx_max, vx_size);

        ddc::init_discrete_space<IDimX>(SplineInterpPointsX::get_sampling<IDimX>());
        ddc::init_discrete_space<IDimVx>(SplineInterpPointsVx::get_sampling<IDimVx>());
    }

    template <class AdvectionOperator>
    double VelocityAdvection(
            AdvectionOperator const& advection_x,
            SplineVxBuilder const& builder_vx)
    {
        // Mesh ----------------------------------------------------------------------------------
        IDomainSpXVx const meshSpXVx(dom_allsp, x_dom, vx_dom);
        IDomainSpX const meshSpX(dom_allsp, x_dom);
        IdxSp const i_elec = dom_allsp.front();
        IdxSp const i_ion = dom_allsp.back();


        // INITIALISATION ------------------------------------------------------------------------
        // Initialization of the charges
        host_t<DFieldMemSp> charges_host(dom_allsp);
        charges_host(i_elec) = -1.;
        charges_host(i_ion) = 1.;
        auto charges_alloc = ddc::
                create_mirror_and_copy(Kokkos::DefaultExecutionSpace(), charges_host.span_view());
        ddc::ChunkSpan charges = charges_alloc.span_view();

        // Initialization of the masses
        host_t<DFieldMemSp> masses_host(dom_allsp);
        masses_host(i_elec) = 1.;
        masses_host(i_ion) = 1.;
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
                    double const v = ddc::coordinate(IndexVx(ispxvx));
                    allfdistribu(ispxvx) = Kokkos::exp(-0.5 * v * v);

                    IndexX const ix(ispxvx);
                    IdxSp const isp(ispxvx);
                    double const electric_field = ddc::distance_at_right(ix);
                    advection_field(ispxvx) = charges(isp)
                                              * Kokkos::sqrt(masses(i_elec) / masses(isp))
                                              * electric_field;
                });


        // Initialization of the derivatives of the advection field
        DDerivFieldSpX advection_field_derivatives_min_alloc(
                builder_vx.batched_derivs_xmin_domain());
        DDerivSpanSpX advection_field_derivatives_min
                = advection_field_derivatives_min_alloc.span_view();

        DDerivFieldSpX advection_field_derivatives_max_alloc(
                builder_vx.batched_derivs_xmax_domain());
        DDerivSpanSpX advection_field_derivatives_max
                = advection_field_derivatives_max_alloc.span_view();

        ddc::parallel_fill(advection_field_derivatives_min_alloc, 1.);
        ddc::parallel_fill(advection_field_derivatives_max_alloc, 1.);


        // SIMULATION ----------------------------------------------------------------------------
        double const timestep = .1;

        advection_x(
                allfdistribu,
                advection_field,
                timestep,
                advection_field_derivatives_min.span_cview(),
                advection_field_derivatives_max.span_cview());


        double const max_advection_error = ddc::parallel_transform_reduce(
                Kokkos::DefaultExecutionSpace(),
                meshSpXVx,
                0.0,
                ddc::reducer::max<double>(),
                KOKKOS_LAMBDA(IndexSpXVx const ispxvx) {
                    IndexX const ix(ispxvx);
                    IndexVx const ivx(ispxvx);
                    double const v = ddc::coordinate(IndexVx(ispxvx));
                    return Kokkos::abs(
                            allfdistribu(ispxvx)
                            - Kokkos::exp(
                                    -0.5 * Kokkos::pow(v - advection_field(ispxvx) * timestep, 2)));
                });
        return max_advection_error;

        return 0;
    }
};

} // namespace


TEST_F(Velocity1DAdvectionTest, BatchedLagrange)
{
    IDomainSpXVx meshSpXVx(dom_allsp, x_dom, vx_dom);

    // Interpolator for the function
    IVectVx static constexpr n_ghosts_vx {0};

    LagrangeInterpolator<IDimVx, BCond::DIRICHLET, BCond::DIRICHLET, Species, IDimX, IDimVx> const
            lagrange_vx_non_preallocatable_interpolator(3, n_ghosts_vx);
    PreallocatableLagrangeInterpolator<
            IDimVx,
            BCond::DIRICHLET,
            BCond::DIRICHLET,
            Species,
            IDimX,
            IDimVx> const lagrange_vx_interpolator(lagrange_vx_non_preallocatable_interpolator);

    // Interpolator of the advection field
    SplineVxBuilder const builder_vx(meshSpXVx);

    CoordVx const vx_min = ddc::coordinate(vx_dom.front());
    CoordVx const vx_max = vx_min + ddcHelper::total_interval_length(vx_dom);
    ddc::ConstantExtrapolationRule<RDimVx> bv_v_min(vx_min);
    ddc::ConstantExtrapolationRule<RDimVx> bv_v_max(vx_max);
    SplineVxEvaluator const spline_vx_evaluator(bv_v_min, bv_v_max);


    Euler<FieldSpXVx<CoordVx>, DFieldSpXVx> euler(meshSpXVx);
    BslAdvection1D<
            IDimVx,
            IDomainSpXVx,
            IDomainSpXVx,
            SplineVxBuilder,
            SplineVxEvaluator,
            Euler<FieldSpXVx<CoordVx>, DFieldSpXVx>> const
            lagrange_advection_vx(lagrange_vx_interpolator, builder_vx, spline_vx_evaluator, euler);


    double const err = VelocityAdvection(lagrange_advection_vx, builder_vx);
    EXPECT_LE(err, 1e-3);
    std::cout << "Max absolute difference to the exact function: " << err << std::endl;
}


TEST_F(Velocity1DAdvectionTest, SplineBatched)
{
    IDomainSpXVx meshSpXVx(dom_allsp, x_dom, vx_dom);

    SplineVxBuilder const builder_vx(meshSpXVx);

    CoordVx const vx_min = ddc::coordinate(vx_dom.front());
    CoordVx const vx_max = vx_min + ddcHelper::total_interval_length(vx_dom);

    ddc::ConstantExtrapolationRule<RDimVx> bv_v_min(vx_min);
    ddc::ConstantExtrapolationRule<RDimVx> bv_v_max(vx_max);
    SplineVxEvaluator const spline_vx_evaluator(bv_v_min, bv_v_max);

    PreallocatableSplineInterpolator const spline_vx_interpolator(builder_vx, spline_vx_evaluator);

    Euler<FieldSpXVx<CoordVx>, DFieldSpXVx> euler(meshSpXVx);
    BslAdvection1D<
            IDimVx,
            IDomainSpXVx,
            IDomainSpXVx,
            SplineVxBuilder,
            SplineVxEvaluator,
            Euler<FieldSpXVx<CoordVx>, DFieldSpXVx>> const
            spline_advection_vx(spline_vx_interpolator, builder_vx, spline_vx_evaluator, euler);


    double const err = VelocityAdvection(spline_advection_vx, builder_vx);
    EXPECT_LE(err, 1.e-5);
    std::cout << "Max absolute difference to the exact function: " << err << std::endl;
}
