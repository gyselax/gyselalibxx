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

// Splines
struct BSplinesX : ddc::UniformBSplines<X, 3>
{
};
struct BSplinesVx : ddc::UniformBSplines<Vx, 3>
{
};

ddc::BoundCond constexpr SplineXBoundary = ddc::BoundCond::PERIODIC;
ddc::BoundCond constexpr SplineVxBoundary = ddc::BoundCond::HERMITE;


// Discrete dimensions
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

using IdxRangeSp = IdxRange<Species>;
using IdxSp = Idx<Species>;
using IdxStepSp = IdxStep<Species>;

using IdxRangeSpX = IdxRange<Species, GridX>;
using IdxSpX = Idx<Species, GridX>;
using IdxStepSpX = IdxStep<Species, GridX>;

using IdxRangeSpXVx = IdxRange<Species, GridX, GridVx>;
using IdxSpXVx = Idx<Species, GridX, GridVx>;

// Chunks, Spans and Views
template <class ElementType>
using FieldMemSpXVx = FieldMem<ElementType, IdxRangeSpXVx>;
using DFieldMemSpXVx = FieldMemSpXVx<double>;

template <class ElementType>
using FieldSp = Field<ElementType, IdxRangeSp>;
using DFieldSp = FieldSp<double>;

template <class ElementType>
using FieldSpXVx = Field<ElementType, IdxRangeSpXVx>;
using DFieldSpXVx = FieldSpXVx<double>;


// For the derivatives
using IdxRangeSpVDerivVx = IdxRange<Species, GridX, ddc::Deriv<Vx>>;

template <class ElementType>
using DerivFieldMemSpX = FieldMem<ElementType, IdxRangeSpVDerivVx>;
using DDerivFieldMemSpX = DerivFieldMemSpX<double>;

template <class ElementType>
using DerivFieldSpX = Field<ElementType, IdxRangeSpVDerivVx>;
using DDerivFieldSpX = DerivFieldSpX<double>;



// Operators
using SplineVxBuilder = ddc::SplineBuilder<
        Kokkos::DefaultExecutionSpace,
        Kokkos::DefaultExecutionSpace::memory_space,
        BSplinesVx,
        GridVx,
        SplineVxBoundary,
        SplineVxBoundary,
        ddc::SplineSolver::LAPACK,
        Species,
        GridX,
        GridVx>;
using SplineVxEvaluator = ddc::SplineEvaluator<
        Kokkos::DefaultExecutionSpace,
        Kokkos::DefaultExecutionSpace::memory_space,
        BSplinesVx,
        GridVx,
        ddc::ConstantExtrapolationRule<Vx>,
        ddc::ConstantExtrapolationRule<Vx>,
        Species,
        GridX,
        GridVx>;


class Velocity1DAdvectionTest : public ::testing::Test
{
protected:
    IdxRangeX const x_dom;
    IdxRangeVx const vx_dom;
    IdxRangeSp const dom_allsp;

    static constexpr IdxStepSp nb_species = IdxStepSp(2);

public:
    Velocity1DAdvectionTest()
        : x_dom(SplineInterpPointsX::get_domain<GridX>())
        , vx_dom(SplineInterpPointsVx::get_domain<GridVx>())
        , dom_allsp(IdxSp(0), nb_species) {};

    ~Velocity1DAdvectionTest() = default;


    static void SetUpTestSuite()
    {
        CoordX const x_min(0);
        CoordX const x_max(2 * M_PI);
        IdxStepX const x_size(50);

        CoordVx const vx_min(-10);
        CoordVx const vx_max(10);
        IdxStepVx const vx_size(100);

        ddc::init_discrete_space<BSplinesX>(x_min, x_max, x_size);
        ddc::init_discrete_space<BSplinesVx>(vx_min, vx_max, vx_size);

        ddc::init_discrete_space<GridX>(SplineInterpPointsX::get_sampling<GridX>());
        ddc::init_discrete_space<GridVx>(SplineInterpPointsVx::get_sampling<GridVx>());
    }

    template <class AdvectionOperator>
    double VelocityAdvection(
            AdvectionOperator const& advection_x,
            SplineVxBuilder const& builder_vx)
    {
        // Mesh ----------------------------------------------------------------------------------
        IdxRangeSpXVx const meshSpXVx(dom_allsp, x_dom, vx_dom);
        IdxRangeSpX const meshSpX(dom_allsp, x_dom);
        IdxSp const i_elec = dom_allsp.front();
        IdxSp const i_ion = dom_allsp.back();


        // INITIALISATION ------------------------------------------------------------------------
        // Initialization of the charges
        host_t<DFieldMemSp> charges_host(dom_allsp);
        charges_host(i_elec) = -1.;
        charges_host(i_ion) = 1.;
        auto charges_alloc = ddc::
                create_mirror_and_copy(Kokkos::DefaultExecutionSpace(), get_field(charges_host));
        ddc::ChunkSpan charges = get_field(charges_alloc);

        // Initialization of the masses
        host_t<DFieldMemSp> masses_host(dom_allsp);
        masses_host(i_elec) = 1.;
        masses_host(i_ion) = 1.;
        auto masses_alloc = ddc::
                create_mirror_and_copy(Kokkos::DefaultExecutionSpace(), get_field(masses_host));
        ddc::ChunkSpan masses = get_field(masses_alloc);


        // Initialization of the distribution function and advection field
        DFieldMemSpXVx allfdistribu_alloc(meshSpXVx);
        DFieldSpXVx allfdistribu = get_field(allfdistribu_alloc);

        DFieldMemSpXVx advection_field_alloc(meshSpXVx);
        DFieldSpXVx advection_field = get_field(advection_field_alloc);

        ddc::parallel_for_each(
                Kokkos::DefaultExecutionSpace(),
                meshSpXVx,
                KOKKOS_LAMBDA(IdxSpXVx const ispxvx) {
                    double const v = ddc::coordinate(IdxVx(ispxvx));
                    allfdistribu(ispxvx) = Kokkos::exp(-0.5 * v * v);

                    IdxX const ix(ispxvx);
                    IdxSp const isp(ispxvx);
                    double const electric_field = ddc::distance_at_right(ix);
                    advection_field(ispxvx) = charges(isp)
                                              * Kokkos::sqrt(masses(i_elec) / masses(isp))
                                              * electric_field;
                });


        // Initialization of the derivatives of the advection field
        DDerivFieldMemSpX advection_field_derivatives_min_alloc(
                builder_vx.batched_derivs_xmin_domain());
        DDerivFieldSpX advection_field_derivatives_min
                = get_field(advection_field_derivatives_min_alloc);

        DDerivFieldMemSpX advection_field_derivatives_max_alloc(
                builder_vx.batched_derivs_xmax_domain());
        DDerivFieldSpX advection_field_derivatives_max
                = get_field(advection_field_derivatives_max_alloc);

        ddc::parallel_fill(advection_field_derivatives_min_alloc, 1.);
        ddc::parallel_fill(advection_field_derivatives_max_alloc, 1.);


        // SIMULATION ----------------------------------------------------------------------------
        double const timestep = .1;

        advection_x(
                allfdistribu,
                advection_field,
                timestep,
                get_const_field(advection_field_derivatives_min),
                get_const_field(advection_field_derivatives_max));


        double const max_advection_error = ddc::parallel_transform_reduce(
                Kokkos::DefaultExecutionSpace(),
                meshSpXVx,
                0.0,
                ddc::reducer::max<double>(),
                KOKKOS_LAMBDA(IdxSpXVx const ispxvx) {
                    IdxX const ix(ispxvx);
                    IdxVx const ivx(ispxvx);
                    double const v = ddc::coordinate(IdxVx(ispxvx));
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
    IdxRangeSpXVx meshSpXVx(dom_allsp, x_dom, vx_dom);

    // Interpolator for the function
    IdxStepVx static constexpr n_ghosts_vx {0};

    LagrangeInterpolator<GridVx, BCond::DIRICHLET, BCond::DIRICHLET, Species, GridX, GridVx> const
            lagrange_vx_non_preallocatable_interpolator(3, n_ghosts_vx);
    PreallocatableLagrangeInterpolator<
            GridVx,
            BCond::DIRICHLET,
            BCond::DIRICHLET,
            Species,
            GridX,
            GridVx> const lagrange_vx_interpolator(lagrange_vx_non_preallocatable_interpolator);

    // Interpolator of the advection field
    SplineVxBuilder const builder_vx(meshSpXVx);

    CoordVx const vx_min = ddc::coordinate(vx_dom.front());
    CoordVx const vx_max = vx_min + ddcHelper::total_interval_length(vx_dom);
    ddc::ConstantExtrapolationRule<Vx> bv_v_min(vx_min);
    ddc::ConstantExtrapolationRule<Vx> bv_v_max(vx_max);
    SplineVxEvaluator const spline_vx_evaluator(bv_v_min, bv_v_max);


    Euler<FieldMemSpXVx<CoordVx>, DFieldMemSpXVx> euler(meshSpXVx);
    BslAdvection1D<
            GridVx,
            IdxRangeSpXVx,
            IdxRangeSpXVx,
            SplineVxBuilder,
            SplineVxEvaluator,
            Euler<FieldMemSpXVx<CoordVx>, DFieldMemSpXVx>> const
            lagrange_advection_vx(lagrange_vx_interpolator, builder_vx, spline_vx_evaluator, euler);


    double const err = VelocityAdvection(lagrange_advection_vx, builder_vx);
    EXPECT_LE(err, 1e-3);
    std::cout << "Max absolute difference to the exact function: " << err << std::endl;
}


TEST_F(Velocity1DAdvectionTest, SplineBatched)
{
    IdxRangeSpXVx meshSpXVx(dom_allsp, x_dom, vx_dom);

    SplineVxBuilder const builder_vx(meshSpXVx);

    CoordVx const vx_min = ddc::coordinate(vx_dom.front());
    CoordVx const vx_max = vx_min + ddcHelper::total_interval_length(vx_dom);

    ddc::ConstantExtrapolationRule<Vx> bv_v_min(vx_min);
    ddc::ConstantExtrapolationRule<Vx> bv_v_max(vx_max);
    SplineVxEvaluator const spline_vx_evaluator(bv_v_min, bv_v_max);

    PreallocatableSplineInterpolator const spline_vx_interpolator(builder_vx, spline_vx_evaluator);

    Euler<FieldMemSpXVx<CoordVx>, DFieldMemSpXVx> euler(meshSpXVx);
    BslAdvection1D<
            GridVx,
            IdxRangeSpXVx,
            IdxRangeSpXVx,
            SplineVxBuilder,
            SplineVxEvaluator,
            Euler<FieldMemSpXVx<CoordVx>, DFieldMemSpXVx>> const
            spline_advection_vx(spline_vx_interpolator, builder_vx, spline_vx_evaluator, euler);


    double const err = VelocityAdvection(spline_advection_vx, builder_vx);
    EXPECT_LE(err, 1.e-5);
    std::cout << "Max absolute difference to the exact function: " << err << std::endl;
}
