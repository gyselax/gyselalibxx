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

using IdxRangeXVx = IdxRange<GridX, GridVx>;

using IdxRangeSp = IdxRange<Species>;
using IdxSp = Idx<Species>;
using IdxStepSp = IdxStep<Species>;

using IdxRangeSpXVx = IdxRange<Species, GridX, GridVx>;
using IdxSpXVx = Idx<Species, GridX, GridVx>;


template <class ElementType>
using FieldMemSpXVx = FieldMem<ElementType, IdxRangeSpXVx>;
using DFieldMemSpXVx = FieldMemSpXVx<double>;

template <class ElementType>
using FieldSp = Field<ElementType, IdxRangeSp>;
using DFieldSp = FieldSp<double>;

template <class ElementType>
using FieldSpXVx = Field<ElementType, IdxRangeSpXVx>;
using DFieldSpXVx = FieldSpXVx<double>;


template <class ElementType>
using FieldMemXVx = FieldMem<ElementType, IdxRangeXVx>;

template <class ElementType>
using FieldXVx = Field<ElementType, IdxRangeXVx>;



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


class Spatial1DAdvectionTest : public ::testing::Test
{
protected:
    IdxRangeX const idx_range_x;
    IdxRangeVx const idx_range_vx;
    IdxRangeSp const idx_range_allsp;

    static constexpr IdxStepSp nb_species = IdxStepSp(2);

public:
    Spatial1DAdvectionTest()
        : idx_range_x(SplineInterpPointsX::get_domain<GridX>())
        , idx_range_vx(SplineInterpPointsVx::get_domain<GridVx>())
        , idx_range_allsp(IdxSp(0), nb_species) {};

    ~Spatial1DAdvectionTest() = default;


    static void SetUpTestSuite()
    {
        CoordX const x_min(-M_PI);
        CoordX const x_max(M_PI);
        IdxStepX const x_size(100);

        CoordVx const vx_min(-6);
        CoordVx const vx_max(6);
        IdxStepVx const vx_size(50);

        ddc::init_discrete_space<BSplinesX>(x_min, x_max, x_size);
        ddc::init_discrete_space<BSplinesVx>(vx_min, vx_max, vx_size);

        ddc::init_discrete_space<GridX>(SplineInterpPointsX::get_sampling<GridX>());
        ddc::init_discrete_space<GridVx>(SplineInterpPointsVx::get_sampling<GridVx>());
    }

    template <class AdvectionOperator>
    double SpatialAdvection(AdvectionOperator const& advection_x)
    {
        // Mesh ----------------------------------------------------------------------------------
        IdxRangeSpXVx const meshSpXVx(idx_range_allsp, idx_range_x, idx_range_vx);
        IdxSp const i_elec = idx_range_allsp.front();
        IdxSp const i_ion = idx_range_allsp.back();


        // INITIALISATION ------------------------------------------------------------------------
        // Initialisation of the masses
        host_t<DFieldMemSp> masses_host(idx_range_allsp);
        masses_host(i_elec) = 1;
        masses_host(i_ion) = 1;
        auto masses_alloc = ddc::
                create_mirror_and_copy(Kokkos::DefaultExecutionSpace(), get_field(masses_host));
        DFieldSp masses = get_field(masses_alloc);


        // Initialisation of the distribution function and advection field
        DFieldMemSpXVx allfdistribu_alloc(meshSpXVx);
        DFieldSpXVx allfdistribu = get_field(allfdistribu_alloc);

        DFieldMemSpXVx advection_field_alloc(meshSpXVx);
        DFieldSpXVx advection_field = get_field(advection_field_alloc);

        ddc::parallel_for_each(
                Kokkos::DefaultExecutionSpace(),
                meshSpXVx,
                KOKKOS_LAMBDA(IdxSpXVx const ispxvx) {
                    IdxX const ix(ispxvx);
                    allfdistribu(ispxvx) = Kokkos::cos(ddc::coordinate(ix));

                    IdxVx const ivx(ispxvx);
                    IdxSp const isp(ispxvx);
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
                KOKKOS_LAMBDA(IdxSpXVx const ispxvx) {
                    IdxX const ix(ispxvx);
                    IdxVx const ivx(ispxvx);
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
    IdxRangeSpXVx meshSpXVx(idx_range_allsp, idx_range_x, idx_range_vx);

    SplineXBuilder const builder_x(meshSpXVx);

    ddc::PeriodicExtrapolationRule<X> bv_x_min;
    ddc::PeriodicExtrapolationRule<X> bv_x_max;
    SplineXEvaluator const spline_x_evaluator(bv_x_min, bv_x_max);

    PreallocatableSplineInterpolator const
            spline_x_interpolator(builder_x, spline_x_evaluator, meshSpXVx);

    Euler<FieldMemSpXVx<CoordX>, DFieldMemSpXVx> euler(meshSpXVx);
    BslAdvection1D<
            GridX,
            IdxRangeSpXVx,
            IdxRangeSpXVx,
            SplineXBuilder,
            SplineXEvaluator,
            Euler<FieldMemSpXVx<CoordX>, DFieldMemSpXVx>> const
            spline_advection_x(spline_x_interpolator, builder_x, spline_x_evaluator, euler);

    double const err = SpatialAdvection(spline_advection_x);
    EXPECT_LE(err, 1.e-6);
    std::cout << "Max absolute difference to the exact function: " << err << std::endl;
}
