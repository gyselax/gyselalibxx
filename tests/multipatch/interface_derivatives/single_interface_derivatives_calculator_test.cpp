// SPDX-License-Identifier: MIT
#include <ddc/ddc.hpp>
#include <ddc/kernels/splines.hpp>

#include <gtest/gtest.h>

#if defined(UNIFORM_MESH)
#include "2patches_2d_onion_shape_uniform.hpp"
#else
#include "2patches_2d_onion_shape_non_uniform.hpp"
#endif
#include "../../test_utils.hpp"

#include "ddc_alias_inline_functions.hpp"
#include "ddc_helper.hpp"
#include "global_2d_onion_shape_non_uniform.hpp"
#include "interface.hpp"
#include "mesh_builder.hpp"
#include "non_uniform_interpolation_points.hpp"
#include "single_interface_derivatives_calculator.hpp"
#include "interface_derivatives_test_utils.hpp"



/*
    Test: 
    * Interpolation with Hermite boundary conditions. 
        * Non-uniform per patch case => test the recursive formula. 
        * Uniform per patch case => test the explicit formula. 
        * Exact formula for 30, 20, 25, 10 and 5 cells. 
        * Approximation formula for 30, 20, 25, 10 and 5 cells. 

    * Interpolation points as closure condition. 
        * Non-uniform per patch case => test the recursive formula. 
        * Exact formula 5|10 cells (we only need to test the boundaries). 
    
    * Different connections between the patches: 
        * Interface East 1 | West 2.  
        * Interface East 1 | East 2.  
        * Interface West 1 | East 2.  
        * Interface West 1 | West 2.  
        * Interface East 1 | South 2.  
*/

namespace {
#if defined(UNIFORM_MESH)
using namespace onion_shape_uniform_2d_2patches;
#else
using namespace onion_shape_non_uniform_2d_2patches;
#endif
using namespace onion_shape_non_uniform_2d_global;


using HostExecSpace = Kokkos::DefaultHostExecutionSpace;


// Rename edge structures for a better readability in the TEST names.
struct EastEdge1 : EastEdge<1>
{
};
struct WestEdge1 : WestEdge<1>
{
};
struct EastEdge2 : EastEdge<2>
{
};
struct WestEdge2 : WestEdge<2>
{
};
struct SouthEdge2 : SouthEdge<2>
{
};

// Rename continuous dimensions of Patch2 to avoid confusion with East1|South2 connection.
using Eta = R;
using Xi = Theta;


// Interpolation points type for the 2 patches
// --- patch 1
template <ddc::BoundCond BoundCondR>
using SplineInterpPointsR1 = ddcHelper::
        NonUniformInterpolationPoints<BSplinesR<1>, BoundCondR, ddc::BoundCond::HERMITE>;
using SplineInterpPointsTheta1 = ddc::KnotsAsInterpolationPoints<
        BSplinesTheta<1>,
        ddc::BoundCond::PERIODIC,
        ddc::BoundCond::PERIODIC>;

// --- patch 2
template <ddc::BoundCond BoundCondR>
using SplineInterpPointsEta2 = ddcHelper::
        NonUniformInterpolationPoints<BSplinesR<2>, ddc::BoundCond::HERMITE, BoundCondR>;
using SplineInterpPointsXi2 = ddc::KnotsAsInterpolationPoints<
        BSplinesTheta<2>,
        ddc::BoundCond::PERIODIC,
        ddc::BoundCond::PERIODIC>;


// Interpolation points type for the equivalent global spline.
template <ddc::BoundCond BoundCondR>
using SplineInterpPointsRg
        = ddcHelper::NonUniformInterpolationPoints<BSplinesRg, BoundCondR, BoundCondR>;
using SplineInterpPointsThetag = ddcHelper::NonUniformInterpolationPoints<
        BSplinesThetag,
        ddc::BoundCond::PERIODIC,
        ddc::BoundCond::PERIODIC>;


// Operators on the equivalent global spline
template <ddc::BoundCond BoundCondR>
using SplineRThetagBuilder = ddc::SplineBuilder2D<
        HostExecSpace,
        typename HostExecSpace::memory_space,
        BSplinesRg,
        BSplinesThetag,
        GridRg,
        GridThetag,
        BoundCondR,
        BoundCondR,
        ddc::BoundCond::PERIODIC,
        ddc::BoundCond::PERIODIC,
        ddc::SplineSolver::LAPACK>;

using SplineRThetagEvaluator = ddc::SplineEvaluator2D<
        HostExecSpace,
        typename HostExecSpace::memory_space,
        BSplinesRg,
        BSplinesThetag,
        GridRg,
        GridThetag,
        ddc::ConstantExtrapolationRule<Rg, Thetag>,
        ddc::ConstantExtrapolationRule<Rg, Thetag>,
        ddc::PeriodicExtrapolationRule<Thetag>,
        ddc::PeriodicExtrapolationRule<Thetag>>;


<<<<<<< HEAD
=======
/**
 *  @brief Get interpolation points from the break points by placing 
 * the interpolation points on the break points and adding one on the
 * left boundary cell at 2/3 of the cell. 
 */
template <class CoordType>
std::vector<CoordType> get_interpolation_points_add_one_on_left(
        std::vector<CoordType> const& break_points)
{
    CoordType additional_point(break_points[0] * 2. / 3. + break_points[1] * 1. / 3.);
    std::vector<CoordType> interpolation_points(break_points);
    interpolation_points.insert(interpolation_points.begin() + 1, additional_point);
    return interpolation_points;
}

/**
 *  @brief Get interpolation points from the break points by placing 
 * the interpolation points on the break points and adding one on the
 * right boundary cell at 1/3 of the cell. 
 */
template <class CoordType>
std::vector<CoordType> get_interpolation_points_add_one_on_right(
        std::vector<CoordType> const& break_points)
{
    int n_bpoints = break_points.size();
    CoordType additional_point(
            break_points[n_bpoints - 1] * 2. / 3. + break_points[n_bpoints - 2] * 1. / 3.);
    std::vector<CoordType> interpolation_points(break_points);
    interpolation_points.insert(interpolation_points.end() - 1, additional_point);
    return interpolation_points;
}

/**
 * @brief Fill in a vector of points for the equivalent global mesh 
 * by conserving the same order of the given points.
 */
template <class CoordTypeG, class CoordTypeP>
void fill_in(std::vector<CoordTypeG>& points_global, std::vector<CoordTypeP> const& points_patch)
{
    for (CoordTypeP pt : points_patch) {
        points_global.push_back(CoordTypeG {double(pt)});
    }
}

/**
 * @brief Fill in a vector of points for the equivalent global mesh
 *  by reversing the order of the given points.
 */
template <class CoordTypeG, class CoordTypeP>
void fill_in_reverse(
        std::vector<CoordTypeG>& points_global,
        std::vector<CoordTypeP> const& points_patch)
{
    std::size_t const n_pt = points_patch.size();
    CoordTypeP const max = points_patch[n_pt - 1];
    CoordTypeP const min = points_patch[0];
    for (int i(0); i < n_pt; ++i) {
        points_global.push_back(Coord<Rg> {double(min + max - points_patch[n_pt - 1 - i])});
    }
}

>>>>>>> main_pvidal
/// @brief Initialise the function with f(r,theta) = r(3-r)sin(theta).
template <class Grid1, class Grid2>
void initialise_2D_function(host_t<DField<IdxRange<Grid1, Grid2>>> function)
{
    ddc::for_each(get_idx_range(function), [&](Idx<Grid1, Grid2> idx) {
        // Get the coordinate on the equivalent global mesh.
        double const rg = ddc::coordinate(Idx<Grid1>(idx));
        double const thetag = ddc::coordinate(Idx<Grid2>(idx));
        function(idx) = rg * (3. - rg) * Kokkos::sin(thetag);
    });
}



template <class T>
struct SingleInterfaceDerivativesCalculatorFixture;

template <class InterpolationType, class Edge_Patch1, class Edge_Patch2>
struct SingleInterfaceDerivativesCalculatorFixture<
        std::tuple<InterpolationType, Edge_Patch1, Edge_Patch2>> : public ::testing::Test
{
    // Get the parameters of the test: patch connection and interpolation type.
    using Interpolation = InterpolationType;
    static constexpr ddc::BoundCond Interpolation_v = InterpolationType::value;
    using Edge1 = Edge_Patch1;
    using Edge2 = Edge_Patch2;


    // DEFINE BOUNDARIES OF THE DOMAINS ----------------------------------------------------------
    // patch 1 -----------------------------------
    static constexpr Coord<R> r1_min = Coord<R>(0.0);
    static constexpr Coord<R> r1_max = Coord<R>(1.0);
    // Select less cells for ddc::BoundCond::GREVILLE to test the boundary.
    static constexpr IdxStep<GridR<1>> r1_ncells = (Interpolation_v == ddc::BoundCond::GREVILLE)
                                                           ? IdxStep<GridR<1>>(5)
                                                           : IdxStep<GridR<1>>(30);

    static constexpr Coord<Theta> theta1_min = Coord<Theta>(0.0);
    static constexpr Coord<Theta> theta1_max = Coord<Theta>(2 * M_PI);
    static constexpr IdxStep<GridTheta<1>> theta1_ncells = IdxStep<GridTheta<1>>(30);

    // patch 2 -----------------------------------
    // Exchange R and Theta values for a connection with SouthEdge2.
    static constexpr Coord<Eta> eta2_min
            = (std::is_same_v<Edge2, SouthEdge2>) ? Coord<Eta>(1.0 * 2 * M_PI) : Coord<Eta>(1.0);
    static constexpr Coord<R> eta2_max
            = (std::is_same_v<Edge2, SouthEdge2>) ? Coord<Eta>(2.0 * 2 * M_PI) : Coord<Eta>(2.0);
    // Select less cells for ddc::BoundCond::GREVILLE to test the boundary.
    static constexpr IdxStep<GridR<2>> eta2_ncells = (Interpolation_v == ddc::BoundCond::GREVILLE)
                                                             ? IdxStep<GridR<2>>(5)
                                                             : IdxStep<GridR<2>>(30);

    // Exchange R and Theta values for a connection with SouthEdge2.
    static constexpr Coord<Theta> xi2_min = Coord<Xi>(0.0);
    static constexpr Coord<Theta> xi2_max
            = (std::is_same_v<Edge2, SouthEdge2>) ? Coord<Xi>(1.0) : Coord<Xi>(2 * M_PI);
    static constexpr IdxStep<GridTheta<2>> xi2_ncells = IdxStep<GridTheta<2>>(30);


    // global ------------------------------------
    static constexpr Coord<Rg> rg_min = Coord<Rg> {double(r1_min)};
    static constexpr Coord<Rg> rg_max = Coord<Rg> {double(eta2_max)};
    static constexpr IdxStep<GridRg> rg_ncells
            = (std::is_same_v<Edge2, SouthEdge2>)
                      ? IdxStep<GridRg>(r1_ncells.value() + xi2_ncells.value())
                      : IdxStep<GridRg>(r1_ncells.value() + eta2_ncells.value());

    static constexpr Coord<Thetag> thetag_min = Coord<Thetag>(0.0);
    static constexpr Coord<Thetag> thetag_max = Coord<Thetag>(2 * M_PI);
    static constexpr IdxStep<GridThetag> thetag_ncells = IdxStep<GridThetag>(30);


    // INITIALISE DOMAINS ------------------------------------------------------------------------
    static void SetUpTestSuite()
    {
        // Creating of meshes and supports .......................................................
        // Patch 1 ...............................................................................
        std::vector<Coord<R>> break_points_r1;
#if defined(UNIFORM_MESH)
        break_points_r1 = build_uniform_break_points(r1_min, r1_max, r1_ncells);
        ddc::init_discrete_space<BSplinesR<1>>(r1_min, r1_max, r1_ncells);
        ddc::init_discrete_space<BSplinesTheta<1>>(theta1_min, theta1_max, theta1_ncells);
#else
        break_points_r1 = build_random_non_uniform_break_points(r1_min, r1_max, r1_ncells);
        ddc::init_discrete_space<BSplinesR<1>>(break_points_r1);
        ddc::init_discrete_space<BSplinesTheta<1>>(
                build_uniform_break_points(theta1_min, theta1_max, theta1_ncells));
#endif

        std::vector<Coord<R>> interpolation_points_r1;
        if constexpr (Interpolation_v == ddc::BoundCond::GREVILLE) {
            // Add an interpolation point in the cell on the boundary of the equivalent global domain.
            // The other interpolation points are the break points.
            if constexpr (std::is_same_v<Edge1, EastEdge1>) {
                interpolation_points_r1 = get_interpolation_points_add_one_on_left(break_points_r1);
            } else if (std::is_same_v<Edge1, WestEdge1>) {
                interpolation_points_r1
                        = get_interpolation_points_add_one_on_right(break_points_r1);
            }
        } else if (Interpolation_v == ddc::BoundCond::HERMITE) {
            // Use the break points as interpolation points.
            interpolation_points_r1 = break_points_r1;
        }
        ddc::init_discrete_space<GridR<1>>(
                SplineInterpPointsR1<Interpolation_v>::template get_sampling<GridR<1>>(
                        interpolation_points_r1));
        ddc::init_discrete_space<GridTheta<1>>(
                SplineInterpPointsTheta1::get_sampling<GridTheta<1>>());


        // Patch 2 ...............................................................................
        std::vector<Coord<R>> break_points_eta2;
        std::vector<Coord<Theta>> break_points_xi2;
#if defined(UNIFORM_MESH)
        break_points_eta2 = build_uniform_break_points(eta2_min, eta2_max, eta2_ncells);
        break_points_xi2 = build_uniform_break_points(xi2_min, xi2_max, xi2_ncells);
        ddc::init_discrete_space<BSplinesR<2>>(eta2_min, eta2_max, eta2_ncells);
        ddc::init_discrete_space<BSplinesTheta<2>>(xi2_min, xi2_max, xi2_ncells);
#else
        if constexpr (std::is_same_v<Edge2, SouthEdge2>) {
            // The break points have to match with the ones of the equivalent global domain.
            break_points_eta2 = build_uniform_break_points(eta2_min, eta2_max, eta2_ncells);
            break_points_xi2 = build_random_non_uniform_break_points(xi2_min, xi2_max, xi2_ncells);
        } else {
            break_points_eta2
                    = build_random_non_uniform_break_points(eta2_min, eta2_max, eta2_ncells);
            break_points_xi2 = build_uniform_break_points(xi2_min, xi2_max, xi2_ncells);
        }
        ddc::init_discrete_space<BSplinesR<2>>(break_points_eta2);
        ddc::init_discrete_space<BSplinesTheta<2>>(break_points_xi2);
#endif

        std::vector<Coord<R>> interpolation_points_eta2;
        std::vector<Coord<Theta>> interpolation_points_xi2;
        if constexpr (Interpolation_v == ddc::BoundCond::GREVILLE) {
            // Add an interpolation point in the cell on the boundary of the equivalent global domain.
            // The other interpolation points are the break points.
            if constexpr (std::is_same_v<Edge2, EastEdge2>) {
                interpolation_points_eta2
                        = get_interpolation_points_add_one_on_left(break_points_eta2);
                interpolation_points_xi2 = break_points_xi2;
            } else if (std::is_same_v<Edge2, WestEdge2>) {
                interpolation_points_eta2
                        = get_interpolation_points_add_one_on_right(break_points_eta2);
                interpolation_points_xi2 = break_points_xi2;
            } else if (std::is_same_v<Edge2, SouthEdge2>) {
                interpolation_points_xi2
                        = get_interpolation_points_add_one_on_right(break_points_xi2);
                interpolation_points_eta2 = break_points_eta2;
            }
        } else if (Interpolation_v == ddc::BoundCond::HERMITE) {
            // Use the break points as interpolation points.
            interpolation_points_eta2 = break_points_eta2;
            interpolation_points_xi2 = break_points_xi2;
        }
        ddc::init_discrete_space<GridR<2>>(
                SplineInterpPointsEta2<Interpolation_v>::template get_sampling<GridR<2>>(
                        interpolation_points_eta2));
        ddc::init_discrete_space<GridTheta<2>>(SplineInterpPointsXi2::get_sampling<GridTheta<2>>());


        // Equivalent global domain ..............................................................
        std::vector<Coord<Rg>> break_points_rg;
        std::vector<Coord<Rg>> interpolation_points_rg;
        // --- fill in with points from patch 1
        if constexpr (std::is_same_v<Edge1, EastEdge1>) {
            // orientation global: ↑→  | local: ↑→
            break_points_r1.pop_back();
            interpolation_points_r1.pop_back();
            fill_in(break_points_rg, break_points_r1);
            fill_in(interpolation_points_rg, interpolation_points_r1);
        } else if (std::is_same_v<Edge1, WestEdge1>) {
            // orientation global: ↑→  | local: ←↓
            fill_in_reverse(break_points_rg, break_points_r1);
            fill_in_reverse(interpolation_points_rg, interpolation_points_r1);
            break_points_rg.pop_back();
            interpolation_points_rg.pop_back();
        }

        // --- fill in with points from patch 2
        if constexpr (std::is_same_v<Edge2, EastEdge2>) {
            // orientation global: ↑→  | local: ←↓
            fill_in_reverse(break_points_rg, break_points_eta2);
            fill_in_reverse(interpolation_points_rg, interpolation_points_eta2);
        } else if (std::is_same_v<Edge2, WestEdge2>) {
            // orientation global: ↑→  | local: ↑→
            fill_in(break_points_rg, break_points_eta2);
            fill_in(interpolation_points_rg, interpolation_points_eta2);
        } else if (std::is_same_v<Edge2, SouthEdge2>) {
            // orientation global: ↑→  | local: ↓→
            for (typename Patch2::Coord2 xi : break_points_xi2) {
                break_points_rg.push_back(Coord<Rg> {double(xi / xi2_max) + double(r1_max)});
            }
            for (typename Patch2::Coord2 xi : interpolation_points_xi2) {
                interpolation_points_rg.push_back(
                        Coord<Rg> {double(xi / xi2_max) + double(r1_max)});
            }
        }
        std::vector<Coord<Thetag>> break_points_thetag
                = build_uniform_break_points(thetag_min, thetag_max, thetag_ncells);
        std::vector<Coord<Thetag>> interpolation_points_thetag = break_points_thetag;
        interpolation_points_thetag.pop_back();

        ddc::init_discrete_space<BSplinesRg>(break_points_rg);
        ddc::init_discrete_space<BSplinesThetag>(break_points_thetag);

        ddc::init_discrete_space<GridRg>(
                SplineInterpPointsRg<Interpolation_v>::template get_sampling<GridRg>(
                        interpolation_points_rg));
        ddc::init_discrete_space<GridThetag>(
                SplineInterpPointsThetag::get_sampling<GridThetag>(interpolation_points_thetag));
    }


    // DEFINE INDEX RANGES -----------------------------------------------------------------------
    using SplineInterpPoints_R1 = SplineInterpPointsR1<Interpolation_v>;
    const typename Patch1::IdxRange1 idx_range_r1 =
            typename Patch1::IdxRange1(SplineInterpPoints_R1::template get_domain<GridR<1>>());
    const typename Patch1::IdxRange2 idx_range_theta1 =
            typename Patch1::IdxRange2(SplineInterpPointsTheta1::get_domain<GridTheta<1>>());
    const typename Patch1::IdxRange12 idx_range_rtheta1 =
            typename Patch1::IdxRange12(idx_range_r1, idx_range_theta1);

    using SplineInterpPoints_R2 = SplineInterpPointsEta2<Interpolation_v>;
    const typename Patch2::IdxRange1 idx_range_eta2 =
            typename Patch2::IdxRange1(SplineInterpPoints_R2::template get_domain<GridR<2>>());
    const typename Patch2::IdxRange2 idx_range_xi2 =
            typename Patch2::IdxRange2(SplineInterpPointsXi2::get_domain<GridTheta<2>>());
    const typename Patch2::IdxRange12 idx_range_etaxi2 =
            typename Patch2::IdxRange12(idx_range_eta2, idx_range_xi2);

    using SplineInterpPoints_Rg = SplineInterpPointsRg<Interpolation_v>;
    const IdxRange<GridRg> idx_range_r_g
            = IdxRange<GridRg>(SplineInterpPoints_Rg::template get_domain<GridRg>());
    const IdxRange<GridThetag> idx_range_theta_g
            = IdxRange<GridThetag>(SplineInterpPointsThetag::get_domain<GridThetag>());
    const IdxRange<GridRg, GridThetag> idx_range_rtheta_g
            = IdxRange<GridRg, GridThetag>(idx_range_r_g, idx_range_theta_g);


    // PATCH CONNECTION --------------------------------------------------------------------------
    static constexpr bool orientation_agree
            = ((std::is_same_v<Edge1, EastEdge1> && std::is_same_v<Edge2, WestEdge2>))
              || ((std::is_same_v<Edge1, WestEdge1> && std::is_same_v<Edge2, EastEdge2>));
    using Interface_1_2 = Interface<Edge1, Edge2, orientation_agree>;


    // TEST OPERATORS ----------------------------------------------------------------------------
    /**
     * @brief Check that the computed interface derivatives match with the equivalent global spline
     * derivatives for the exact formula. 
     * Check that error of the computed interface derivatives with respect to the equivalent global 
     * spline derivatives for the approximation are smaller than a given bound.
     * Test with different number of cells. 
     */
    void check_exact_and_approximation(
            int const n_cells,
            double const approximation_error_bound,
            host_t<DField<Patch1::IdxRange12>> const& function_1,
            host_t<DField<Patch2::IdxRange12>> const& function_2,
            SplineRThetagEvaluator const& evaluator_g,
            host_t<DField<IdxRange<BSplinesRg, BSplinesThetag>>> const& function_g_coef,
            EdgeTransformation<Interface_1_2> const& idx_convertor_12)
    {
        std::size_t n_points_1 = Kokkos::min(n_cells + 1, int(r1_ncells.value()));
        std::size_t n_points_2 = Kokkos::min(n_cells + 1, int(eta2_ncells.value()));
        // --- select the cells in patch 1
        Patch1::IdxRange1 reduced_idx_range_perp1;
        if constexpr (std::is_same_v<Edge1, EastEdge1>) {
            reduced_idx_range_perp1 = idx_range_r1.take_last(Patch1::IdxStep1(n_points_1));
        } else if (std::is_same_v<Edge1, WestEdge1>) {
            reduced_idx_range_perp1 = idx_range_r1.take_first(Patch1::IdxStep1(n_points_1));
        }

        // --- select the cells in patch 2
        using IdxRangePerp2 = std::conditional_t<
                std::is_same_v<Edge2, SouthEdge2>,
                Patch2::IdxRange2,
                Patch2::IdxRange1>;
        IdxRangePerp2 reduced_idx_range_perp2;
        if constexpr (std::is_same_v<Edge2, EastEdge2>) {
            reduced_idx_range_perp2 = idx_range_eta2.take_last(Patch2::IdxStep1(n_points_2));
        } else if constexpr (std::is_same_v<Edge2, WestEdge2>) {
            reduced_idx_range_perp2 = idx_range_eta2.take_first(Patch2::IdxStep1(n_points_2));
        } else if constexpr (std::is_same_v<Edge2, SouthEdge2>) {
            reduced_idx_range_perp2 = idx_range_xi2.take_first(Patch2::IdxStep2(n_points_2));
        }

        SingleInterfaceDerivativesCalculator<Interface_1_2> const
                derivatives_calculator(reduced_idx_range_perp1, reduced_idx_range_perp2);

        // Coefficients a and b
        double const coeff_deriv_patch_1 = derivatives_calculator.get_coeff_deriv_patch_1();
        double const coeff_deriv_patch_2 = derivatives_calculator.get_coeff_deriv_patch_2();

        EXPECT_EQ(
                coeff_deriv_patch_1,
                derivatives_calculator.template get_coeff_deriv_on_patch<Patch1>());
        EXPECT_EQ(
                coeff_deriv_patch_2,
                derivatives_calculator.template get_coeff_deriv_on_patch<Patch2>());

        using IdxPar2
                = std::conditional_t<std::is_same_v<Edge2, SouthEdge2>, Patch2::Idx1, Patch2::Idx2>;

        ddc::for_each(idx_range_theta1, [&](Patch1::Idx2 const& idx_par_1) {
            IdxPar2 idx_par_2 = idx_convertor_12(idx_par_1);

            // Coordinates to evaluate the derivative of the equivalent global spline.
            // --- coordinate of the boundary on the patch 1 (left of the interface)
            double coord_theta;
            Coord<Rg, Thetag> interface_minus_coord;
            if constexpr (std::is_same_v<Edge1, EastEdge1>) {
                coord_theta = double(ddc::coordinate(idx_par_1));
                interface_minus_coord = Coord<Rg, Thetag>(
                        double(ddc::coordinate(reduced_idx_range_perp1.front())),
                        coord_theta);
            } else if (std::is_same_v<Edge1, WestEdge1>) {
                coord_theta = double(theta1_max - ddc::coordinate(idx_par_1));
                interface_minus_coord = Coord<Rg, Thetag>(
                        double(r1_max - ddc::coordinate(reduced_idx_range_perp1.back())),
                        coord_theta);
            }

            // --- coordinate at the interface
            Coord<Rg, Thetag> interface_coord(double(r1_max), coord_theta);

            // --- coordinate of the boundary on the patch 2 (right of the interface)
            Coord<Rg, Thetag> interface_plus_coord;
            if constexpr (std::is_same_v<Edge2, EastEdge2>) {
                interface_plus_coord = Coord<Rg, Thetag>(
                        double(eta2_max - ddc::coordinate(reduced_idx_range_perp2.front())
                               + r1_max),
                        coord_theta);
            } else if (std::is_same_v<Edge2, WestEdge2>) {
                interface_plus_coord = Coord<Rg, Thetag>(
                        double(ddc::coordinate(reduced_idx_range_perp2.back())),
                        coord_theta);
            } else if (std::is_same_v<Edge2, SouthEdge2>) {
                interface_plus_coord = Coord<Rg, Thetag>(
                        double(ddc::coordinate(reduced_idx_range_perp2.back()) / xi2_max)
                                + double(r1_max),
                        coord_theta);
            }

            // Coefficient c
            double sum_values = derivatives_calculator.get_function_coefficients(
                    get_const_field(function_1[idx_par_1][reduced_idx_range_perp1]),
                    get_const_field(function_2[idx_par_2][reduced_idx_range_perp2]));

            double global_deriv
                    = evaluator_g.deriv_dim_1(interface_coord, get_const_field(function_g_coef));

            // Exact formula ---------------------------------------------------------------------
            double const deriv_patch_1
                    = evaluator_g
                              .deriv_dim_1(interface_minus_coord, get_const_field(function_g_coef));
            double const deriv_patch_2
                    = evaluator_g
                              .deriv_dim_1(interface_plus_coord, get_const_field(function_g_coef));
            double const local_deriv = sum_values + coeff_deriv_patch_1 * deriv_patch_1
                                       + coeff_deriv_patch_2 * deriv_patch_2;
            EXPECT_NEAR(local_deriv, global_deriv, 5e-13);

            // Approximation ---------------------------------------------------------------------
            EXPECT_NEAR(sum_values, global_deriv, approximation_error_bound);
        });
    };
};

} // end namespace


// Tuple of all the combinations of parameters we want to test.
template <typename... input_t>
using tuple_cat_t = decltype(std::tuple_cat(std::declval<input_t>()...));

using Cases = tuple_to_types_t<tuple_cat_t<
        cartesian_product_t<
                std::tuple<
                        std::integral_constant<ddc::BoundCond, ddc::BoundCond::HERMITE>
#if defined(NON_UNIFORM_MESH)
                        /*
                            Test with additional interpolation points as closure.
                            They are not especially Greville points. 
                            The mesh cannot be uniform to use the additional interpolation 
                            points as closure. 
                        */
                        ,
                        std::integral_constant<ddc::BoundCond, ddc::BoundCond::GREVILLE>
#endif
                        >,
                std::tuple<EastEdge1, WestEdge1>,
                std::tuple<EastEdge2, WestEdge2>>,
        std::tuple<std::tuple<
                std::integral_constant<ddc::BoundCond, ddc::BoundCond::HERMITE>,
                EastEdge1,
                SouthEdge2>>>>;


TYPED_TEST_SUITE(SingleInterfaceDerivativesCalculatorFixture, Cases);



// Check that the local grids and the equivalent global grid match together.
TYPED_TEST(SingleInterfaceDerivativesCalculatorFixture, InterpolationPointsCheck)
{
    // Get parameters of the test.
    constexpr ddc::BoundCond Interpolation_v = TestFixture::Interpolation_v;
    using Edge1 = typename TestFixture::Edge1;
    using Edge2 = typename TestFixture::Edge2;

    // --- Check with the patch 1 ----------------------------------------------------------------
    if constexpr (std::is_same_v<Edge1, EastEdge1>) {
        // Orientation first patch:   ↑→
        // Orientation global domain: ↑→
        ddc::for_each(TestFixture::idx_range_rtheta1, [&](Patch1::Idx12 const& idx) {
            Patch1::IdxStep1 idx_r(Patch1::Idx1(idx) - TestFixture::idx_range_r1.front());
            Patch1::IdxStep2 idx_theta(Patch1::Idx2(idx) - TestFixture::idx_range_theta1.front());
            Idx<GridRg, GridThetag> idx_g(idx_r.value(), idx_theta.value());
            EXPECT_NEAR(
                    ddc::coordinate(Patch1::Idx1(idx)),
                    ddc::coordinate(Idx<GridRg>(idx_g)),
                    1e-15);
            EXPECT_NEAR(
                    ddc::coordinate(Patch1::Idx2(idx)),
                    ddc::coordinate(Idx<GridThetag>(idx_g)),
                    1e-15);
        });
    } else if (std::is_same_v<Edge1, WestEdge1>) {
        // Orientation first patch:   ←↓
        // Orientation global domain: ↑→
        ddc::for_each(TestFixture::idx_range_rtheta1, [&](Patch1::Idx12 const& idx) {
            Patch1::IdxStep1 idx_r(Patch1::Idx1(idx) - TestFixture::idx_range_r1.front());
            Patch1::IdxStep2 idx_theta(Patch1::Idx2(idx) - TestFixture::idx_range_theta1.front());
            Idx<GridRg, GridThetag> idx_g;
            if constexpr (Interpolation_v == ddc::BoundCond::GREVILLE) {
                idx_g = Idx<GridRg, GridThetag>(
                        TestFixture::r1_ncells.value() + 1 - idx_r.value(),
                        TestFixture::theta1_ncells.value() - 1 - idx_theta.value());
            } else {
                idx_g = Idx<GridRg, GridThetag>(
                        TestFixture::r1_ncells.value() - idx_r.value(),
                        TestFixture::theta1_ncells.value() - 1 - idx_theta.value());
            }
            EXPECT_NEAR(
                    TestFixture::r1_max - ddc::coordinate(Patch1::Idx1(idx)),
                    ddc::coordinate(Idx<GridRg>(idx_g)),
                    1e-15);
            EXPECT_NEAR(
                    ddc::coordinate(Patch1::Idx2(idx)),
                    ddc::coordinate(TestFixture::idx_range_theta_g.back())
                            - ddc::coordinate(Idx<GridThetag>(idx_g)),
                    1e-15);
        });
    }

    // --- Check with the patch 2 ----------------------------------------------------------------
    if constexpr (std::is_same_v<Edge2, EastEdge2>) {
        // Orientation second patch:  ←↓
        // Orientation global domain: ↑→
        ddc::for_each(TestFixture::idx_range_etaxi2, [&](Patch2::Idx12 const& idx) {
            Patch2::IdxStep1 idx_eta(Patch2::Idx1(idx) - TestFixture::idx_range_eta2.front());
            Patch2::IdxStep2 idx_xi(Patch2::Idx2(idx) - TestFixture::idx_range_xi2.front());
            Idx<GridRg, GridThetag> idx_g;
            if constexpr (Interpolation_v == ddc::BoundCond::GREVILLE) {
                idx_g = Idx<GridRg, GridThetag>(
                        TestFixture::eta2_ncells.value() - idx_eta.value()
                                + TestFixture::r1_ncells.value() + 2,
                        idx_xi.value());
            } else {
                idx_g = Idx<GridRg, GridThetag>(
                        TestFixture::eta2_ncells.value() - idx_eta.value()
                                + TestFixture::r1_ncells.value(),
                        idx_xi.value());
            }
            EXPECT_NEAR(
                    TestFixture::eta2_max - ddc::coordinate(Patch2::Idx1(idx))
                            + TestFixture::r1_max,
                    ddc::coordinate(Idx<GridRg>(idx_g)),
                    5e-15);
            EXPECT_NEAR(
                    ddc::coordinate(Patch2::Idx2(idx)),
                    ddc::coordinate(Idx<GridThetag>(idx_g)),
                    1e-15);
        });
    } else if (std ::is_same_v<Edge2, WestEdge2>) {
        // Orientation second patch:  ↑→
        // Orientation global domain: ↑→
        ddc::for_each(TestFixture::idx_range_etaxi2, [&](Patch2::Idx12 const& idx) {
            Patch2::IdxStep1 idx_eta(Patch2::Idx1(idx) - TestFixture::idx_range_eta2.front());
            Patch2::IdxStep2 idx_xi(Patch2::Idx2(idx) - TestFixture::idx_range_xi2.front());
            Idx<GridRg, GridThetag> idx_g;
            if constexpr (Interpolation_v == ddc::BoundCond::GREVILLE) {
                idx_g = Idx<GridRg, GridThetag>(
                        idx_eta.value() + TestFixture::r1_ncells.value() + 1,
                        idx_xi.value());
            } else {
                idx_g = Idx<GridRg, GridThetag>(
                        idx_eta.value() + TestFixture::r1_ncells.value(),
                        idx_xi.value());
            }
            EXPECT_NEAR(
                    ddc::coordinate(Patch2::Idx1(idx)),
                    ddc::coordinate(Idx<GridRg>(idx_g)),
                    5e-15);
            EXPECT_NEAR(
                    ddc::coordinate(Patch2::Idx2(idx)),
                    ddc::coordinate(Idx<GridThetag>(idx_g)),
                    1e-15);
        });
    } else if (std ::is_same_v<Edge2, SouthEdge2>) {
        // Orientation second patch:  ↓→
        // Orientation global domain: ↑→
        if constexpr (Interpolation_v == ddc::BoundCond::GREVILLE) {
            Kokkos::abort("Test case not implemented."
                          "We cannot make the Edge1 on Theta with N points to match with the "
                          "Edge2 on R with N+1 points.");
        }
        ddc::for_each(TestFixture::idx_range_etaxi2, [&](Patch2::Idx12 const& idx) {
            Patch2::IdxStep1 idx_eta(Patch2::Idx1(idx) - TestFixture::idx_range_eta2.front());
            Patch2::IdxStep2 idx_xi(Patch2::Idx2(idx) - TestFixture::idx_range_xi2.front());
            Idx<GridRg, GridThetag>
                    idx_g(idx_xi.value() + TestFixture::r1_ncells.value(),
                          TestFixture::eta2_ncells.value() - idx_eta.value());
            EXPECT_NEAR(
                    double(ddc::coordinate(Patch2::Idx2(idx))) + double(TestFixture::r1_max),
                    ddc::coordinate(Idx<GridRg>(idx_g)),
                    5e-15);
            // On periodic mesh, the last interpolation point is not the last break point.
            if (idx_eta.value() != 0) {
                EXPECT_NEAR(
                        double(TestFixture::eta2_max - ddc::coordinate(Patch2::Idx1(idx))),
                        ddc::coordinate(Idx<GridThetag>(idx_g)),
                        1e-14);
            }
        });
    }
}


// Check the values of the computed interface derivatives.
TYPED_TEST(
        SingleInterfaceDerivativesCalculatorFixture,
        InterfaceDerivativesExactAndApproximationFormulae)
{
    // Get parameters of the test.
    constexpr ddc::BoundCond Interpolation_v = TestFixture::Interpolation_v;
    using Edge1 = typename TestFixture::Edge1;
    using Edge2 = typename TestFixture::Edge2;
    using Interface_1_2 = typename TestFixture::Interface_1_2;

    // Define EdgeTransformation operator.
    using IdxRangePar2 = std::conditional_t<
            (std::is_same_v<Edge2, SouthEdge2>),
            typename Patch2::IdxRange1,
            typename Patch2::IdxRange2>;
    IdxRangePar2 idx_range_par2;
    if constexpr (std::is_same_v<Edge2, SouthEdge2>) {
        idx_range_par2 = TestFixture::idx_range_eta2;
    } else {
        idx_range_par2 = TestFixture::idx_range_xi2;
    }

    EdgeTransformation<Interface_1_2>
            idx_convertor_12(TestFixture::idx_range_theta1, idx_range_par2);

    // Initialise functions values ===============================================================
    // --- patch 1
    host_t<DFieldMem<Patch1::IdxRange12>> function_1_alloc(TestFixture::idx_range_rtheta1);
    host_t<DField<Patch1::IdxRange12>> function_1 = get_field(function_1_alloc);

    // --- patch 2
    host_t<DFieldMem<Patch2::IdxRange12>> function_2_alloc(TestFixture::idx_range_etaxi2);
    host_t<DField<Patch2::IdxRange12>> function_2 = get_field(function_2_alloc);

    // --- global
    host_t<DFieldMem<IdxRange<GridRg, GridThetag>>> function_g_alloc(
            TestFixture::idx_range_rtheta_g);
    host_t<DField<IdxRange<GridRg, GridThetag>>> function_g = get_field(function_g_alloc);

    // Fill in with the correct value.
    // --- patch 1
    if constexpr (std::is_same_v<Edge1, EastEdge1>) {
        // Same orientation than the equivalent global domain.
        initialise_2D_function<Patch1::Grid1, Patch1::Grid2>(function_1);
    } else {
        // Different orientation: global: ↑→ | local: ←↓
        ddc::for_each(get_idx_range(function_1), [&](Idx<Patch1::Grid1, Patch1::Grid2> idx) {
            // Get the coordinate on the equivalent global domain.
            double const rg = TestFixture::r1_max - ddc::coordinate(Idx<Patch1::Grid1>(idx));
            double const thetag
                    = TestFixture::theta1_max - ddc::coordinate(Idx<Patch1::Grid2>(idx));
            function_1(idx) = rg * (3. - rg) * Kokkos::sin(thetag);
        });
    }

    // --- patch 2
    if constexpr (std::is_same_v<Edge2, WestEdge2>) {
        // Same orientation than the equivalent global domain.
        initialise_2D_function<Patch2::Grid1, Patch2::Grid2>(function_2);
    } else if (std::is_same_v<Edge2, EastEdge2>) {
        // Different orientation: global: ↑→ | local: ←↓
        ddc::for_each(get_idx_range(function_2), [&](Idx<Patch2::Grid1, Patch2::Grid2> idx) {
            // Get the coordinate on the equivalent global domain.
            double const rg = TestFixture::eta2_max - ddc::coordinate(Idx<Patch2::Grid1>(idx))
                              + TestFixture::r1_max;
            double const thetag = TestFixture::xi2_max - ddc::coordinate(Idx<Patch2::Grid2>(idx));
            function_2(idx) = rg * (3. - rg) * Kokkos::sin(thetag);
        });
    } else if (std::is_same_v<Edge2, SouthEdge2>) {
        // Different orientation: global: ↑→ | local: ↓→
        ddc::for_each(get_idx_range(function_2), [&](Idx<Patch2::Grid1, Patch2::Grid2> idx) {
            // Get the coordinate on the equivalent global domain.
            double const rg
                    = double(ddc::coordinate(Idx<Patch2::Grid2>(idx))) + TestFixture::r1_max;
            double const thetag
                    = double((TestFixture::eta2_max - ddc::coordinate(Idx<Patch2::Grid1>(idx))));
            function_2(idx) = rg * (3. - rg) * Kokkos::sin(thetag);
        });
    }

    // --- global
    initialise_2D_function<GridRg, GridThetag>(function_g);


    // Build an equivalent global spline =========================================================
    SplineRThetagBuilder<Interpolation_v> builder_g(TestFixture::idx_range_rtheta_g);

    host_t<DFieldMem<IdxRange<BSplinesRg, BSplinesThetag>>> function_g_coef_alloc(
            builder_g.batched_spline_domain(TestFixture::idx_range_rtheta_g));
    host_t<DField<IdxRange<BSplinesRg, BSplinesThetag>>> function_g_coef
            = get_field(function_g_coef_alloc);

    if constexpr (Interpolation_v == ddc::BoundCond::GREVILLE) {
        // --- Spline builder
        builder_g(function_g_coef, get_const_field(function_g));
    } else {
        // --- Set derivatives
        Idx<ddc::Deriv<Rg>> first_deriv_rg(1);
        IdxStep<ddc::Deriv<Rg>> n_deriv_rg(1);
        IdxRange<ddc::Deriv<Rg>> deriv_rg_idx_range(first_deriv_rg, n_deriv_rg);

        IdxRange<ddc::Deriv<Rg>, GridThetag>
                derivs_rg_idx_range(deriv_rg_idx_range, TestFixture::idx_range_theta_g);

        host_t<DFieldMem<IdxRange<ddc::Deriv<Rg>, GridThetag>>> derivs_rgmin_alloc(
                derivs_rg_idx_range);
        host_t<DField<IdxRange<ddc::Deriv<Rg>, GridThetag>>> derivs_rgmin
                = get_field(derivs_rgmin_alloc);

        host_t<DFieldMem<IdxRange<ddc::Deriv<Rg>, GridThetag>>> derivs_rgmax_alloc(
                derivs_rg_idx_range);
        host_t<DField<IdxRange<ddc::Deriv<Rg>, GridThetag>>> derivs_rgmax
                = get_field(derivs_rgmax_alloc);

        ddc::for_each(TestFixture::idx_range_theta_g, [&](Idx<GridThetag> const& idx_thetag) {
            derivs_rgmin(first_deriv_rg, idx_thetag) = 3 * Kokkos::sin(ddc::coordinate(idx_thetag));
            derivs_rgmax(first_deriv_rg, idx_thetag)
                    = -3 * Kokkos::sin(ddc::coordinate(idx_thetag));
        });

        // --- Spline builder
        builder_g(
                function_g_coef,
                get_const_field(function_g),
                std::optional(get_const_field(derivs_rgmin)),
                std::optional(get_const_field(derivs_rgmax)));
    }

    ddc::ConstantExtrapolationRule<Rg, Thetag> bc_rmin_g(TestFixture::rg_min);
    ddc::ConstantExtrapolationRule<Rg, Thetag> bc_rmax_g(TestFixture::rg_max);
    ddc::PeriodicExtrapolationRule<Thetag> bc_theta_g;
    SplineRThetagEvaluator evaluator_g(bc_rmin_g, bc_rmax_g, bc_theta_g, bc_theta_g);


    // Check the local derivatives with the global ones at the interfaces ========================
    if constexpr (Interpolation_v == ddc::BoundCond::GREVILLE) {
        // We test if the boundaries are well treated => only work with 5 cells to better identify an error.
        // 5 cells -------------------------------------------------------------------------------
        SingleInterfaceDerivativesCalculator<
                Interface_1_2,
                ddc::BoundCond::GREVILLE,
                ddc::BoundCond::GREVILLE> const
                derivatives_calculator(
                        TestFixture::idx_range_rtheta1,
                        TestFixture::idx_range_etaxi2);


        ddc::for_each(TestFixture::idx_range_theta1, [&](Patch1::Idx2 const& idx2_1) {
            Patch2::Idx2 idx2_2 = idx_convertor_12(idx2_1);

            // Coordinate at the interface.
            Coord<Rg, Thetag> interface_coord;
            if constexpr (std::is_same_v<Edge1, EastEdge1>) {
                interface_coord = Coord<
                        Rg,
                        Thetag>(double(TestFixture::r1_max), double(ddc::coordinate(idx2_1)));
            } else {
                interface_coord = Coord<Rg, Thetag>(
                        double(TestFixture::r1_max),
                        double(TestFixture::theta1_max - ddc::coordinate(idx2_1)));
            }

            // Coefficient c.
            double const sum_values = derivatives_calculator.get_function_coefficients(
                    get_const_field(function_2[idx2_2]),
                    get_const_field(function_1[idx2_1]));

            double const global_deriv
                    = evaluator_g.deriv_dim_1(interface_coord, get_const_field(function_g_coef));

            // Exact formula ---------------------------------------------------------------------
            double const local_deriv = sum_values;
            EXPECT_NEAR(local_deriv, global_deriv, 5e-13);
        });
    } else {
        // 30 cells ------------------------------------------------------------------------------
        TestFixture::check_exact_and_approximation(
                30,
                5e-13,
                function_1,
                function_2,
                evaluator_g,
                function_g_coef,
                idx_convertor_12);

        // 20 cells ------------------------------------------------------------------------------
        TestFixture::check_exact_and_approximation(
                20,
                1e-10,
                function_1,
                function_2,
                evaluator_g,
                function_g_coef,
                idx_convertor_12);

        // 15 cells ------------------------------------------------------------------------------
        TestFixture::check_exact_and_approximation(
                15,
                1e-7,
                function_1,
                function_2,
                evaluator_g,
                function_g_coef,
                idx_convertor_12);

        // 10 cells ------------------------------------------------------------------------------
        TestFixture::check_exact_and_approximation(
                10,
                5e-5,
                function_1,
                function_2,
                evaluator_g,
                function_g_coef,
                idx_convertor_12);

        // 5 cells -------------------------------------------------------------------------------
        TestFixture::check_exact_and_approximation(
                5,
                2e-2,
                function_1,
                function_2,
                evaluator_g,
                function_g_coef,
                idx_convertor_12);
    }
}