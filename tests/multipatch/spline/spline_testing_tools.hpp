// SPDX-License-Identifier: MIT
#pragma once

#include <ddc/ddc.hpp>
#include <ddc/kernels/splines.hpp>

#include <sll/mapping/cartesian_to_circular.hpp>
#include <sll/mapping/circular_to_cartesian.hpp>

#include <gtest/gtest.h>

#include "2patches_2d_onion_shape_non_uniform.hpp"
#include "ddc_helper.hpp"
#include "mesh_builder.hpp"
#include "multipatch_field.hpp"
#include "multipatch_type.hpp"
#include "onion_patch_locator.hpp"
#include "physical_geometry.hpp"
#include "types.hpp"



using namespace onion_shape_non_uniform_2d_2patches;
using namespace physical_geometry;

using HostExecSpace = Kokkos::DefaultHostExecutionSpace;
using DeviceExecSpace = Kokkos::DefaultExecutionSpace;

template <int PatchIdx>
using SplineInterpPointsR = ddc::GrevilleInterpolationPoints<
        BSplinesR<PatchIdx>,
        ddc::BoundCond::GREVILLE,
        ddc::BoundCond::GREVILLE>;
template <int PatchIdx>
using SplineInterpPointsTheta = ddc::GrevilleInterpolationPoints<
        BSplinesTheta<PatchIdx>,
        ddc::BoundCond::PERIODIC,
        ddc::BoundCond::PERIODIC>;


template <int PatchIdx, class ExecSpace>
using SplineRThetaBuilder = ddc::SplineBuilder2D<
        ExecSpace,
        typename ExecSpace::memory_space,
        BSplinesR<PatchIdx>,
        BSplinesTheta<PatchIdx>,
        GridR<PatchIdx>,
        GridTheta<PatchIdx>,
        ddc::BoundCond::GREVILLE,
        ddc::BoundCond::GREVILLE,
        ddc::BoundCond::PERIODIC,
        ddc::BoundCond::PERIODIC,
        ddc::SplineSolver::LAPACK,
        GridR<PatchIdx>,
        GridTheta<PatchIdx>>;


template <int PatchIdx, class ExecSpace>
using SplineRThetaEvaluator = ddc::SplineEvaluator2D<
        ExecSpace,
        typename ExecSpace::memory_space,
        BSplinesR<PatchIdx>,
        BSplinesTheta<PatchIdx>,
        GridR<PatchIdx>,
        GridTheta<PatchIdx>,
        ddc::ConstantExtrapolationRule<R, Theta>,
        ddc::ConstantExtrapolationRule<R, Theta>,
        ddc::PeriodicExtrapolationRule<Theta>,
        ddc::PeriodicExtrapolationRule<Theta>,
        GridR<PatchIdx>,
        GridTheta<PatchIdx>>;


using LogicalToPhysicalMapping = CircularToCartesian<R, Theta, X, Y>;
using PhysicalToLogicalMapping = CartesianToCircular<X, Y, R, Theta>;
using MultipatchIdxRange = MultipatchType<IdxRangeOnPatch, Patch1, Patch2>;

template <class ExecSpace>
using PatchLocator = OnionPatchLocator<
        MultipatchIdxRange,
        LogicalToPhysicalMapping,
        PhysicalToLogicalMapping,
        ExecSpace>;


/**
 * @brief Google test class to define index ranges and splines 
 * on an Onion shape domain with 2 patches. 
 */
class MultipatchSplineOnionShapeTest : public ::testing::Test
{
protected:
    /// @brief Nomber of cells in the R direction on the Patch1.
    static constexpr Patch1::IdxStep1 r1_ncells = Patch1::IdxStep1(16);
    /// @brief Nomber of cells in the Theta direction on the Patch1.
    static constexpr Patch1::IdxStep2 theta1_ncells = Patch1::IdxStep2(10);

    /// @brief Nomber of cells in the R direction on the Patch2.
    static constexpr Patch2::IdxStep1 r2_ncells = Patch2::IdxStep1(8);
    /// @brief Nomber of cells in the Theta direction on the Patch2.
    static constexpr Patch2::IdxStep2 theta2_ncells = Patch2::IdxStep2(12);

    // Coordinates delimiting the patches
    /// @brief Minimum value of R of the domain on the Patch1.
    static constexpr typename Patch1::Coord1 r1_min = Patch1::Coord1(0.25);
    /// @brief Maximum value of R of the domain on the Patch1.
    static constexpr typename Patch1::Coord1 r1_max = Patch1::Coord1(1.0);

    /// @brief Minimum value of Theta of the domain on the Patch1.
    static constexpr typename Patch1::Coord2 theta1_min = Patch1::Coord2(0.0);
    /// @brief Maximum value of Theta of the domain on the Patch1.
    static constexpr typename Patch1::Coord2 theta1_max = Patch1::Coord2(2 * M_PI);

    /// @brief Minimum value of R of the domain on the Patch2.
    static constexpr typename Patch2::Coord1 r2_min = Patch2::Coord1(1.0);
    /// @brief Maximum value of R of the domain on the Patch2.
    static constexpr typename Patch2::Coord1 r2_max = Patch2::Coord1(2.0);

    /// @brief Minimum value of Theta of the domain on the Patch2.
    static constexpr Patch2::Coord2 theta2_min = Patch2::Coord2(0.0);
    /// @brief Maximum value of Theta of the domain on the Patch2.
    static constexpr Patch2::Coord2 theta2_max = Patch2::Coord2(2 * M_PI);

    // Index ranges
    /// @brief 2D index range for the grids on Patch1.
    Patch1::IdxRange12 const idx_range_rtheta1;
    /// @brief 2D index range for the B-splines on Patch1.
    Patch1::IdxRangeBS12 const spline_idx_range_rtheta1;
    /// @brief 2D index range for the grids on Patch2.
    Patch2::IdxRange12 const idx_range_rtheta2;
    /// @brief 2D index range for the B-splines on Patch2.
    Patch2::IdxRangeBS12 const spline_idx_range_rtheta2;

    // Spline representations
    /// @brief Memory allocation of the spline coefficients of the function on Patch1.
    DFieldMem<IdxRange<BSplinesR<1>, BSplinesTheta<1>>> function_1_coef_alloc;
    /// @brief Memory allocation of the spline coefficients of the function on Patch2.
    DFieldMem<IdxRange<BSplinesR<2>, BSplinesTheta<2>>> function_2_coef_alloc;

    /// @brief MultipatchType of const fields of the spline coefficients of the function
    // on the patches.
    MultipatchField<ConstSplineCoeffOnPatch_2D, Patch1, Patch2> const splines;

    /// @brief MultipatchType of index ranges on the patches.
    MultipatchIdxRange const all_idx_ranges;


public:
    MultipatchSplineOnionShapeTest()
        : idx_range_rtheta1(
                SplineInterpPointsR<1>::get_domain<GridR<1>>(),
                SplineInterpPointsTheta<1>::get_domain<GridTheta<1>>())
        , spline_idx_range_rtheta1(
                  ddc::discrete_space<BSplinesR<1>>().full_domain(),
                  ddc::discrete_space<BSplinesTheta<1>>().full_domain())
        , idx_range_rtheta2(
                  SplineInterpPointsR<2>::get_domain<GridR<2>>(),
                  SplineInterpPointsTheta<2>::get_domain<GridTheta<2>>())
        , spline_idx_range_rtheta2(
                  ddc::discrete_space<BSplinesR<2>>().full_domain(),
                  ddc::discrete_space<BSplinesTheta<2>>().full_domain())
        , function_1_coef_alloc(spline_idx_range_rtheta1)
        , function_2_coef_alloc(spline_idx_range_rtheta2)
        , splines(get_const_field(set_spline_1(get_field(function_1_coef_alloc))),
                  get_const_field(set_spline_2(get_field(function_2_coef_alloc))))
        , all_idx_ranges(idx_range_rtheta1, idx_range_rtheta2)
    {
    }

    /// @brief Initialisation the discrete spaces on the patches.
    static void SetUpTestSuite()
    {
        // Creating of meshes and supports ...........................................................
        // Patch 1
        ddc::init_discrete_space<BSplinesR<1>>(
                build_uniform_break_points(r1_min, r1_max, r1_ncells));
        ddc::init_discrete_space<BSplinesTheta<1>>(
                build_uniform_break_points(theta1_min, theta1_max, theta1_ncells));

        ddc::init_discrete_space<GridR<1>>(SplineInterpPointsR<1>::get_sampling<GridR<1>>());
        ddc::init_discrete_space<GridTheta<1>>(
                SplineInterpPointsTheta<1>::get_sampling<GridTheta<1>>());

        // Patch 2
        ddc::init_discrete_space<BSplinesR<2>>(
                build_uniform_break_points(r2_min, r2_max, r2_ncells));
        ddc::init_discrete_space<BSplinesTheta<2>>(
                build_uniform_break_points(theta2_min, theta2_max, theta2_ncells));

        ddc::init_discrete_space<GridR<2>>(SplineInterpPointsR<2>::get_sampling<GridR<2>>());
        ddc::init_discrete_space<GridTheta<2>>(
                SplineInterpPointsTheta<2>::get_sampling<GridTheta<2>>());
    }


    /**
     * @brief Initialise the function on Patch1. 
     * Here the function is @f$ f(r,\theta)= r \sin(\theta)@f$.
     * @param[inout] function_1_coef Field of function value on Patch1.
     * @return A field of function value initiliased on Patch1.
     */
    DField<Patch1::IdxRangeBS12> set_spline_1(DField<Patch1::IdxRangeBS12> const& function_1_coef);

    /**
     * @brief Initialise the function on Patch2.
     * Here the function is @f$ f(r,\theta)= r^2 \sin(\theta)@f$.
     * @param[inout] function_2_coef Field of function value on Patch2.
     * @return A field of function value initiliased on Patch2.
     */
    DField<Patch2::IdxRangeBS12> set_spline_2(DField<Patch2::IdxRangeBS12> const& function_2_coef);
};
