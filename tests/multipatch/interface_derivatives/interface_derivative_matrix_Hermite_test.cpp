// SPDX-License-Identifier: MIT

#include <ddc/ddc.hpp>
#include <ddc/kernels/splines.hpp>

#include <gtest/gtest.h>

#include "3patches_2d_non_periodic_non_uniform.hpp"
#include "ddc_alias_inline_functions.hpp"
#include "ddc_helper.hpp"
#include "derivative_field.hpp"
#include "derivative_field_mem.hpp"
#include "interface.hpp"
#include "interface_derivative_matrix.hpp"
#include "interface_derivatives_matrix_test_utils.hpp"
#include "interface_derivatives_test_utils.hpp"
#include "linear_coord_transform.hpp"
#include "mesh_builder.hpp"
#include "non_uniform_interpolation_points.hpp"
#include "orthogonal_coord_transforms.hpp"
#include "single_interface_derivatives_calculator.hpp"
#include "single_interface_derivatives_calculator_collection.hpp"
#include "types.hpp"
#include "view.hpp"


/*
    Test InterfaceDerivativeMatrix on the following geometry:

    #if defined(REVERSE_PATCH1)
        |  1  |  2  |  3  | 
           ←↓   ↑→     ↑→

    #elif defined(REVERSE_PATCH2)
        |  1  |  2  |  3  | 
          ↑→     ←↓   ↑→

    #elif defined(REVERSE_PATCH3)
        |  1  |  2  |  3  | 
           ↑→    ↑→    ←↓  
    
    #elif defined(CHANGE_BOUND1)
        |  1  |  2  |  3  | 
           ←↑    ↑→    ↑→  
           
    #elif defined(CHANGE_BOUND3)
        |  1  |  2  |  3  | 
           ↑→    ↑→    ↓→  

        with the global X dimension with Hermite boundary conditions 
        and the global Y spline with additional points as closure condition 
        (ddc::BoundCond::GREVILLE).

    > test ddc::BoundCond::HERMITE boundary conditions. 
    > test application on the X direction. 
    > test application to compute first derivatives and cross-derivatives. 
    > test on a non-uniform patches. 
    > test with a middle ill-oriented interface. 
    > test with not matching directions of the patches.
    > test with approximation formulation in SingleInterfaceDerivativeCalculator. 
    > test agreement between computed and global spline derivatives. 
    > test agreement between local and global splines.
*/

namespace {
// Multi-patch tags ---
using namespace non_periodic_non_uniform_2d_3patches;

//enum TestCase { REVERSE_PATCH1, REVERSE_PATCH2, REVERSE_PATCH13, CHANGE_BOUND1, CHANGE_BOUND3 };

// Equivalent global mesh tags ---
struct Xg
{
    static bool constexpr PERIODIC = false;
    static bool constexpr IS_COVARIANT = true;
    static bool constexpr IS_CONTRAVARIANT = true;
    using Dual = Xg;
};
struct Yg
{
    static bool constexpr PERIODIC = false;
    static bool constexpr IS_COVARIANT = true;
    static bool constexpr IS_CONTRAVARIANT = true;
    using Dual = Yg;
};

struct GridXg : NonUniformGridBase<Xg>
{
};
struct GridYg : NonUniformGridBase<Yg>
{
};

struct BSplinesXg : ddc::NonUniformBSplines<Xg, 3>
{
};
struct BSplinesYg : ddc::NonUniformBSplines<Yg, 3>
{
};

using DerivXg = ddc::Deriv<Xg>;
using DerivYg = ddc::Deriv<Yg>;

using HostExecSpace = Kokkos::DefaultHostExecutionSpace;

// Interpolation points type for the patches.
template <std::size_t PatchIdx>
using SplineInterpPointsX = ddcHelper::NonUniformInterpolationPoints<
        BSplinesX<PatchIdx>,
        ddc::BoundCond::HERMITE,
        ddc::BoundCond::HERMITE>;

template <std::size_t PatchIdx>
using SplineInterpPointsY = ddcHelper::NonUniformInterpolationPoints<
        BSplinesY<PatchIdx>,
        ddc::BoundCond::HERMITE,
        ddc::BoundCond::HERMITE>;

// Interpolation points type for the equivalent global spline.
using SplineInterpPointsXg = ddcHelper::
        NonUniformInterpolationPoints<BSplinesXg, ddc::BoundCond::HERMITE, ddc::BoundCond::HERMITE>;
using SplineInterpPointsYg = ddcHelper::
        NonUniformInterpolationPoints<BSplinesYg, ddc::BoundCond::HERMITE, ddc::BoundCond::HERMITE>;

// Operators on the equivalent global spline.
using SplineRThetagBuilder = ddc::SplineBuilder2D<
        HostExecSpace,
        typename HostExecSpace::memory_space,
        BSplinesXg,
        BSplinesYg,
        GridXg,
        GridYg,
        ddc::BoundCond::HERMITE,
        ddc::BoundCond::HERMITE,
        ddc::BoundCond::HERMITE,
        ddc::BoundCond::HERMITE,
        ddc::SplineSolver::LAPACK>;

using SplineRThetagBuilderDerivField = SplineBuliderDerivField2D<
        HostExecSpace,
        BSplinesXg,
        BSplinesYg,
        GridXg,
        GridYg,
        ddc::BoundCond::HERMITE,
        ddc::BoundCond::HERMITE,
        ddc::BoundCond::HERMITE,
        ddc::BoundCond::HERMITE>;

using SplineRThetagEvaluator = ddc::SplineEvaluator2D<
        HostExecSpace,
        typename HostExecSpace::memory_space,
        BSplinesXg,
        BSplinesYg,
        GridXg,
        GridYg,
        ddc::ConstantExtrapolationRule<Xg, Yg>,
        ddc::ConstantExtrapolationRule<Xg, Yg>,
        ddc::ConstantExtrapolationRule<Yg, Xg>,
        ddc::ConstantExtrapolationRule<Yg, Xg>>;

// INTERFACES ------------------------------------------------------------------------------------
using NorthInterface1 = Interface<NorthEdge<1>, OutsideEdge, true>;
using NorthInterface2 = Interface<NorthEdge<2>, OutsideEdge, true>;
using NorthInterface3 = Interface<NorthEdge<3>, OutsideEdge, true>;

using SouthInterface1 = Interface<OutsideEdge, SouthEdge<1>, true>;
using SouthInterface2 = Interface<OutsideEdge, SouthEdge<2>, true>;
using SouthInterface3 = Interface<OutsideEdge, SouthEdge<3>, true>;

using EastInterface1 = Interface<OutsideEdge, EastEdge<1>, true>;
using EastInterface3 = Interface<OutsideEdge, EastEdge<3>, true>;

using WestInterface1 = Interface<OutsideEdge, WestEdge<1>, true>;
using WestInterface3 = Interface<OutsideEdge, WestEdge<3>, true>;

using AllGrids = ddc::detail::TypeSeq<GridX<1>, GridX<2>, GridX<3>, GridY<1>, GridY<2>, GridY<3>>;
using AllBSpls = ddc::detail::
        TypeSeq<BSplinesX<1>, BSplinesX<2>, BSplinesX<3>, BSplinesY<1>, BSplinesY<2>, BSplinesY<3>>;

template <int I>
struct MatchingPatchTransform
{
    using XTransform = LinearCoordTransform<Xg, X<I>>;
    using YTransform = LinearCoordTransform<Yg, Y<I>>;
    XTransform x_transform;
    YTransform y_transform;
    OrthogonalCoordTransforms<
            Coord<Xg, Yg>,
            Coord<X<I>, Y<I>>,
            Coord<Xg, Yg>,
            XTransform,
            YTransform>
            coord_transform;
    MatchingPatchTransform(XTransform x_transform, YTransform y_transform)
        : x_transform(x_transform)
        , y_transform(y_transform)
        , coord_transform(x_transform, y_transform)
    {
    }
};

template <int I>
struct AlignedPatchTransform : public MatchingPatchTransform<I>
{
    AlignedPatchTransform(Coord<Xg, Yg> origin)
        : MatchingPatchTransform<I>(
                typename MatchingPatchTransform<I>::
                        XTransform(Coord<Xg>(origin), convert_dim<X<I>>(Coord<Xg>(origin)), 1.0),
                typename MatchingPatchTransform<I>::
                        YTransform(Coord<Yg>(origin), convert_dim<Y<I>>(Coord<Yg>(origin)), 1.0))
    {
    }
};

template <int I>
struct ReversePatchTransform : public MatchingPatchTransform<I>
{
    ReversePatchTransform(Coord<Xg, Yg> origin)
        : MatchingPatchTransform<I>(
                typename MatchingPatchTransform<I>::XTransform(
                        Coord<Xg>(origin),
                        convert_dim<X<I>>(Coord<Xg>(origin) + 1.0),
                        -1.0),
                typename MatchingPatchTransform<I>::YTransform(
                        Coord<Yg>(origin),
                        convert_dim<Y<I>>(Coord<Yg>(origin) + 1.0),
                        -1.0))
    {
    }
};

template <int I>
struct BoundChangePatchTransform
{
    using XTransform = LinearCoordTransform<Xg, Y<I>>;
    using YTransform = LinearCoordTransform<Yg, X<I>>;
    XTransform x_transform;
    YTransform y_transform;
    OrthogonalCoordTransforms<
            Coord<Yg, Xg>,
            Coord<X<I>, Y<I>>,
            Coord<Xg, Yg>,
            XTransform,
            YTransform>
            coord_transform;
    BoundChangePatchTransform(XTransform x_transform, YTransform y_transform)
        : x_transform(x_transform)
        , y_transform(y_transform)
        , coord_transform(x_transform, y_transform)
    {
    }
};

template <int I>
struct BoundChangeXPatchTransform : public BoundChangePatchTransform<I>
{
    BoundChangeXPatchTransform(Coord<Xg, Yg> origin)
        : BoundChangePatchTransform<I>(
                typename BoundChangePatchTransform<I>::XTransform(
                        Coord<Xg>(origin),
                        convert_dim<Y<I>>(Coord<Xg>(origin) + 1.0),
                        -1.0),
                typename BoundChangePatchTransform<I>::
                        YTransform(Coord<Yg>(origin), convert_dim<X<I>>(Coord<Yg>(origin)), -1.0))
    {
    }
};

template <int I>
struct BoundChangeYPatchTransform : public BoundChangePatchTransform<I>
{
    BoundChangeYPatchTransform(Coord<Xg, Yg> origin)
        : BoundChangePatchTransform<I>(
                typename BoundChangePatchTransform<I>::
                        XTransform(Coord<Xg>(origin), convert_dim<Y<I>>(Coord<Xg>(origin)), -1.0),
                typename BoundChangePatchTransform<I>::YTransform(
                        Coord<Yg>(origin),
                        convert_dim<X<I>>(Coord<Yg>(origin) + 1.0),
                        -1.0))
    {
    }
};

struct RevPatch1
{
    using Interface_1_2 = Interface<WestEdge<1>, WestEdge<2>, false>;
    using Interface_2_3 = Interface<EastEdge<2>, WestEdge<3>, true>;

    using OutsideInterface1 = Interface<OutsideEdge, EastEdge<1>, true>;
    using OutsideInterface3 = Interface<EastEdge<3>, OutsideEdge, true>;

    using CoordTransform1 = ReversePatchTransform<1>;
    using CoordTransform2 = AlignedPatchTransform<2>;
    using CoordTransform3 = AlignedPatchTransform<3>;

    using Connectivity = MultipatchConnectivity<
            NorthInterface1,
            SouthInterface1,
            NorthInterface2,
            SouthInterface2,
            NorthInterface3,
            SouthInterface3,
            OutsideInterface1,
            OutsideInterface3,
            Interface_1_2,
            Interface_2_3>;
};

struct RevPatch2
{
    using Interface_1_2 = Interface<EastEdge<1>, EastEdge<2>, false>;
    using Interface_2_3 = Interface<WestEdge<2>, WestEdge<3>, false>;

    using OutsideInterface1 = Interface<OutsideEdge, WestEdge<1>, true>;
    using OutsideInterface3 = Interface<EastEdge<3>, OutsideEdge, true>;

    using CoordTransform1 = AlignedPatchTransform<1>;
    using CoordTransform2 = ReversePatchTransform<2>;
    using CoordTransform3 = AlignedPatchTransform<3>;

    using Connectivity = MultipatchConnectivity<
            NorthInterface1,
            SouthInterface1,
            NorthInterface2,
            SouthInterface2,
            NorthInterface3,
            SouthInterface3,
            OutsideInterface1,
            OutsideInterface3,
            Interface_1_2,
            Interface_2_3>;
};

struct RevPatch3
{
    using Interface_1_2 = Interface<EastEdge<1>, WestEdge<2>, true>;
    using Interface_2_3 = Interface<EastEdge<3>, EastEdge<2>, false>;

    using OutsideInterface1 = Interface<OutsideEdge, WestEdge<1>, true>;
    using OutsideInterface3 = Interface<WestEdge<3>, OutsideEdge, true>;

    using CoordTransform1 = AlignedPatchTransform<1>;
    using CoordTransform2 = AlignedPatchTransform<2>;
    using CoordTransform3 = ReversePatchTransform<3>;

    using Connectivity = MultipatchConnectivity<
            NorthInterface1,
            SouthInterface1,
            NorthInterface2,
            SouthInterface2,
            NorthInterface3,
            SouthInterface3,
            OutsideInterface1,
            OutsideInterface3,
            Interface_1_2,
            Interface_2_3>;
};

struct ChangeBound1
{
    using Interface_1_2 = Interface<SouthEdge<1>, WestEdge<2>, true>;
    using Interface_2_3 = Interface<EastEdge<2>, WestEdge<3>, true>;

    using OutsideInterface1 = Interface<OutsideEdge, NorthEdge<1>, true>;
    using OutsideInterface3 = Interface<EastEdge<3>, OutsideEdge, true>;

    using CoordTransform1 = BoundChangeXPatchTransform<1>;
    using CoordTransform2 = AlignedPatchTransform<2>;
    using CoordTransform3 = AlignedPatchTransform<3>;

    using Connectivity = MultipatchConnectivity<
            EastInterface1,
            WestInterface1,
            NorthInterface2,
            SouthInterface2,
            NorthInterface3,
            SouthInterface3,
            OutsideInterface1,
            OutsideInterface3,
            Interface_1_2,
            Interface_2_3>;
};

struct ChangeBound3
{
    using Interface_1_2 = Interface<EastEdge<1>, WestEdge<2>, true>;
    using Interface_2_3 = Interface<EastEdge<2>, SouthEdge<3>, false>;

    using OutsideInterface1 = Interface<OutsideEdge, WestEdge<1>, true>;
    using OutsideInterface3 = Interface<NorthEdge<3>, OutsideEdge, true>;

    using CoordTransform1 = AlignedPatchTransform<1>;
    using CoordTransform2 = AlignedPatchTransform<2>;
    using CoordTransform3 = BoundChangeYPatchTransform<3>;

    using Connectivity = MultipatchConnectivity<
            NorthInterface1,
            SouthInterface1,
            NorthInterface2,
            SouthInterface2,
            EastInterface3,
            WestInterface3,
            OutsideInterface1,
            OutsideInterface3,
            Interface_1_2,
            Interface_2_3>;
};


// CONNECTIVITY ----------------------------------------------------------------------------------

template <class PatchLayout_>
struct InterfaceDerivativeMatrixHermiteFixture : public ::testing::Test
{
    static constexpr int ncells_per_patch = 5;

    // global ------------------------------------
    static constexpr Coord<Xg> xg_min = Coord<Xg> {double(0.0)};
    static constexpr Coord<Xg> xg_max = Coord<Xg> {double(3.0)};

    static constexpr Coord<Yg> yg_min = Coord<Yg> {double(0.0)};
    static constexpr Coord<Yg> yg_max = Coord<Yg> {double(1.0)};

    const IdxRange<GridX<1>> idx_range_x1;
    const IdxRange<GridX<2>> idx_range_x2;
    const IdxRange<GridX<3>> idx_range_x3;

    const IdxRange<GridY<1>> idx_range_y1;
    const IdxRange<GridY<2>> idx_range_y2;
    const IdxRange<GridY<3>> idx_range_y3;

    const IdxRange<GridX<1>, GridY<1>> idx_range_xy1;
    const IdxRange<GridX<2>, GridY<2>> idx_range_xy2;
    const IdxRange<GridX<3>, GridY<3>> idx_range_xy3;

    const IdxRange<GridXg> idx_range_xg;
    const IdxRange<GridYg> idx_range_yg;
    const IdxRange<GridXg, GridYg> idx_range_xy_g;

    using PatchLayout = PatchLayout_;

    typename PatchLayout::CoordTransform1 coord_transform_1;
    typename PatchLayout::CoordTransform2 coord_transform_2;
    typename PatchLayout::CoordTransform3 coord_transform_3;

public:
    InterfaceDerivativeMatrixHermiteFixture()
        : idx_range_x1(SplineInterpPointsX<1>::template get_domain<GridX<1>>())
        , idx_range_x2(SplineInterpPointsX<2>::template get_domain<GridX<2>>())
        , idx_range_x3(SplineInterpPointsX<3>::template get_domain<GridX<3>>())
        , idx_range_y1(SplineInterpPointsY<1>::template get_domain<GridY<1>>())
        , idx_range_y2(SplineInterpPointsY<2>::template get_domain<GridY<2>>())
        , idx_range_y3(SplineInterpPointsY<3>::template get_domain<GridY<3>>())
        , idx_range_xy1(idx_range_x1, idx_range_y1)
        , idx_range_xy2(idx_range_x2, idx_range_y2)
        , idx_range_xy3(idx_range_x3, idx_range_y3)
        , idx_range_xg(SplineInterpPointsXg::get_domain<GridXg>())
        , idx_range_yg(SplineInterpPointsYg::get_domain<GridYg>())
        , idx_range_xy_g(idx_range_xg, idx_range_yg)
        , coord_transform_1(Coord<Xg, Yg>(xg_min, yg_min))
        , coord_transform_2(Coord<Xg, Yg>(xg_min + (xg_max - xg_min) / 3, yg_min))
        , coord_transform_3(Coord<Xg, Yg>(xg_min + 2 * (xg_max - xg_min) / 3, yg_min))
    {
    }

    // INITIALISE DOMAINS ------------------------------------------------------------------------
    static void SetUpTestSuite()
    {
        init_global_mesh();
        IdxRange<GridXg> idx_range_xg(SplineInterpPointsXg::get_domain<GridXg>());
        IdxRange<GridYg> idx_range_yg(SplineInterpPointsYg::get_domain<GridYg>());
        IdxStep<GridXg> nbreaks_per_patch(ncells_per_patch + 1);

        init_patch(
                typename PatchLayout::CoordTransform1(Coord<Xg, Yg>(xg_min, yg_min)),
                idx_range_xg.take_first(nbreaks_per_patch),
                idx_range_yg);
        init_patch(
                typename PatchLayout::CoordTransform2(Coord<Xg, Yg>(
                        ddc::coordinate(idx_range_xg.front() + IdxStep<GridXg>(ncells_per_patch)),
                        yg_min)),
                idx_range_xg.remove(nbreaks_per_patch - 1, nbreaks_per_patch - 1),
                idx_range_yg);
        init_patch(
                typename PatchLayout::CoordTransform3(Coord<Xg, Yg>(
                        ddc::coordinate(
                                idx_range_xg.front() + IdxStep<GridXg>(2 * ncells_per_patch)),
                        yg_min)),
                idx_range_xg.take_last(nbreaks_per_patch),
                idx_range_yg);
    }

    static void init_global_mesh()
    {
        std::vector<Coord<Xg>> break_points_xg = build_random_non_uniform_break_points(
                xg_min,
                xg_max,
                IdxStep<GridXg>(ncells_per_patch * 3));
        std::vector<Coord<Yg>> break_points_yg = build_random_non_uniform_break_points(
                yg_min,
                yg_max,
                IdxStep<GridYg>(ncells_per_patch));

        ddc::init_discrete_space<BSplinesXg>(break_points_xg);
        ddc::init_discrete_space<BSplinesYg>(break_points_yg);

        ddc::init_discrete_space<GridXg>(break_points_xg);
        ddc::init_discrete_space<GridYg>(break_points_yg);
    }

    template <class PatchTransform>
    static void init_patch(
            PatchTransform const& patch_transform,
            IdxRange<GridXg> idx_range_xg,
            IdxRange<GridYg> idx_range_yg)
    {
        using XLinearTransform = typename PatchTransform::XTransform;
        using YLinearTransform = typename PatchTransform::YTransform;
        std::vector<typename XLinearTransform::CoordResult> break_points_along_xg;
        std::vector<typename YLinearTransform::CoordResult> break_points_along_yg;
        ddc::host_for_each(idx_range_xg, [&](Idx<GridXg> idx) {
            break_points_along_xg.push_back(patch_transform.x_transform(ddc::coordinate(idx)));
        });
        ddc::host_for_each(idx_range_yg, [&](Idx<GridYg> idx) {
            break_points_along_yg.push_back(patch_transform.y_transform(ddc::coordinate(idx)));
        });

        using LocalXg = ddc::
                type_seq_element_t<0, ddc::to_type_seq_t<typename XLinearTransform::CoordResult>>;
        using LocalYg = ddc::
                type_seq_element_t<0, ddc::to_type_seq_t<typename YLinearTransform::CoordResult>>;

        using GridLocAlongXg = find_grid_t<LocalXg, AllGrids>;
        using GridLocAlongYg = find_grid_t<LocalYg, AllGrids>;
        using BSplAlongXg = find_grid_t<LocalXg, AllBSpls>;
        using BSplAlongYg = find_grid_t<LocalYg, AllBSpls>;


        if (break_points_along_xg.back() > break_points_along_xg.front()) {
            ddc::init_discrete_space<BSplAlongXg>(break_points_along_xg);
            ddc::init_discrete_space<GridLocAlongXg>(break_points_along_xg);
        } else {
            ddc::init_discrete_space<
                    BSplAlongXg>(break_points_along_xg.rbegin(), break_points_along_xg.rend());
            ddc::init_discrete_space<
                    GridLocAlongXg>(break_points_along_xg.rbegin(), break_points_along_xg.rend());
        }

        if (break_points_along_yg.back() > break_points_along_yg.front()) {
            ddc::init_discrete_space<BSplAlongYg>(break_points_along_yg);
            ddc::init_discrete_space<GridLocAlongYg>(break_points_along_yg);
        } else {
            ddc::init_discrete_space<
                    BSplAlongYg>(break_points_along_yg.rbegin(), break_points_along_yg.rend());
            ddc::init_discrete_space<
                    GridLocAlongYg>(break_points_along_yg.rbegin(), break_points_along_yg.rend());
        }
    }
};

using TestTypes = ::testing::Types<RevPatch1, RevPatch2, RevPatch3, ChangeBound1, ChangeBound3>;
TYPED_TEST_SUITE(InterfaceDerivativeMatrixHermiteFixture, TestTypes);
} // end namespace



TYPED_TEST(InterfaceDerivativeMatrixHermiteFixture, CheckForHermiteBc)
{
    using PatchLayout = typename TestFixture::PatchLayout;
    using Interface_1_2 = typename PatchLayout::Interface_1_2;
    using Interface_2_3 = typename PatchLayout::Interface_2_3;
    std::tuple coord_transforms(
            this->coord_transform_1,
            this->coord_transform_2,
            this->coord_transform_3);

    // Instantiate the derivatives calculators ---------------------------------------------------
    // SingleInterfaceDerivativesCalculators for interfaces along y (periodic).
    SingleInterfaceDerivativesCalculator<Interface_1_2> const
            derivatives_calculator_1_2(this->idx_range_xy1, this->idx_range_xy2);
    SingleInterfaceDerivativesCalculator<Interface_2_3> const
            derivatives_calculator_2_3(this->idx_range_xy2, this->idx_range_xy3);

    // Collect the derivative calculators --------------------------------------------------------
    // We do not follow the physical order to test the operator.
    SingleInterfaceDerivativesCalculatorCollection
            deriv_calculators_collect(derivatives_calculator_2_3, derivatives_calculator_1_2);

    // Collect the index ranges ------------------------------------------------------------------
    MultipatchType<IdxRangeOnPatch, Patch1, Patch2, Patch3>
            idx_ranges(this->idx_range_xy1, this->idx_range_xy2, this->idx_range_xy3);

    using Connectivity = typename PatchLayout::Connectivity;

    // Instantiate the matrix calculators --------------------------------------------------------
    using Grid1AlongXg = find_grid_t<
            ddc::type_seq_element_t<
                    0,
                    ddc::to_type_seq_t<
                            typename PatchLayout::CoordTransform1::XTransform::CoordResult>>,
            AllGrids>;

    InterfaceDerivativeMatrix<
            Connectivity,
            Grid1AlongXg,
            ddc::detail::TypeSeq<Patch1, Patch2, Patch3>,
            SingleInterfaceDerivativesCalculatorCollection<Interface_2_3, Interface_1_2>>
            matrix(idx_ranges, deriv_calculators_collect);

    // Instantiate DerivField ====================================================================
    // Instantiate index range slices ------------------------------------------------------------
    IdxRangeSlice<GridX<1>> idx_range_slice_dx1 = get_bound_idx_range_slice(this->idx_range_x1);
    IdxRangeSlice<GridX<2>> idx_range_slice_dx2 = get_bound_idx_range_slice(this->idx_range_x2);
    IdxRangeSlice<GridX<3>> idx_range_slice_dx3 = get_bound_idx_range_slice(this->idx_range_x3);

    IdxRangeSlice<GridY<1>> idx_range_slice_dy1 = get_bound_idx_range_slice(this->idx_range_y1);
    IdxRangeSlice<GridY<2>> idx_range_slice_dy2 = get_bound_idx_range_slice(this->idx_range_y2);
    IdxRangeSlice<GridY<3>> idx_range_slice_dy3 = get_bound_idx_range_slice(this->idx_range_y3);

    // Instantiate DerivField --------------------------------------------------------------------
    DerivFieldMemOnPatch_host<Patch1> function_and_derivs_1_alloc(
            this->idx_range_xy1,
            idx_range_slice_dx1,
            idx_range_slice_dy1);
    DerivFieldMemOnPatch_host<Patch2> function_and_derivs_2_alloc(
            this->idx_range_xy2,
            idx_range_slice_dx2,
            idx_range_slice_dy2);
    DerivFieldMemOnPatch_host<Patch3> function_and_derivs_3_alloc(
            this->idx_range_xy3,
            idx_range_slice_dx3,
            idx_range_slice_dy3);

    DerivFieldOnPatch_host<Patch1> function_and_derivs_1(function_and_derivs_1_alloc);
    DerivFieldOnPatch_host<Patch2> function_and_derivs_2(function_and_derivs_2_alloc);
    DerivFieldOnPatch_host<Patch3> function_and_derivs_3(function_and_derivs_3_alloc);

    // Collect the fields with derivatives.
    MultipatchField<DerivFieldOnPatch_host, Patch1, Patch2, Patch3> functions_and_derivs(
            function_and_derivs_1,
            function_and_derivs_2,
            function_and_derivs_3);

    // Instantiate the global function.
    IdxRangeSlice<GridXg> idx_range_slice_dxg = get_bound_idx_range_slice(this->idx_range_xg);
    IdxRangeSlice<GridYg> idx_range_slice_dyg = get_bound_idx_range_slice(this->idx_range_yg);

    DerivFieldMem<double, IdxRange<DerivXg, GridXg, DerivYg, GridYg>, 1>
            function_and_derivs_g_alloc(
                    this->idx_range_xy_g,
                    idx_range_slice_dxg,
                    idx_range_slice_dyg);
    DerivField<double, IdxRange<DerivXg, GridXg, DerivYg, GridYg>> function_and_derivs_g(
            function_and_derivs_g_alloc);

    // Initialise the data =======================================================================
    // --- the function values.
    initialise_all_functions<Xg, Yg>(functions_and_derivs, coord_transforms);
    initialise_2D_function(function_and_derivs_g.get_values_field());

    // --- the derivatives of the equivalent global spline.
    Idx<DerivXg> first_dxg(1);
    Idx<DerivYg> first_dyg(1);

    Idx<GridXg> idx_xg_min(this->idx_range_xg.front());
    Idx<GridXg> idx_xg_max(this->idx_range_xg.back());
    Idx<GridYg> idx_yg_min(this->idx_range_yg.front());
    Idx<GridYg> idx_yg_max(this->idx_range_yg.back());

    Idx<DerivXg, GridXg> idx_dxg_min(first_dxg, idx_xg_min);
    Idx<DerivXg, GridXg> idx_dxg_max(first_dxg, idx_xg_max);
    Idx<DerivYg, GridYg> idx_dyg_min(first_dyg, idx_yg_min);
    Idx<DerivYg, GridYg> idx_dyg_max(first_dyg, idx_yg_max);

    Idx<DerivXg, GridXg, DerivYg, GridYg> idx_dxgdyg_min_min(idx_dxg_min, idx_dyg_min);
    Idx<DerivXg, GridXg, DerivYg, GridYg> idx_dxgdyg_max_min(idx_dxg_max, idx_dyg_min);
    Idx<DerivXg, GridXg, DerivYg, GridYg> idx_dxgdyg_min_max(idx_dxg_min, idx_dyg_max);
    Idx<DerivXg, GridXg, DerivYg, GridYg> idx_dxgdyg_max_max(idx_dxg_max, idx_dyg_max);

    double const coef_a = 2. / 3 * M_PI;
    double const coef_b = 0.25;

    double const xg_min = this->xg_min;
    double const xg_max = this->xg_max;
    double const yg_min = this->yg_min;
    double const yg_max = this->yg_max;
    ddc::host_for_each(this->idx_range_yg, [&](Idx<GridYg> const idx) {
        double const yg = ddc::coordinate(idx);
        function_and_derivs_g[idx_dxg_min](idx)
                = -coef_a * std::sin(coef_a * xg_min + coef_b) * std::sin(yg);
        function_and_derivs_g[idx_dxg_max](idx)
                = -coef_a * std::sin(coef_a * xg_max + coef_b) * std::sin(yg);
    });
    ddc::host_for_each(this->idx_range_xg, [&](Idx<GridXg> const idx) {
        double const xg = ddc::coordinate(Idx<GridXg>(idx));
        function_and_derivs_g[idx_dyg_min](idx)
                = std::cos(coef_a * xg + coef_b) * std ::cos(yg_min);
        function_and_derivs_g[idx_dyg_max](idx)
                = std::cos(coef_a * xg + coef_b) * std ::cos(yg_max);
    });
    function_and_derivs_g(idx_dxgdyg_min_min)
            = -coef_a * std::sin(coef_a * xg_min + coef_b) * std::sin(yg_min);
    function_and_derivs_g(idx_dxgdyg_max_min)
            = -coef_a * std::sin(coef_a * xg_max + coef_b) * std::sin(yg_min);
    function_and_derivs_g(idx_dxgdyg_min_max)
            = -coef_a * std::sin(coef_a * xg_min + coef_b) * std::sin(yg_max);
    function_and_derivs_g(idx_dxgdyg_max_max)
            = -coef_a * std::sin(coef_a * xg_max + coef_b) * std::sin(yg_max);

    // --- the local derivatives from an equivalent global spline.
    // ------- build global spline representation
    SplineRThetagBuilder builder_g(this->idx_range_xy_g);
    SplineRThetagBuilderDerivField apply_builder_g(builder_g);

    host_t<DFieldMem<IdxRange<BSplinesXg, BSplinesYg>>> function_g_coef_alloc(
            builder_g.batched_spline_domain(this->idx_range_xy_g));
    host_t<DField<IdxRange<BSplinesXg, BSplinesYg>>> function_g_coef
            = get_field(function_g_coef_alloc);

    apply_builder_g(function_g_coef, function_and_derivs_g);

    host_t<DConstField<IdxRange<BSplinesXg, BSplinesYg>>> const_function_g_coef
            = get_const_field(function_g_coef);

    // ------ global spline evaluator
    ddc::ConstantExtrapolationRule<Yg, Xg> bc_ymin_g(this->yg_min, this->xg_min, this->xg_max);
    ddc::ConstantExtrapolationRule<Yg, Xg> bc_ymax_g(this->yg_max, this->xg_min, this->xg_max);
    ddc::ConstantExtrapolationRule<Xg, Yg> bc_xmin_g(this->xg_min, this->yg_min, this->yg_max);
    ddc::ConstantExtrapolationRule<Xg, Yg> bc_xmax_g(this->xg_max, this->yg_min, this->yg_max);
    SplineRThetagEvaluator evaluator_g(bc_xmin_g, bc_xmax_g, bc_ymin_g, bc_ymax_g);

    IdxRange<GridXg> idx_range_x_global_patch_1(
            this->idx_range_xg.take_first(IdxStep<GridXg>(this->ncells_per_patch + 1)));
    IdxRange<GridXg> idx_range_y_global_patch_2(
            this->idx_range_xg
                    .remove(IdxStep<GridXg>(this->ncells_per_patch),
                            IdxStep<GridXg>(this->ncells_per_patch)));
    IdxRange<GridXg> idx_range_y_global_patch_3(
            this->idx_range_xg.take_last(IdxStep<GridXg>(this->ncells_per_patch + 1)));

    // ------ initialise the boundary first derivatives from the global spline
    // Left X bound ---
    initialise_derivatives<Patch1>(
            this->idx_range_xg.front(),
            function_and_derivs_1,
            evaluator_g,
            const_function_g_coef,
            this->coord_transform_1);

    // lower Y bound ---
    initialise_derivatives<Patch1>(
            this->idx_range_yg.front(),
            function_and_derivs_1,
            evaluator_g,
            const_function_g_coef,
            this->coord_transform_1);
    initialise_derivatives<Patch2>(
            this->idx_range_yg.front(),
            function_and_derivs_2,
            evaluator_g,
            const_function_g_coef,
            this->coord_transform_2);
    initialise_derivatives<Patch3>(
            this->idx_range_yg.front(),
            function_and_derivs_3,
            evaluator_g,
            const_function_g_coef,
            this->coord_transform_3);

    // Right X bound ---
    initialise_derivatives<Patch3>(
            this->idx_range_xg.back(),
            function_and_derivs_3,
            evaluator_g,
            const_function_g_coef,
            this->coord_transform_3);

    // upper Y bound ---
    initialise_derivatives<Patch1>(
            this->idx_range_yg.back(),
            function_and_derivs_1,
            evaluator_g,
            const_function_g_coef,
            this->coord_transform_1);
    initialise_derivatives<Patch2>(
            this->idx_range_yg.back(),
            function_and_derivs_2,
            evaluator_g,
            const_function_g_coef,
            this->coord_transform_2);
    initialise_derivatives<Patch3>(
            this->idx_range_yg.back(),
            function_and_derivs_3,
            evaluator_g,
            const_function_g_coef,
            this->coord_transform_3);

    // ------ initialise the cross-derivatives from the global spline
    initialise_all_cross_derivatives(
            functions_and_derivs,
            evaluator_g,
            const_function_g_coef,
            coord_transforms);

    // --- the first derivatives (on inner interfaces) from the function values.
    matrix.solve_deriv(functions_and_derivs);

    // --- the cross-derivatives from the first derivatives.
    /*
        Here, it is not needed to compute the cross-derivatives because
        they are given by the boundary conditions. However, we want to 
        check that the matrix computes correctly the values. 
    */
    matrix.solve_cross_deriv(functions_and_derivs);

    // Test the values of the derivatives ========================================================
    using EmptyPatchSeq = ddc::detail::TypeSeq<>;

    // Check each derivatives ---
    check_all_x_derivatives(
            functions_and_derivs,
            evaluator_g,
            const_function_g_coef,
            coord_transforms,
            5e-2);
    check_all_y_derivatives<EmptyPatchSeq, EmptyPatchSeq>(
            functions_and_derivs,
            evaluator_g,
            const_function_g_coef,
            coord_transforms,
            5e-2);
    check_all_xy_derivatives<EmptyPatchSeq, EmptyPatchSeq>(
            functions_and_derivs,
            evaluator_g,
            const_function_g_coef,
            coord_transforms,
            5e-2);

    // Check the whole spline representations ---
    check_all_spline_representation_agreement<EmptyPatchSeq, EmptyPatchSeq>(
            functions_and_derivs,
            evaluator_g,
            const_function_g_coef,
            coord_transforms,
            1e-6);
}
