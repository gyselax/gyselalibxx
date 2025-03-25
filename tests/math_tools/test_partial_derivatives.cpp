// SPDX-License-Identifier: MIT
#include <ddc/ddc.hpp>
#include <ddc/kernels/splines.hpp>

#include <gtest/gtest.h>

#include "central_fdm_partial_derivatives.hpp"
#include "ddc_aliases.hpp"
#include "ddc_helper.hpp"
#include "math_tools.hpp"
#include "mesh_builder.hpp"
#include "spline_1d_partial_derivative.hpp"
#include "spline_2d_partial_derivative.hpp"

namespace {

struct X
{
    static bool constexpr PERIODIC = false;
};

struct Y
{
    static bool constexpr PERIODIC = false;
};

auto static constexpr SplineBoundary = ddc::BoundCond::GREVILLE;

std::size_t static constexpr spline_degree = 3;

using CoordX = Coord<X>;
using CoordY = Coord<Y>;
using CoordXY = Coord<X, Y>;


/**
 * @brief A class that represents a test function for computing partial derivatives.
 * The test function is defined as the product of two cosine functions.
 */
class FunctionToDifferentiateCosine
{
public:
    /**
     * @brief Default constructor.
     */
    KOKKOS_DEFAULTED_FUNCTION FunctionToDifferentiateCosine() = default;

    /**
     * @brief Default copy constructor.
     */
    KOKKOS_DEFAULTED_FUNCTION FunctionToDifferentiateCosine(FunctionToDifferentiateCosine const&)
            = default;

    /**
     * @brief Default destructor.
     */
    KOKKOS_DEFAULTED_FUNCTION ~FunctionToDifferentiateCosine() = default;

    /**
     * @brief Get the value of the function at given coordinate.
     *
     * @param[in] coord_xy The coordinate where we want to evaluate
     * the function.
     *
     * @return The value of the function at the coordinate.
     */
    KOKKOS_FUNCTION double operator()(CoordXY const coord_xy) const
    {
        double const x = ddc::get<X>(coord_xy);
        double const y = ddc::get<Y>(coord_xy);
        return Kokkos::cos(x) * Kokkos::cos(y);
    }

    /**
     * @brief Get the value of the partial derivative of the function
     * in the DerivativeDimension direction at a given coordinate.
     *
     * @param[in] coord_xy The coordinate where we want to evaluate 
     * the partial derivative.
     *
     * @return The value of the partial derivative of the function
     * at the coordinate.
     */
    template <class DerivativeDimension>
    KOKKOS_FUNCTION double differentiate(CoordXY const coord_xy) const
    {
        static_assert(
                std::is_same_v<DerivativeDimension, X> || std::is_same_v<DerivativeDimension, Y>);

        double const x = ddc::get<X>(coord_xy);
        double const y = ddc::get<Y>(coord_xy);

        if constexpr (std::is_same_v<DerivativeDimension, X>) {
            return -Kokkos::sin(x) * Kokkos::cos(y);
        } else {
            return -Kokkos::cos(x) * Kokkos::sin(y);
        }
    }
};

/**
 * @brief A class that represents a test for partial derivatives.
 * The test can be used with several implementations for computing
 * partial derivatives.
 */
template <std::size_t ncells_x, std::size_t ncells_y>
class PartialDerivativeTest
{
public:
    struct GridX : NonUniformGridBase<X>
    {
    };

    struct GridY : NonUniformGridBase<Y>
    {
    };

    using IdxRangeX = IdxRange<GridX>;
    using IdxRangeY = IdxRange<GridY>;
    using IdxRangeXY = IdxRange<GridX, GridY>;

    using IdxX = Idx<GridX>;
    using IdxY = Idx<GridY>;
    using IdxXY = Idx<GridX, GridY>;

    using IdxStepX = IdxStep<GridX>;
    using IdxStepY = IdxStep<GridY>;

    using DFieldMemType = DFieldMem<IdxRangeXY>;
    using DFieldType = DField<IdxRangeXY>;

protected:
    CoordX const m_xmin;
    CoordX const m_xmax;
    CoordY const m_ymin;
    CoordY const m_ymax;

    IdxStepX const m_ncells_x;
    IdxStepY const m_ncells_y;

public:
    PartialDerivativeTest(
            double const xmin,
            double const xmax,
            double const ymin,
            double const ymax)
        : m_xmin(xmin)
        , m_xmax(xmax)
        , m_ymin(ymin)
        , m_ymax(ymax)
        , m_ncells_x(ncells_x)
        , m_ncells_y(ncells_y)
    {
    }

    template <class DerivativeDimension, class FunctionToDifferentiate>
    double compute_max_error(
            IdxRangeXY const& idxrange_xy,
            FunctionToDifferentiate const& function_to_differentiate,
            IPartialDerivativeCreator<IdxRangeXY, DerivativeDimension> const&
                    partial_derivative_creator,
            double& max_distance) const
    {
        using GridDDim = std::conditional_t<std::is_same_v<DerivativeDimension, X>, GridX, GridY>;

        max_distance = ddcHelper::maximum_distance_between_adjacent_points(
                ddc::select<GridDDim>(idxrange_xy));

        // field to be differentiated
        DFieldMemType field_to_differentiate(idxrange_xy);
        DFieldType field_to_differentiate_proxy = get_field(field_to_differentiate);
        ddc::parallel_for_each(
                Kokkos::DefaultExecutionSpace(),
                idxrange_xy,
                KOKKOS_LAMBDA(IdxXY const idx) {
                    field_to_differentiate_proxy(idx)
                            = function_to_differentiate(ddc::coordinate(idx));
                });

        std::unique_ptr<IPartialDerivative<IdxRangeXY, DerivativeDimension>> const
                partial_derivative_creator_pointer
                = partial_derivative_creator.create_instance(
                        get_const_field(field_to_differentiate));

        IPartialDerivative<IdxRangeXY, DerivativeDimension> const& partial_derivative
                = *partial_derivative_creator_pointer;

        DFieldMemType field_differentiated_alloc(idxrange_xy);
        DFieldType field_differentiated = get_field(field_differentiated_alloc);
        partial_derivative(field_differentiated);

        double const max_error = ddc::parallel_transform_reduce(
                Kokkos::DefaultExecutionSpace(),
                idxrange_xy,
                0.,
                ddc::reducer::max<double>(),
                KOKKOS_LAMBDA(IdxXY const idx) {
                    return Kokkos::abs(
                            field_differentiated(idx)
                            - function_to_differentiate.template differentiate<DerivativeDimension>(
                                    ddc::coordinate(idx)));
                });

        return max_error;
    }
};

/**
 * @brief A class that represents a test for partial derivatives.
 * The test can be used with 1d splines for computing partial 
 * derivatives.
 */
template <class DerivativeDimension, std::size_t ncells_x, std::size_t ncells_y>
class PartialDerivativeTestSpline1D : public PartialDerivativeTest<ncells_x, ncells_y>
{
    static_assert(std::is_same_v<DerivativeDimension, X> || std::is_same_v<DerivativeDimension, Y>);

public:
    using base_type = PartialDerivativeTest<ncells_x, ncells_y>;
    using DDim = DerivativeDimension;

    using typename base_type::GridX;
    using typename base_type::GridY;

    using typename base_type::IdxRangeX;
    using typename base_type::IdxRangeXY;
    using typename base_type::IdxRangeY;

    struct BSplinesX : ddc::NonUniformBSplines<X, spline_degree>
    {
    };
    using SplineInterpPointsX
            = ddc::GrevilleInterpolationPoints<BSplinesX, SplineBoundary, SplineBoundary>;


    struct BSplinesY : ddc::NonUniformBSplines<Y, spline_degree>
    {
    };
    using SplineInterpPointsY
            = ddc::GrevilleInterpolationPoints<BSplinesY, SplineBoundary, SplineBoundary>;

    using BSplinesDDim = std::conditional_t<std::is_same_v<DDim, X>, BSplinesX, BSplinesY>;
    using BSplinesODim = std::conditional_t<std::is_same_v<DDim, X>, BSplinesY, BSplinesX>;

    using GridDDim = std::conditional_t<std::is_same_v<DDim, X>, GridX, GridY>;
    using GridODim = std::conditional_t<std::is_same_v<DDim, X>, GridY, GridX>;

    using SplineInterpPointsDDim
            = std::conditional_t<std::is_same_v<DDim, X>, SplineInterpPointsX, SplineInterpPointsY>;
    using SplineInterpPointsODim
            = std::conditional_t<std::is_same_v<DDim, X>, SplineInterpPointsY, SplineInterpPointsX>;

    using SplineDDimBuilder = ddc::SplineBuilder<
            Kokkos::DefaultExecutionSpace,
            Kokkos::DefaultExecutionSpace::memory_space,
            BSplinesDDim,
            GridDDim,
            SplineBoundary,
            SplineBoundary,
            ddc::SplineSolver::LAPACK,
            GridX,
            GridY>;

    using SplineDDimEvaluator = ddc::SplineEvaluator<
            Kokkos::DefaultExecutionSpace,
            Kokkos::DefaultExecutionSpace::memory_space,
            BSplinesDDim,
            GridDDim,
            ddc::ConstantExtrapolationRule<DDim>,
            ddc::ConstantExtrapolationRule<DDim>,
            GridX,
            GridY>;

public:
    PartialDerivativeTestSpline1D(
            double const xmin,
            double const xmax,
            double const ymin,
            double const ymax)
        : base_type(xmin, xmax, ymin, ymax)
    {
        std::vector<CoordX> point_sampling_x = build_random_non_uniform_break_points(
                base_type::m_xmin,
                base_type::m_xmax,
                base_type::m_ncells_x,
                0.2);

        ddc::init_discrete_space<BSplinesX>(point_sampling_x);
        ddc::init_discrete_space<GridX>(SplineInterpPointsX::template get_sampling<GridX>());

        std::vector<CoordY> point_sampling_y = build_random_non_uniform_break_points(
                base_type::m_ymin,
                base_type::m_ymax,
                base_type::m_ncells_y,
                0.2);

        ddc::init_discrete_space<BSplinesY>(point_sampling_y);
        ddc::init_discrete_space<GridY>(SplineInterpPointsY::template get_sampling<GridY>());
    }

    double compute_error(double& max_distance) const
    {
        IdxRangeX const idxrange_x = SplineInterpPointsX::template get_domain<GridX>();
        IdxRangeY const idxrange_y = SplineInterpPointsY::template get_domain<GridY>();
        IdxRangeXY const idxrange = IdxRangeXY(idxrange_x, idxrange_y);

        double dmin, dmax;
        if constexpr (std::is_same_v<DDim, X>) {
            dmin = base_type::m_xmin;
            dmax = base_type::m_xmax;
        } else {
            dmin = base_type::m_ymin;
            dmax = base_type::m_ymax;
        }
        ddc::ConstantExtrapolationRule<DDim> const bv_min(Coord<DDim> {dmin});
        ddc::ConstantExtrapolationRule<DDim> const bv_max(Coord<DDim> {dmax});

        SplineDDimBuilder const spline_builder(idxrange);
        SplineDDimEvaluator const spline_evaluator(bv_min, bv_max);

        Spline1DPartialDerivativeCreator<SplineDDimBuilder, SplineDDimEvaluator> const
                derivative_creator(spline_builder, spline_evaluator);

        FunctionToDifferentiateCosine function_to_differentiate;
        double const max_error = base_type::
                template compute_max_error<DerivativeDimension, FunctionToDifferentiateCosine>(
                        idxrange,
                        function_to_differentiate,
                        derivative_creator,
                        max_distance);

        return max_error;
    }
};

/**
 * @brief A class that represents a test for partial derivatives.
 * The test can be used with 2d splines for computing partial 
 * derivatives.
 */
template <class DerivativeDimension, std::size_t ncells_x, std::size_t ncells_y>
class PartialDerivativeTestSpline2D : public PartialDerivativeTest<ncells_x, ncells_y>
{
    static_assert(std::is_same_v<DerivativeDimension, X> || std::is_same_v<DerivativeDimension, Y>);

public:
    using base_type = PartialDerivativeTest<ncells_x, ncells_y>;
    using DDim = DerivativeDimension;

    using typename base_type::GridX;
    using typename base_type::GridY;

    using typename base_type::IdxRangeX;
    using typename base_type::IdxRangeXY;
    using typename base_type::IdxRangeY;

    struct BSplinesX : ddc::NonUniformBSplines<X, spline_degree>
    {
    };
    using SplineInterpPointsX
            = ddc::GrevilleInterpolationPoints<BSplinesX, SplineBoundary, SplineBoundary>;


    struct BSplinesY : ddc::NonUniformBSplines<Y, spline_degree>
    {
    };
    using SplineInterpPointsY
            = ddc::GrevilleInterpolationPoints<BSplinesY, SplineBoundary, SplineBoundary>;

    using SplineBuilder2D = ddc::SplineBuilder2D<
            Kokkos::DefaultExecutionSpace,
            Kokkos::DefaultExecutionSpace::memory_space,
            BSplinesX,
            BSplinesY,
            GridX,
            GridY,
            SplineBoundary,
            SplineBoundary,
            SplineBoundary,
            SplineBoundary,
            ddc::SplineSolver::LAPACK,
            GridX,
            GridY>;

    using SplineEvaluator2D = ddc::SplineEvaluator2D<
            Kokkos::DefaultExecutionSpace,
            Kokkos::DefaultExecutionSpace::memory_space,
            BSplinesX,
            BSplinesY,
            GridX,
            GridY,
            ddc::ConstantExtrapolationRule<X, Y>,
            ddc::ConstantExtrapolationRule<X, Y>,
            ddc::ConstantExtrapolationRule<Y, X>,
            ddc::ConstantExtrapolationRule<Y, X>,
            GridX,
            GridY>;

    ddc::ConstantExtrapolationRule<X, Y> const m_bv_xmin;
    ddc::ConstantExtrapolationRule<X, Y> const m_bv_xmax;
    ddc::ConstantExtrapolationRule<Y, X> const m_bv_ymin;
    ddc::ConstantExtrapolationRule<Y, X> const m_bv_ymax;

public:
    PartialDerivativeTestSpline2D(
            double const xmin,
            double const xmax,
            double const ymin,
            double const ymax)
        : base_type(xmin, xmax, ymin, ymax)
        , m_bv_xmin(CoordX(xmin), CoordY(ymin), CoordY(ymax))
        , m_bv_xmax(CoordX(xmax), CoordY(ymin), CoordY(ymax))
        , m_bv_ymin(CoordY(ymin), CoordX(xmin), CoordX(xmax))
        , m_bv_ymax(CoordY(ymax), CoordX(xmin), CoordX(xmax))
    {
        std::vector<CoordX> point_sampling_x = build_random_non_uniform_break_points(
                base_type::m_xmin,
                base_type::m_xmax,
                base_type::m_ncells_x,
                0.2);

        ddc::init_discrete_space<BSplinesX>(point_sampling_x);
        ddc::init_discrete_space<GridX>(SplineInterpPointsX::template get_sampling<GridX>());

        std::vector<CoordY> point_sampling_y = build_random_non_uniform_break_points(
                base_type::m_ymin,
                base_type::m_ymax,
                base_type::m_ncells_y,
                0.2);

        ddc::init_discrete_space<BSplinesY>(point_sampling_y);
        ddc::init_discrete_space<GridY>(SplineInterpPointsY::template get_sampling<GridY>());
    }

    double compute_error(double& max_distance) const
    {
        IdxRangeX const idxrange_x = SplineInterpPointsX::template get_domain<GridX>();
        IdxRangeY const idxrange_y = SplineInterpPointsY::template get_domain<GridY>();
        IdxRangeXY const idxrange = IdxRangeXY(idxrange_x, idxrange_y);

        SplineBuilder2D builder(idxrange);
        SplineBuilder2DCache<SplineBuilder2D> builder_cache(builder);
        SplineEvaluator2D evaluator(m_bv_xmin, m_bv_xmax, m_bv_ymin, m_bv_ymax);

        Spline2DPartialDerivativeCreator<
                SplineBuilder2DCache<SplineBuilder2D>,
                SplineEvaluator2D,
                DDim> const derivative_creator(builder_cache, evaluator);

        FunctionToDifferentiateCosine function_to_differentiate;
        double const max_error = base_type::
                template compute_max_error<DerivativeDimension, FunctionToDifferentiateCosine>(
                        idxrange,
                        function_to_differentiate,
                        derivative_creator,
                        max_distance);

        return max_error;
    }
};


/**
 * @brief A class that represents a test for partial derivatives.
 * The test can be used with finite difference method for
 * computing partial derivatives.
 */
template <std::size_t ncells_x, std::size_t ncells_y>
class PartialDerivativeTestFDM : public PartialDerivativeTest<ncells_x, ncells_y>
{
private:
    using base_type = PartialDerivativeTest<ncells_x, ncells_y>;

    using typename base_type::IdxRangeX;
    using typename base_type::IdxRangeXY;
    using typename base_type::IdxRangeY;

public:
    PartialDerivativeTestFDM(
            double const xmin,
            double const xmax,
            double const ymin,
            double const ymax)
        : base_type(xmin, xmax, ymin, ymax)
    {
        std::vector<CoordX> point_sampling_x = build_random_non_uniform_break_points(
                base_type::m_xmin,
                base_type::m_xmax,
                base_type::m_ncells_x,
                0.2);
        ddc::init_discrete_space<typename base_type::GridX>(point_sampling_x);

        std::vector<CoordY> point_sampling_y = build_random_non_uniform_break_points(
                base_type::m_ymin,
                base_type::m_ymax,
                base_type::m_ncells_y,
                0.2);
        ddc::init_discrete_space<typename base_type::GridY>(point_sampling_y);
    }

    template <class DerivativeDimension>
    double compute_error(double& max_distance) const
    {
        IdxRangeX idxrange_x(typename base_type::IdxX(0), base_type::m_ncells_x + 1);
        IdxRangeY idxrange_y(typename base_type::IdxY(0), base_type::m_ncells_y + 1);
        IdxRangeXY idxrange_xy(idxrange_x, idxrange_y);

        CentralFDMPartialDerivativeCreator<IdxRangeXY, DerivativeDimension> const
                derivative_creator;

        FunctionToDifferentiateCosine function_to_differentiate;
        double const max_error = base_type::
                template compute_max_error<DerivativeDimension, FunctionToDifferentiateCosine>(
                        idxrange_xy,
                        function_to_differentiate,
                        derivative_creator,
                        max_distance);

        return max_error;
    }
};


/** 
 * We expect a convergence of the error following error ~ (dx)^d
 * with d the degree of the splines. 
 * We approximate d using two values for dx, dx1 and dx2
 * following @f$ log(error1/error2) / log(dx1/dx2) ~ d @f$. 
 * We check that this relation is satisfied with a given tolerance. 
 * Since we use non uniform interpolation points, dx depends on the 
 * position. We take dx as the maximum distance between two points.
 */
TEST(PartialDerivative, Spline1DPartialDerivative)
{
    double const xmin(0.1);
    double const xmax(1.1);
    double const ymin(0.2);
    double const ymax(1.2);
    // relative error of convergence should be less than 15%
    double const TOL = 0.15;

    // Partial Derivative in X direction
    double delta_low_x,
            delta_high_x; // the maximum distance between points in the derivative direction
    PartialDerivativeTestSpline1D<X, 10, 5> const test_low_x(xmin, xmax, ymin, ymax);
    PartialDerivativeTestSpline1D<X, 100, 5> const test_high_x(xmin, xmax, ymin, ymax);
    double const error_low_x = test_low_x.compute_error(delta_low_x);
    double const error_high_x = test_high_x.compute_error(delta_high_x);


    double const order_x
            = std::log(error_high_x / error_low_x) / std::log(delta_high_x / delta_low_x);
    double const relative_error_order_x = std::fabs((spline_degree - order_x) / spline_degree);

    EXPECT_LE(relative_error_order_x, TOL);

    // Partial Derivative in Y direction
    double delta_low_y,
            delta_high_y; // the maximum distance between points in the derivative direction
    PartialDerivativeTestSpline1D<Y, 5, 10> const test_low_y(xmin, xmax, ymin, ymax);
    PartialDerivativeTestSpline1D<Y, 5, 100> const test_high_y(xmin, xmax, ymin, ymax);
    double const error_low_y = test_low_y.compute_error(delta_low_y);
    double const error_high_y = test_high_y.compute_error(delta_high_y);

    double const order_y
            = std::log(error_high_y / error_low_y) / std::log(delta_high_y / delta_low_y);
    double const relative_error_order_y = std::fabs((spline_degree - order_y) / spline_degree);

    EXPECT_LE(relative_error_order_y, TOL);
}


TEST(PartialDerivative, Spline2DPartialDerivative)
{
    double const xmin(0.1);
    double const xmax(1.1);
    double const ymin(0.2);
    double const ymax(1.2);
    // relative error of convergence should be less than 15%
    double const TOL = 0.15;

    // Partial Derivative in X direction
    double delta_low_x,
            delta_high_x; // the maximum distance between points in the derivative direction
    PartialDerivativeTestSpline2D<X, 10, 5> const test_low_x(xmin, xmax, ymin, ymax);
    PartialDerivativeTestSpline2D<X, 100, 5> const test_high_x(xmin, xmax, ymin, ymax);
    double const error_low_x = test_low_x.compute_error(delta_low_x);
    double const error_high_x = test_high_x.compute_error(delta_high_x);


    double const order_x
            = std::log(error_high_x / error_low_x) / std::log(delta_high_x / delta_low_x);
    double const relative_error_order_x = std::fabs((spline_degree - order_x) / spline_degree);

    EXPECT_LE(relative_error_order_x, TOL);

    // Partial Derivative in Y direction
    double delta_low_y,
            delta_high_y; // the maximum distance between points in the derivative direction
    PartialDerivativeTestSpline2D<Y, 5, 10> const test_low_y(xmin, xmax, ymin, ymax);
    PartialDerivativeTestSpline2D<Y, 5, 100> const test_high_y(xmin, xmax, ymin, ymax);
    double const error_low_y = test_low_y.compute_error(delta_low_y);
    double const error_high_y = test_high_y.compute_error(delta_high_y);

    double const order_y
            = std::log(error_high_y / error_low_y) / std::log(delta_high_y / delta_low_y);
    double const relative_error_order_y = std::fabs((spline_degree - order_y) / spline_degree);

    EXPECT_LE(relative_error_order_y, TOL);
}


TEST(PartialDerivative, CentralFDMPartialDerivative)
{
    double const xmin(0.1);
    double const xmax(1.1);
    double const ymin(0.2);
    double const ymax(1.2);
    // relative error of convergence should be less than 15%
    double const TOL = 0.15;
    int const expected_order(2);

    PartialDerivativeTestFDM<10, 10> const test_low_res(xmin, xmax, ymin, ymax);
    PartialDerivativeTestFDM<100, 100> const test_high_res(xmin, xmax, ymin, ymax);

    // Partial Derivative in X direction
    double delta_low_x, delta_high_x; // max distance between points in the derivative direction
    double const error_low_x = test_low_res.compute_error<X>(delta_low_x);
    double const error_high_x = test_high_res.compute_error<X>(delta_high_x);

    double const order_x
            = std::log(error_high_x / error_low_x) / std::log(delta_high_x / delta_low_x);
    double const relative_error_order_x = std::fabs((expected_order - order_x) / expected_order);

    EXPECT_LE(relative_error_order_x, TOL);

    // Partial Derivative in Y direction
    double delta_low_y, delta_high_y;
    double const error_low_y = test_low_res.compute_error<Y>(delta_low_y);
    double const error_high_y = test_high_res.compute_error<Y>(delta_high_y);

    double const order_y
            = std::log(error_high_y / error_low_y) / std::log(delta_high_y / delta_low_y);
    double const relative_error_order_y = std::fabs((expected_order - order_y) / expected_order);

    EXPECT_LE(relative_error_order_y, TOL);
}
} // namespace
