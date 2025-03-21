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

/**
 * @brief A class that represents a polynomial test function for computing partial derivatives.
 * The polynomial depends on two variables.
 */
template <class DerivativeDimension, std::size_t spline_degree>
class FunctionToDifferentiatePolynomial
{
    static_assert(std::is_same_v<DerivativeDimension, X> || std::is_same_v<DerivativeDimension, Y>);

    using DDim = DerivativeDimension;
    using ODim = std::conditional_t<std::is_same_v<DDim, X>, Y, X>;

    using CoordFull = Coord<DDim, ODim>;

    int const m_degree;

public:
    FunctionToDifferentiatePolynomial() : m_degree(spline_degree + 1) {}
    /**
     * @brief Get the value of the function at given coordinate.
     *
     * @param[in] coord_xy The coordinate where we want to evaluate
     * the function.
     *
     * @return The value of the function at the coordinate.
     */
    KOKKOS_FUNCTION double operator()(CoordFull const coord) const
    {
        double const x = ddc::get<X>(coord);
        double const y = ddc::get<Y>(coord);
        return (ipow(x, m_degree + 1) + ipow(y, m_degree)) * y;
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
    KOKKOS_FUNCTION double differentiate(CoordFull const coord) const
    {
        double const x = ddc::get<X>(coord);
        double const y = ddc::get<Y>(coord);

        if constexpr (std::is_same_v<DerivativeDimension, X>) {
            return (1. + m_degree) * ipow(x, m_degree) * y;
        } else {
            return ipow(x, m_degree + 1) + (m_degree + 1) * ipow(y, m_degree);
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
                            - function_to_differentiate.differentiate(ddc::coordinate(idx)));
                });

        return max_error;
    }
};

#if 0
/**
 * @brief A class that represents a test for partial derivatives.
 * The test can be used with 1d splines for computing partial 
 * derivatives.
 */
template <class DerivativeDimension, std::size_t N_ddim, std::size_t N_odim>
class PartialDerivativeTestSpline1D
    : public PartialDerivativeTest<DerivativeDimension, N_ddim, N_odim>
{
private:
    using base_type = PartialDerivativeTest<DerivativeDimension, N_ddim, N_odim>;
    using typename base_type::DDim;

    using typename base_type::IdxRangeDDim;
    using typename base_type::IdxRangeXY;
    using typename base_type::IdxRangeODim;

    using typename base_type::CoordDDim;
    using typename base_type::CoordODim;

    using typename base_type::GridDDim;
    using typename base_type::GridODim;

    using typename base_type::IdxStepDDim;
    using typename base_type::IdxStepODim;

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
            GridDDim,
            GridODim>;

    using SplineDDimEvaluator = ddc::SplineEvaluator<
            Kokkos::DefaultExecutionSpace,
            Kokkos::DefaultExecutionSpace::memory_space,
            BSplinesDDim,
            GridDDim,
            ddc::ConstantExtrapolationRule<DDim>,
            ddc::ConstantExtrapolationRule<DDim>,
            GridDDim,
            GridODim>;

    ddc::ConstantExtrapolationRule<DDim> const m_bv_min;
    ddc::ConstantExtrapolationRule<DDim> const m_bv_max;

public:
    PartialDerivativeTestSpline1D(
            double const ddim_min,
            double const ddim_max,
            double const odim_min,
            double const odim_max)
        : base_type(ddim_min, ddim_max, odim_min, odim_max)
        , m_bv_min(CoordDDim(ddim_min))
        , m_bv_max(CoordDDim(ddim_max))
    {
        std::vector<CoordDDim> point_sampling_ddim = build_random_non_uniform_break_points(
                base_type::m_ddim_min,
                base_type::m_ddim_max,
                base_type::m_ncells_ddim);

        ddc::init_discrete_space<BSplinesDDim>(point_sampling_ddim);
        ddc::init_discrete_space<GridDDim>(
                SplineInterpPointsDDim::template get_sampling<GridDDim>());


        std::vector<CoordODim> point_sampling_odim = build_random_non_uniform_break_points(
                base_type::m_odim_min,
                base_type::m_odim_max,
                base_type::m_ncells_odim);

        ddc::init_discrete_space<BSplinesODim>(point_sampling_odim);
        ddc::init_discrete_space<GridODim>(
                SplineInterpPointsODim::template get_sampling<GridODim>());
    }

    double const operator()(double& delta_ddim) const
    {
        IdxRangeDDim const idxrange_ddim = SplineInterpPointsDDim::template get_domain<GridDDim>();
        IdxRangeODim const idxrange_odim = SplineInterpPointsODim::template get_domain<GridODim>();
        IdxRangeXY const idxrange = IdxRangeXY(idxrange_ddim, idxrange_odim);

        SplineDDimBuilder const spline_builder(idxrange);
        SplineDDimEvaluator const spline_evaluator(m_bv_min, m_bv_max);

        Spline1DPartialDerivativeCreator<SplineDDimBuilder, SplineDDimEvaluator> const
                derivative_creator(spline_builder, spline_evaluator);

        using FunDiff = FunctionToDifferentiatePolynomial<typename base_type::DDim, spline_degree>;
        FunDiff function_to_differentiate;
        double const max_error = base_type::template compute_max_error<
                FunDiff>(idxrange, function_to_differentiate, derivative_creator, delta_ddim);

        return max_error;
    }
};
#endif

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
                base_type::m_ncells_x);
        ddc::init_discrete_space<typename base_type::GridX>(point_sampling_x);

        std::vector<CoordY> point_sampling_y = build_random_non_uniform_break_points(
                base_type::m_ymin,
                base_type::m_ymax,
                base_type::m_ncells_y);
        ddc::init_discrete_space<typename base_type::GridY>(point_sampling_y);
    }

    template <class DerivativeDimension>
    double const compute_error(double& max_distance) const
    {
        typename base_type::IdxRangeX
                idxrange_x(typename base_type::IdxX(0), base_type::m_ncells_x + 1);
        typename base_type::IdxRangeY
                idxrange_y(typename base_type::IdxY(0), base_type::m_ncells_y + 1);
        typename base_type::IdxRangeXY idxrange_xy(idxrange_x, idxrange_y);

        CentralFDMPartialDerivativeCreator<
                typename base_type::IdxRangeXY,
                DerivativeDimension> const derivative_creator;

        using FunDiff = FunctionToDifferentiatePolynomial<DerivativeDimension, 4>;
        FunDiff function_to_differentiate;
        double const max_error = base_type::template compute_max_error<
                DerivativeDimension,
                FunDiff>(idxrange_xy, function_to_differentiate, derivative_creator, max_distance);

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
#if 0
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
    PartialDerivativeTestSpline1D<X, 10, 10> const test_low_x(xmin, xmax, ymin, ymax);
    PartialDerivativeTestSpline1D<X, 100, 10> const test_high_x(xmin, xmax, ymin, ymax);
    double const error_low_x = test_low_x(delta_low_x);
    double const error_high_x = test_high_x(delta_high_x);


    double const order_x
            = std::log(error_high_x / error_low_x) / std::log(delta_high_x / delta_low_x);
    double const relative_error_order_x = std::fabs((spline_degree - order_x) / spline_degree);

    EXPECT_LE(relative_error_order_x, TOL);

    // Partial Derivative in Y direction
    double delta_low_y,
            delta_high_y; // the maximum distance between points in the derivative direction
    PartialDerivativeTestSpline1D<Y, 10, 10> const test_low_y(ymin, ymax, xmin, xmax);
    PartialDerivativeTestSpline1D<Y, 100, 10> const test_high_y(ymin, ymax, xmin, xmax);
    double const error_low_y = test_low_y(delta_low_y);
    double const error_high_y = test_high_y(delta_high_y);

    double const order_y
            = std::log(error_high_y / error_low_y) / std::log(delta_high_y / delta_low_y);
    double const relative_error_order_y = std::fabs((spline_degree - order_y) / spline_degree);

    EXPECT_LE(relative_error_order_y, TOL);
}
#endif

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

#if 0
    // Partial Derivative in Y direction
    double delta_low_y, delta_high_y; 
    double const error_low_y = test_low_res.compute_error<Y>(delta_low_y);
    double const error_high_y = test_low_res.compute_error<Y>(delta_high_y);

    double const order_y
            = std::log(error_high_y / error_low_y) / std::log(delta_high_y / delta_low_y);
    double const relative_error_order_y = std::fabs((expected_order - order_y) / expected_order);

    EXPECT_LE(relative_error_order_y, TOL);
#endif
}
} // namespace
