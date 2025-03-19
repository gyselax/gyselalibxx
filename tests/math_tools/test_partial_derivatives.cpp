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
template <
        class DerivativeDimension,
        std::size_t N_ddim,
        std::size_t N_odim,
        std::size_t spline_degree>
class PartialDerivativeTest
{
    static_assert(std::is_same_v<DerivativeDimension, X> || std::is_same_v<DerivativeDimension, Y>);

public:
    struct BSplinesX : ddc::NonUniformBSplines<X, spline_degree>
    {
    };

    auto static constexpr SplineBoundary = ddc::BoundCond::GREVILLE;

    using SplineInterpPointsX
            = ddc::GrevilleInterpolationPoints<BSplinesX, SplineBoundary, SplineBoundary>;

    struct GridX : NonUniformGridBase<X>
    {
    };

    struct BSplinesY : ddc::NonUniformBSplines<Y, spline_degree>
    {
    };

    using SplineInterpPointsY
            = ddc::GrevilleInterpolationPoints<BSplinesY, SplineBoundary, SplineBoundary>;
    struct GridY : NonUniformGridBase<Y>
    {
    };

    using DDim = DerivativeDimension;
    using ODim = std::conditional_t<std::is_same_v<DDim, X>, Y, X>;

    using CoordDDim = Coord<DDim>;
    using CoordODim = Coord<ODim>;

    using GridDDim = std::conditional_t<std::is_same_v<DDim, X>, GridX, GridY>;
    using GridODim = std::conditional_t<std::is_same_v<DDim, X>, GridY, GridX>;

    using IdxRangeDDim = IdxRange<GridDDim>;
    using IdxRangeODim = IdxRange<GridODim>;
    using IdxRangeFull = IdxRange<GridDDim, GridODim>;

    using IdxFull = Idx<GridDDim, GridODim>;

    using IdxStepDDim = IdxStep<GridDDim>;
    using IdxStepODim = IdxStep<GridODim>;
    using IdxStepFull = IdxStep<GridDDim, GridODim>;

    using BSplinesDDim = std::conditional_t<std::is_same_v<DDim, X>, BSplinesX, BSplinesY>;
    using BSplinesODim = std::conditional_t<std::is_same_v<DDim, X>, BSplinesY, BSplinesX>;

    using SplineInterpPointsDDim
            = std::conditional_t<std::is_same_v<DDim, X>, SplineInterpPointsX, SplineInterpPointsY>;
    using SplineInterpPointsODim
            = std::conditional_t<std::is_same_v<DDim, X>, SplineInterpPointsY, SplineInterpPointsX>;

    using DFieldMemType = DFieldMem<IdxRangeFull>;
    using DFieldType = DField<IdxRangeFull>;

    CoordDDim const m_ddim_min;
    CoordDDim const m_ddim_max;
    CoordODim const m_odim_min;
    CoordODim const m_odim_max;
    IdxStepDDim const m_ncells_ddim;
    IdxStepODim const m_ncells_odim;

    PartialDerivativeTest(
            double const ddim_min,
            double const ddim_max,
            double const odim_min,
            double const odim_max)
        : m_ddim_min(ddim_min)
        , m_ddim_max(ddim_max)
        , m_odim_min(odim_min)
        , m_odim_max(odim_max)
        , m_ncells_ddim(N_ddim)
        , m_ncells_odim(N_odim)
    {
        std::vector<CoordDDim> point_sampling_ddim
                = build_random_non_uniform_break_points(m_ddim_min, m_ddim_max, m_ncells_ddim);

        ddc::init_discrete_space<BSplinesDDim>(point_sampling_ddim);
        ddc::init_discrete_space<GridDDim>(
                SplineInterpPointsDDim::template get_sampling<GridDDim>());


        std::vector<CoordODim> point_sampling_odim
                = build_random_non_uniform_break_points(m_odim_min, m_odim_max, m_ncells_odim);

        ddc::init_discrete_space<BSplinesODim>(point_sampling_odim);
        ddc::init_discrete_space<GridODim>(
                SplineInterpPointsODim::template get_sampling<GridODim>());
    }

    template <class FunctionToDifferentiate>
    double compute_max_error(
            FunctionToDifferentiate const& function_to_differentiate,
            IPartialDerivativeCreator<IdxRangeFull, DerivativeDimension> const&
                    partial_derivative_creator,
            double& delta_ddim) const
    {
        IdxRangeDDim idxrange_ddim(SplineInterpPointsDDim::template get_domain<GridDDim>());
        IdxRangeODim idxrange_odim(SplineInterpPointsODim::template get_domain<GridODim>());
        IdxRangeFull idxrange(idxrange_ddim, idxrange_odim);

        delta_ddim = ddcHelper::maximum_distance_between_adjacent_points(idxrange_ddim);

        // field to be differentiated
        DFieldMemType field_to_differentiate(idxrange);
        DFieldType field_to_differentiate_proxy = get_field(field_to_differentiate);
        ddc::parallel_for_each(
                Kokkos::DefaultExecutionSpace(),
                idxrange,
                KOKKOS_LAMBDA(IdxFull const idx) {
                    field_to_differentiate_proxy(idx)
                            = function_to_differentiate(ddc::coordinate(idx));
                });

        std::unique_ptr<IPartialDerivative<IdxRangeFull, DerivativeDimension>> const
                partial_derivative_creator_pointer
                = partial_derivative_creator.create_instance(
                        get_const_field(field_to_differentiate));

        IPartialDerivative<IdxRangeFull, DerivativeDimension> const& partial_derivative
                = *partial_derivative_creator_pointer;

        DFieldMemType field_differentiated_alloc(idxrange);
        DFieldType field_differentiated = get_field(field_differentiated_alloc);
        partial_derivative(field_differentiated);

        double const max_error = ddc::parallel_transform_reduce(
                Kokkos::DefaultExecutionSpace(),
                idxrange,
                0.,
                ddc::reducer::max<double>(),
                KOKKOS_LAMBDA(IdxFull const idx) {
                    return Kokkos::abs(
                            field_differentiated(idx)
                            - function_to_differentiate.differentiate(ddc::coordinate(idx)));
                });

        return max_error;
    }
};

/**
 * @brief A class that represents a test for partial derivatives.
 * The test can be used with 1d splines for computing partial 
 * derivatives.
 */
template <
        class DerivativeDimension,
        std::size_t N_ddim,
        std::size_t N_odim,
        std::size_t spline_degree>
class PartialDerivativeTestSpline1D
    : public PartialDerivativeTest<DerivativeDimension, N_ddim, N_odim, spline_degree>
{
private:
    using base_type = PartialDerivativeTest<DerivativeDimension, N_ddim, N_odim, spline_degree>;
    using DDim = DerivativeDimension;

    using SplineDDimBuilder = ddc::SplineBuilder<
            Kokkos::DefaultExecutionSpace,
            Kokkos::DefaultExecutionSpace::memory_space,
            typename base_type::BSplinesDDim,
            typename base_type::GridDDim,
            base_type::SplineBoundary,
            base_type::SplineBoundary,
            ddc::SplineSolver::LAPACK,
            typename base_type::GridDDim,
            typename base_type::GridODim>;

    using SplineDDimEvaluator = ddc::SplineEvaluator<
            Kokkos::DefaultExecutionSpace,
            Kokkos::DefaultExecutionSpace::memory_space,
            typename base_type::BSplinesDDim,
            typename base_type::GridDDim,
            ddc::ConstantExtrapolationRule<DDim>,
            ddc::ConstantExtrapolationRule<DDim>,
            typename base_type::GridDDim,
            typename base_type::GridODim>;

    using IdxRangeDDim = typename base_type::IdxRangeDDim;
    using IdxRangeODim = typename base_type::IdxRangeODim;
    using IdxRangeFull = typename base_type::IdxRangeFull;

    using SplineInterpPointsDDim = typename base_type::SplineInterpPointsDDim;
    using SplineInterpPointsODim = typename base_type::SplineInterpPointsODim;

    using GridDDim = typename base_type::GridDDim;
    using GridODim = typename base_type::GridODim;

    ddc::ConstantExtrapolationRule<DDim> const m_bv_min;
    ddc::ConstantExtrapolationRule<DDim> const m_bv_max;

public:
    PartialDerivativeTestSpline1D(
            double const ddim_min,
            double const ddim_max,
            double const odim_min,
            double const odim_max)
        : base_type(ddim_min, ddim_max, odim_min, odim_max)
        , m_bv_min(typename base_type::CoordDDim(ddim_min))
        , m_bv_max(typename base_type::CoordDDim(ddim_max))
    {
    }

    double const operator()(double& delta_ddim) const
    {
        IdxRangeDDim const idxrange_ddim = SplineInterpPointsDDim::template get_domain<GridDDim>();
        IdxRangeODim const idxrange_odim = SplineInterpPointsODim::template get_domain<GridODim>();
        IdxRangeFull const idxrange = IdxRangeFull(idxrange_ddim, idxrange_odim);
        SplineDDimBuilder const builder(idxrange);

        SplineDDimEvaluator const spline_evaluator(m_bv_min, m_bv_max);

        Spline1DPartialDerivativeCreator<SplineDDimBuilder, SplineDDimEvaluator> const
                derivative_creator(builder, spline_evaluator);

        using FunDiff = FunctionToDifferentiatePolynomial<typename base_type::DDim, spline_degree>;
        FunDiff function_to_differentiate;
        double const max_error = base_type::template compute_max_error<
                FunDiff>(function_to_differentiate, derivative_creator, delta_ddim);

        return max_error;
    }
};

/**
 * @brief A class that represents a test for partial derivatives.
<<<<<<< HEAD
 * The test can be used with 2d splines for computing partial 
 * derivatives.
 */
template <
        class DerivativeDimension,
        std::size_t N_ddim,
        std::size_t N_odim,
        std::size_t spline_degree>
class PartialDerivativeTestSpline2D
    : public PartialDerivativeTest<DerivativeDimension, N_ddim, N_odim, spline_degree>
{
private:
    using base_type = PartialDerivativeTest<DerivativeDimension, N_ddim, N_odim, spline_degree>;

    using DDim = DerivativeDimension;
    using ODim = typename base_type::ODim;

    using CoordDDim = typename base_type::CoordDDim;
    using CoordODim = typename base_type::CoordODim;

    using GridDDim = typename base_type::GridDDim;
    using GridODim = typename base_type::GridODim;

    using SplineBuilder2D = ddc::SplineBuilder2D<
            Kokkos::DefaultExecutionSpace,
            Kokkos::DefaultExecutionSpace::memory_space,
            typename base_type::BSplinesDDim,
            typename base_type::BSplinesODim,
            GridDDim,
            GridODim,
            base_type::SplineBoundary,
            base_type::SplineBoundary,
            base_type::SplineBoundary,
            base_type::SplineBoundary,
            ddc::SplineSolver::LAPACK,
            GridDDim,
            GridODim>;

    using SplineEvaluator2D = ddc::SplineEvaluator2D<
            Kokkos::DefaultExecutionSpace,
            Kokkos::DefaultExecutionSpace::memory_space,
            typename base_type::BSplinesDDim,
            typename base_type::BSplinesODim,
            GridDDim,
            GridODim,
            ddc::ConstantExtrapolationRule<DDim, ODim>,
            ddc::ConstantExtrapolationRule<DDim, ODim>,
            ddc::ConstantExtrapolationRule<ODim, DDim>,
            ddc::ConstantExtrapolationRule<ODim, DDim>,
            GridDDim,
            GridODim>;

    using IdxRangeDDim = typename base_type::IdxRangeDDim;
    using IdxRangeODim = typename base_type::IdxRangeODim;
    using IdxRangeFull = typename base_type::IdxRangeFull;

    using SplineInterpPointsDDim = typename base_type::SplineInterpPointsDDim;
    using SplineInterpPointsODim = typename base_type::SplineInterpPointsODim;

    ddc::ConstantExtrapolationRule<DDim, ODim> const m_bv_ddim_min;
    ddc::ConstantExtrapolationRule<DDim, ODim> const m_bv_ddim_max;
    ddc::ConstantExtrapolationRule<ODim, DDim> const m_bv_odim_min;
    ddc::ConstantExtrapolationRule<ODim, DDim> const m_bv_odim_max;

public:
    PartialDerivativeTestSpline2D(
            double const ddim_min,
            double const ddim_max,
            double const odim_min,
            double const odim_max)
        : base_type(ddim_min, ddim_max, odim_min, odim_max)
        , m_bv_ddim_min(CoordDDim(ddim_min), CoordODim(odim_min), CoordODim(odim_max))
        , m_bv_ddim_max(CoordDDim(ddim_max), CoordODim(odim_min), CoordODim(odim_max))
        , m_bv_odim_min(CoordODim(odim_min), CoordDDim(ddim_min), CoordDDim(ddim_max))
        , m_bv_odim_max(CoordODim(odim_max), CoordDDim(ddim_min), CoordDDim(ddim_max))
    {
    }

    double const operator()(double& delta_ddim) const
    {
        IdxRangeDDim const idxrange_ddim = SplineInterpPointsDDim::template get_domain<GridDDim>();
        IdxRangeODim const idxrange_odim = SplineInterpPointsODim::template get_domain<GridODim>();
        IdxRangeFull const idxrange = IdxRangeFull(idxrange_ddim, idxrange_odim);

        SplineBuilder2D builder(idxrange);
        SplineBuilder2DCache<SplineBuilder2D> builder_cache(builder);
        SplineEvaluator2D evaluator(m_bv_ddim_min, m_bv_ddim_max, m_bv_odim_min, m_bv_odim_max);

        Spline2DPartialDerivativeCreator<
                SplineBuilder2DCache<SplineBuilder2D>,
                SplineEvaluator2D,
                DDim> const derivative_creator(builder_cache, evaluator);

        using FunDiff = FunctionToDifferentiatePolynomial<typename base_type::DDim, spline_degree>;
        FunDiff function_to_differentiate;
        double const max_error = base_type::template compute_max_error<
                FunDiff>(function_to_differentiate, derivative_creator, delta_ddim);

        return max_error;
    }
};


/**
 * @brief A class that represents a test for partial derivatives.
 * The test can be used with finite difference method for
 * computing partial derivatives.
 */
template <
        class DerivativeDimension,
        std::size_t N_ddim,
        std::size_t N_odim,
        std::size_t spline_degree>
class PartialDerivativeTestFDM
    : public PartialDerivativeTest<DerivativeDimension, N_ddim, N_odim, spline_degree>
{
private:
    using base_type = PartialDerivativeTest<DerivativeDimension, N_ddim, N_odim, spline_degree>;
    using DDim = DerivativeDimension;

public:
    PartialDerivativeTestFDM(
            double const ddim_min,
            double const ddim_max,
            double const odim_min,
            double const odim_max)
        : base_type(ddim_min, ddim_max, odim_min, odim_max)
    {
    }

    double const operator()(double& delta_ddim) const
    {
        CentralFDMPartialDerivativeCreator<typename base_type::IdxRangeFull, DDim> const
                derivative_creator;

        using FunDiff = FunctionToDifferentiatePolynomial<typename base_type::DDim, spline_degree>;
        FunDiff function_to_differentiate;
        double const max_error = base_type::template compute_max_error<
                FunDiff>(function_to_differentiate, derivative_creator, delta_ddim);

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
    std::size_t constexpr spline_degree = 4;
    // relative error of convergence should be less than 15%
    double const TOL = 0.15;

    // Partial Derivative in X direction
    double delta_low_x,
            delta_high_x; // the maximum distance between points in the derivative direction
    PartialDerivativeTestSpline1D<X, 10, 10, spline_degree> const
            test_low_x(xmin, xmax, ymin, ymax);
    PartialDerivativeTestSpline1D<X, 100, 10, spline_degree> const
            test_high_x(xmin, xmax, ymin, ymax);
    double const error_low_x = test_low_x(delta_low_x);
    double const error_high_x = test_high_x(delta_high_x);


    double const order_x
            = std::log(error_high_x / error_low_x) / std::log(delta_high_x / delta_low_x);
    double const relative_error_order_x = std::fabs((spline_degree - order_x) / spline_degree);

    EXPECT_LE(relative_error_order_x, TOL);

    // Partial Derivative in Y direction
    double delta_low_y,
            delta_high_y; // the maximum distance between points in the derivative direction
    PartialDerivativeTestSpline1D<Y, 10, 10, spline_degree> const
            test_low_y(ymin, ymax, xmin, xmax);
    PartialDerivativeTestSpline1D<Y, 100, 10, spline_degree> const
            test_high_y(ymin, ymax, xmin, xmax);
    double const error_low_y = test_low_y(delta_low_y);
    double const error_high_y = test_high_y(delta_high_y);

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
    std::size_t constexpr spline_degree = 4;
    // relative error of convergence should be less than 15%
    double const TOL = 0.15;

    // Partial Derivative in X direction
    double delta_low_x,
            delta_high_x; // the maximum distance between points in the derivative direction
    PartialDerivativeTestSpline2D<X, 10, 10, spline_degree> const
            test_low_x(xmin, xmax, ymin, ymax);
    PartialDerivativeTestSpline2D<X, 100, 10, spline_degree> const
            test_high_x(xmin, xmax, ymin, ymax);
    double const error_low_x = test_low_x(delta_low_x);
    double const error_high_x = test_high_x(delta_high_x);


    double const order_x
            = std::log(error_high_x / error_low_x) / std::log(delta_high_x / delta_low_x);
    double const relative_error_order_x = std::fabs((spline_degree - order_x) / spline_degree);

    EXPECT_LE(relative_error_order_x, TOL);

    // Partial Derivative in Y direction
    double delta_low_y,
            delta_high_y; // the maximum distance between points in the derivative direction
    PartialDerivativeTestSpline2D<Y, 10, 10, spline_degree> const
            test_low_y(ymin, ymax, xmin, xmax);
    PartialDerivativeTestSpline2D<Y, 100, 10, spline_degree> const
            test_high_y(ymin, ymax, xmin, xmax);
    double const error_low_y = test_low_y(delta_low_y);
    double const error_high_y = test_high_y(delta_high_y);

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
    std::size_t constexpr spline_degree = 4;
    // relative error of convergence should be less than 15%
    double const TOL = 0.15;

    // Partial Derivative in X direction
    double delta_low_x,
            delta_high_x; // the maximum distance between points in the derivative direction
    PartialDerivativeTestFDM<X, 10, 10, spline_degree> const test_low_x(xmin, xmax, ymin, ymax);
    PartialDerivativeTestFDM<X, 100, 10, spline_degree> const test_high_x(xmin, xmax, ymin, ymax);
    double const error_low_x = test_low_x(delta_low_x);
    double const error_high_x = test_high_x(delta_high_x);

    int const expected_order(2);
    double const order_x
            = std::log(error_high_x / error_low_x) / std::log(delta_high_x / delta_low_x);
    double const relative_error_order_x = std::fabs((expected_order - order_x) / spline_degree);

    EXPECT_LE(relative_error_order_x, TOL);

    // Partial Derivative in Y direction
    double delta_low_y,
            delta_high_y; // the maximum distance between points in the derivative direction
    PartialDerivativeTestFDM<Y, 10, 10, spline_degree> const test_low_y(ymin, ymax, xmin, xmax);
    PartialDerivativeTestFDM<Y, 100, 10, spline_degree> const test_high_y(ymin, ymax, xmin, xmax);
    double const error_low_y = test_low_y(delta_low_y);
    double const error_high_y = test_high_y(delta_high_y);

    double const order_y
            = std::log(error_high_y / error_low_y) / std::log(delta_high_y / delta_low_y);
    double const relative_error_order_y = std::fabs((expected_order - order_y) / spline_degree);

    EXPECT_LE(relative_error_order_y, TOL);
}
} // namespace
