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
using CoordX = Coord<X>;
using CoordY = Coord<Y>;
using CoordXY = Coord<X, Y>;

/**
 * @brief A class that represents a polynomial test function for computing partial derivatives.
 * The polynomial depends on two variables.
 */
template<class DerivativeDimension>
class FunctionToDifferentiatePolynomial
{
    static_assert(std::is_same_v<DerivativeDimension,X> || std::is_same_v<DerivativeDimension,Y>);

public:
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
        return (ipow(x, 6) + ipow(y, 5)) * y;
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
    KOKKOS_FUNCTION double differentiate(CoordXY const coord_xy) const
    {
        double const x = ddc::get<X>(coord_xy);
        double const y = ddc::get<Y>(coord_xy);

        if constexpr (std::is_same_v<DerivativeDimension, X>) {
            return 6. * ipow(x, 5) * y;
        } else {
            return ipow(x, 6) + 5. * ipow(y, 4);
        }
    }
};


template<class DerivativeDimension, std::size_t N_ddim, std::size_t spline_degree>
class PartialDerivativeTest
{
public:

struct BSplinesX : ddc::UniformBSplines<X, spline_degree>
{
};
auto static constexpr SplineBoundary = ddc::BoundCond::GREVILLE;

using SplineInterpPointsX
        = ddc::GrevilleInterpolationPoints<BSplinesX, SplineBoundary, SplineBoundary>;

struct GridX : NonUniformGridBase<X>
{
};


struct BSplinesY : ddc::UniformBSplines<Y, spline_degree>
{
};
using SplineInterpPointsY
        = ddc::GrevilleInterpolationPoints<BSplinesY, SplineBoundary, SplineBoundary>;
struct GridY : SplineInterpPointsY::interpolation_discrete_dimension_type
{
};

using IdxRangeX = IdxRange<GridX>;
using IdxRangeY = IdxRange<GridY>;
using IdxRangeXY = IdxRange<GridX, GridY>;

using IdxXY = Idx<GridX, GridY>;

using DFieldMemXY = DFieldMem<IdxRangeXY>;
using DFieldXY = DField<IdxRangeXY>;


  template <class FunctionToDifferentiate>
   double compute_max_error(IdxRangeXY const& idxrange_xy, FunctionToDifferentiate const& function_to_differentiate,
    IPartialDerivativeCreator<IdxRangeXY, DerivativeDimension> const& partial_derivative_creator) const
{
    // field to be differentiated
    DFieldMemXY field_to_differentiate(idxrange_xy);
    DFieldXY field_to_differentiate_proxy = get_field(field_to_differentiate);
    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            idxrange_xy,
            KOKKOS_LAMBDA(IdxXY const idxy) {
                field_to_differentiate_proxy(idxy) = function_to_differentiate(ddc::coordinate(idxy));
            });

    std::unique_ptr<IPartialDerivative<IdxRangeXY, DerivativeDimension>> const
            partial_derivative_creator_pointer
            = partial_derivative_creator.create_instance(get_const_field(field_to_differentiate));

    IPartialDerivative<IdxRangeXY, DerivativeDimension> const& partial_derivative
            = *partial_derivative_creator_pointer;

    DFieldMemXY field_differentiated_alloc(idxrange_xy);
    DFieldXY field_differentiated = get_field(field_differentiated_alloc);
    partial_derivative(field_differentiated);

    double const max_error = ddc::parallel_transform_reduce(
            Kokkos::DefaultExecutionSpace(),
            idxrange_xy,
            0.,
            ddc::reducer::max<double>(),
            KOKKOS_LAMBDA(IdxXY const idxy) {
                return Kokkos::abs(field_differentiated(idxy) - function_to_differentiate.differentiate(ddc::coordinate(idxy)));
            });

    return max_error;
}
};

template <class DerivativeDimension, std::size_t N_ddim, std::size_t N_odim, std::size_t spline_degree>
class PartialDerivativeTestSpline1D : public PartialDerivativeTest<DerivativeDimension, N_ddim, spline_degree>
{
private:
  using base_type = PartialDerivativeTest<DerivativeDimension, N_ddim, spline_degree>;
  using DDim = DerivativeDimension;
  using ODim = std::conditional_t<std::is_same_v<DDim, X>, Y, X>;

  using CoordDDim = Coord<DDim>;
  using CoordODim = Coord<ODim>;

  using GridDDim = std::conditional_t<std::is_same_v<DDim, X>, typename base_type::GridX, typename base_type::GridY>;
  using GridODim = std::conditional_t<std::is_same_v<DDim, X>, typename base_type::GridY, typename base_type::GridX>;

  using IdxRangeDDim = IdxRange<GridDDim>;
  using IdxRangeODim = IdxRange<GridODim>;
  using IdxRangeFull = IdxRange<GridDDim, GridODim>;

  using IdxStepDDim = IdxStep<GridDDim>;
  using IdxStepODim = IdxStep<GridODim>;

  using BSplinesDDim = std::conditional_t<std::is_same_v<DDim, X>, typename base_type::BSplinesX, typename base_type::BSplinesY>;
  using BSplinesODim = std::conditional_t<std::is_same_v<DDim, X>, typename base_type::BSplinesY, typename base_type::BSplinesX>;

  using SplineInterpPointsDDim = std::conditional_t<std::is_same_v<DDim, X>, typename base_type::SplineInterpPointsX, 
  typename base_type::SplineInterpPointsY>;
  using SplineInterpPointsODim = std::conditional_t<std::is_same_v<DDim, X>, typename base_type::SplineInterpPointsY, 
  typename base_type::SplineInterpPointsX>;

using SplineDDimBuilder = ddc::SplineBuilder<
        Kokkos::DefaultExecutionSpace,
        Kokkos::DefaultExecutionSpace::memory_space,
        BSplinesDDim,
        GridDDim,
        ddc::BoundCond::GREVILLE,
        ddc::BoundCond::GREVILLE,
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

    CoordDDim const m_ddim_min;
    CoordDDim const m_ddim_max;
    CoordODim const m_odim_min;
    CoordODim const m_odim_max; 
    IdxStepDDim const m_ncells_ddim;
    IdxStepODim const m_ncells_odim;

public:

PartialDerivativeTestSpline1D(double const ddim_min, double const ddim_max, double const odim_min, double const odim_max)
    : m_ddim_min(ddim_min)
    , m_ddim_max(ddim_max)
    , m_odim_min(odim_min)
    , m_odim_max(odim_max)
    , m_ncells_ddim(N_ddim)
    , m_ncells_odim(N_odim)
  {
    ddc::init_discrete_space<BSplinesDDim>(m_ddim_min, m_ddim_max, m_ncells_ddim);
    ddc::init_discrete_space<GridDDim>(SplineInterpPointsDDim::template get_sampling<GridDDim>());

    ddc::init_discrete_space<BSplinesODim>(m_odim_min, m_odim_max, m_ncells_odim);
    ddc::init_discrete_space<GridODim>(SplineInterpPointsODim::template get_sampling<GridODim>());
}

double const operator()(double& delta_ddim) const {

    IdxRangeDDim idxrange_ddim(SplineInterpPointsDDim::template get_domain<GridDDim>());
    IdxRangeODim idxrange_odim(SplineInterpPointsODim::template get_domain<GridODim>());

    
    IdxRangeFull const idxrange(idxrange_ddim, idxrange_odim);
    IdxRangeDDim const idxrange_derivative(idxrange); 
    delta_ddim = ddcHelper::maximum_distance_between_two_points(idxrange_derivative);
    
    SplineDDimBuilder const builder(idxrange);

    ddc::ConstantExtrapolationRule<DDim> bv_min(m_ddim_min);
    ddc::ConstantExtrapolationRule<DDim> bv_max(m_ddim_max);
    SplineDDimEvaluator const spline_evaluator(bv_min, bv_max);
    
    Spline1DPartialDerivativeCreator<SplineDDimBuilder, SplineDDimEvaluator> const
            derivative_creator(builder, spline_evaluator);

    FunctionToDifferentiatePolynomial<DDim> function_to_differentiate;
    double const max_error = base_type::template compute_max_error<FunctionToDifferentiatePolynomial<DDim>>(
    idxrange, function_to_differentiate, derivative_creator);

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
TEST(PartialDerivative, Spline1DPartialDerivative1D)
{
    double const xmin(0.1);
    double const xmax(1.1);
    double const ymin(0.2);
    double const ymax(2.);
    std::size_t constexpr spline_degree = 4;
    // relative error of convergence should be less than 5%
    double const TOL = 0.05;

    // Partial Derivative in X direction 
    double delta_low_x, delta_high_x; // the maximum distance between points in the derivative direction
    PartialDerivativeTestSpline1D<X, 10, 10, spline_degree> const test_low_x(xmin, xmax, ymin, ymax);
    PartialDerivativeTestSpline1D<X, 100, 10, spline_degree> const test_high_x(xmin, xmax, ymin, ymax);
    double const error_low_x = test_low_x(delta_low_x);
    double const error_high_x = test_high_x(delta_high_x);

    
    double const order_x = std::log(error_high_x / error_low_x) / std::log(delta_high_x / delta_low_x);
    double const relative_error_order_x = std::fabs((spline_degree-order_x) / spline_degree);

    EXPECT_LE(relative_error_order_x, TOL);
    

    // Partial Derivative in Y direction 
    double delta_low_y, delta_high_y; // the maximum distance between points in the derivative direction
    PartialDerivativeTestSpline1D<Y, 10, 10, spline_degree> const test_low_y(ymin, ymax, xmin, xmax);
    PartialDerivativeTestSpline1D<Y, 100, 10, spline_degree> const test_high_y(ymin, ymax, xmin, xmax);
    double const error_low_y = test_low_y(delta_low_y);
    double const error_high_y = test_low_y(delta_high_y);
    
    double const order_y = std::log(error_high_y / error_low_y) / std::log(delta_high_y / delta_low_y);
    double const relative_error_order_y = std::fabs((spline_degree-order_y) / spline_degree);

    EXPECT_LE(relative_error_order_y, TOL);
}

//TEST(PartialDerivative, CentralFDMPartialDerivativeDx)
//{
//    int n_elems_x(10);
//    int n_elems_y(20);
//
//    Coord<X> const x_min(0.0);
//    Coord<X> const x_max(1.0);
//    IdxStepX x_ncells(n_elems_x);
//
//    Coord<Y> const y_min(0.0);
//    Coord<Y> const y_max(2.0);
//    IdxStepY y_ncells(n_elems_y);
//
//    ddc::init_discrete_space<GridX>(build_random_non_uniform_break_points(x_min, x_max, x_ncells));
//    IdxRangeX idxrange_x(IdxX {0}, x_ncells);
//
//    ddc::init_discrete_space<GridY>(build_random_non_uniform_break_points(y_min, y_max, y_ncells));
//    IdxRangeY idxrange_y(IdxY {0}, y_ncells);
//
//    IdxRangeXY idxrange_xy(idxrange_x, idxrange_y);
//
//    CentralFDMPartialDerivativeCreator<IdxRangeXY, X> const partial_dx_creator;
//
//    FunctionToDifferentiatePolynomial function_to_differentiate;
//}
} // namespace
