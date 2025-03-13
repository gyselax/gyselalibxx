// SPDX-License-Identifier: MIT
#include <ddc/ddc.hpp>
#include <ddc/kernels/splines.hpp>

#include <gtest/gtest.h>

#include "central_fdm_partial_derivatives.hpp"
#include "ddc_aliases.hpp"
#include "math_tools.hpp"
#include "mesh_builder.hpp"
#include "spline_1d_partial_derivative.hpp"

namespace {
/**
 * @brief A class that represents a polynomial test function for computing partial derivatives.
 * The polynomial depends on two variables.
 */
template<class DimA, class DimB>
class FunctionToDifferentiatePolynomial
{
private:
  using CoordAB = Coord<DimA, DimB>;
public:
    /**
     * @brief Get the value of the function at given coordinate.
     *
     * @param[in] coord_ab The coordinate where we want to evaluate
     * the function.
     *
     * @return The value of the function at the coordinate.
     */
    KOKKOS_FUNCTION double operator()(CoordAB const coord_ab) const
    {
        double const xa = ddc::get<DimA>(coord_ab);
        double const xb = ddc::get<DimB>(coord_ab);
        return (ipow(xa, 3) + ipow(xb, 2)) * xb;
    }

    /**
     * @brief Get the value of the partial derivative of the function
     * in the DerivativeDimension direction at a given coordinate.
     *
     * @param[in] coord_ab The coordinate where we want to evaluate 
     * the partial derivative.
     *
     * @return The value of the partial derivative of the function
     * at the coordinate.
     */
    template <class DerivativeDimension>
    KOKKOS_FUNCTION double differentiate(CoordAB const coord_ab) const
    {
        static_assert(
                std::is_same_v<DerivativeDimension, DimA> || std::is_same_v<DerivativeDimension, DimB>);
        double const xa = ddc::get<DimA>(coord_ab);
        double const xb = ddc::get<DimB>(coord_ab);

        if constexpr (std::is_same_v<DerivativeDimension, DimA>) {
            return 3. * ipow(xa, 2) * xb;
        } else {
            return ipow(xa, 3) + 3. * ipow(xb, 2);
        }
    }
};


template<std::size_t Nx, std::size_t Ny>
class PartialDerivativeTest
{
public:

  struct GridX;
struct X
{
    static bool constexpr PERIODIC = false;
};
using CoordX = Coord<X>;

struct BSplinesX : ddc::UniformBSplines<X, 3>
{
};
auto static constexpr SplineXBoundary = ddc::BoundCond::GREVILLE;

using SplineInterpPointsX
        = ddc::GrevilleInterpolationPoints<BSplinesX, SplineXBoundary, SplineXBoundary>;

struct GridX : NonUniformGridBase<X>
{
};

using IdxX = Idx<GridX>;
using IdxStepX = IdxStep<GridX>;
using IdxRangeX = IdxRange<GridX>;

struct Y
{
    static bool constexpr PERIODIC = false;
};
using CoordY = Coord<Y>;
struct BSplinesY : ddc::UniformBSplines<Y, 3>
{
};
auto static constexpr SplineYBoundary = ddc::BoundCond::GREVILLE;
using SplineInterpPointsY
        = ddc::GrevilleInterpolationPoints<BSplinesY, SplineYBoundary, SplineYBoundary>;
struct GridY : SplineInterpPointsY::interpolation_discrete_dimension_type
{
};
using IdxY = Idx<GridY>;
using IdxStepY = IdxStep<GridY>;
using IdxRangeY = IdxRange<GridY>;

using IdxXY = Idx<GridX, GridY>;
using IdxRangeXY = IdxRange<GridX, GridY>;
using DFieldMemXY = DFieldMem<IdxRangeXY>;
using DFieldXY = DField<IdxRangeXY>;
using DConstFieldXY = DConstField<IdxRangeXY>;

using CoordXY = Coord<X, Y>;

// --- Operators ---
using SplineXBuilder = ddc::SplineBuilder<
        Kokkos::DefaultExecutionSpace,
        Kokkos::DefaultExecutionSpace::memory_space,
        BSplinesX,
        GridX,
        SplineXBoundary,
        SplineXBoundary,
        ddc::SplineSolver::LAPACK,
        GridX,
        GridY>;
using SplineXEvaluator = ddc::SplineEvaluator<
        Kokkos::DefaultExecutionSpace,
        Kokkos::DefaultExecutionSpace::memory_space,
        BSplinesX,
        GridX,
        ddc::ConstantExtrapolationRule<X>,
        ddc::ConstantExtrapolationRule<X>,
        GridX,
        GridY>;

// --- Operators ---
using SplineYBuilder = ddc::SplineBuilder<
        Kokkos::DefaultExecutionSpace,
        Kokkos::DefaultExecutionSpace::memory_space,
        BSplinesY,
        GridY,
        SplineYBoundary,
        SplineYBoundary,
        ddc::SplineSolver::LAPACK,
        GridX,
        GridY>;
using SplineYEvaluator = ddc::SplineEvaluator<
        Kokkos::DefaultExecutionSpace,
        Kokkos::DefaultExecutionSpace::memory_space,
        BSplinesY,
        GridY,
        ddc::ConstantExtrapolationRule<Y>,
        ddc::ConstantExtrapolationRule<Y>,
        GridX,
        GridY>;


private:
    CoordX const m_xmin;
    CoordX const m_xmax;
    CoordY const m_ymin;
    CoordY const m_ymax;

    IdxStepX const m_ncellsx;
    IdxStepY const m_ncellsy;
public:

    IdxRangeXY m_idxrange_xy; 

  PartialDerivativeTest(CoordX xmin, CoordX xmax, CoordY ymin, CoordY ymax)
  : m_xmin(xmin)
  , m_xmax(xmax)
  , m_ymin(ymin)
  , m_ymax(ymax)
  , m_ncellsx(Nx)
  , m_ncellsy(Ny)
  {
    ddc::init_discrete_space<BSplinesX>(m_xmin, m_xmax, m_ncellsx);
    ddc::init_discrete_space<GridX>(SplineInterpPointsX::template get_sampling<GridX>());
    IdxRangeX idxrange_x(SplineInterpPointsX::template get_domain<GridX>());

    ddc::init_discrete_space<BSplinesY>(m_ymin, m_ymax, m_ncellsy);
    ddc::init_discrete_space<GridY>(SplineInterpPointsY::template get_sampling<GridY>());
    IdxRangeY idxrange_y(SplineInterpPointsY::template get_domain<GridY>());

    IdxRangeXY idxrange_xy(idxrange_x, idxrange_y);
    m_idxrange_xy = idxrange_xy;
      }

  template <class DerivativeDimension, class FunctionToDifferentiate>
   double operator()(FunctionToDifferentiate const& function_to_differentiate,
    IPartialDerivativeCreator<IdxRangeXY, DerivativeDimension> const& partial_derivative_creator) const
{
    // field to be differentiated
    DFieldMemXY field_to_differentiate(m_idxrange_xy);
    DFieldXY field_to_differentiate_proxy = get_field(field_to_differentiate);
    ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            m_idxrange_xy,
            KOKKOS_LAMBDA(IdxXY const idxy) {
                field_to_differentiate_proxy(idxy) = function_to_differentiate(ddc::coordinate(idxy));
            });

    std::unique_ptr<IPartialDerivative<IdxRangeXY, DerivativeDimension>> const
            partial_derivative_creator_pointer
            = partial_derivative_creator.create_instance(get_const_field(field_to_differentiate));

    IPartialDerivative<IdxRangeXY, DerivativeDimension> const& partial_derivative
            = *partial_derivative_creator_pointer;

    DFieldMemXY field_differentiated_alloc(m_idxrange_xy);
    DFieldXY field_differentiated = get_field(field_differentiated_alloc);
    partial_derivative(field_differentiated);

    double const max_error = ddc::parallel_transform_reduce(
            Kokkos::DefaultExecutionSpace(),
            m_idxrange_xy,
            0.,
            ddc::reducer::max<double>(),
            KOKKOS_LAMBDA(IdxXY const idxy) {
                return Kokkos::abs(field_differentiated(idxy) - function_to_differentiate.template differentiate<DerivativeDimension>(ddc::coordinate(idxy)));
            });

    return max_error;
}
};

double test_spline_partial)_derivative(std::size_t ncells_x, std::size_t ncells_y)


TEST(PartialDerivative, Spline1DPartialDerivative)
{
    std::vector<std::size_t> nbcells_x_list = {10, 100};
    std::vector<double> error_inf;

    for ( std::size_t ncells_x : nbcells_x_list) {

  
    std::size_t ncells_y(10);
    using PTest = PartialDerivativeTest<ncells_x, ncells_y>;
    PTest::CoordX const xmin(0.1);
    PTest::CoordX const xmax(1);
    PTest::CoordY const ymin(0.2);
    PTest::CoordY const ymax(2);
    PTest const derivative_test_x_spline1d(xmin, xmax, ymin, ymax);

    // partial derivatives evaluated with 1d splines in X direction
    PTest::SplineXBuilder const builder_x(derivative_test_x_spline1d.m_idxrange_xy);

    ddc::ConstantExtrapolationRule<PTest::X> bv_x_min(xmin);
    ddc::ConstantExtrapolationRule<PTest::X> bv_x_max(xmax);
    PTest::SplineXEvaluator const spline_evaluator_x(bv_x_min, bv_x_max);

    Spline1DPartialDerivativeCreator<typename PTest::SplineXBuilder, typename PTest::SplineXEvaluator> const
            partial_dx_creator(builder_x, spline_evaluator_x);

    FunctionToDifferentiatePolynomial<PTest::X, PTest::Y> function_to_differentiate;
    error_inf.emplace_back(derivative_test_x_spline1d.operator()<PTest::X, 
    FunctionToDifferentiatePolynomial<PTest::X, PTest::Y>>(function_to_differentiate, partial_dx_creator));
  }
    int const order = std::log(error_inf.front() / error_inf.back());
    double const relative_error = std::fabs(
            (error_inf.back() * pow(10., order) - error_inf.front()) / error_inf.front());
    // relative error should be less than 5%
    double const TOL = 0.05;
    EXPECT_LE(relative_error, TOL);

    // partial derivatives evaluated with 1d splines in Y direction
    //SplineYBuilder const builder_y(idxrange_xy);
    //ddc::ConstantExtrapolationRule<Y> bv_y_min(y_min);
    //ddc::ConstantExtrapolationRule<Y> bv_y_max(y_max);
    //SplineYEvaluator const spline_evaluator_y(bv_y_min, bv_y_max);

    //Spline1DPartialDerivativeCreator<SplineYBuilder, SplineYEvaluator> const
    //        partial_dy_creator(builder_y, spline_evaluator_y);
    //PartialDerivativeTest<10, 10, FunctionToDifferentiatePolynomial, Y> const derivative_test_y_spline1d(xmin, xmax, ymin, ymax, 
    //                                                                                                     partial_dy_creator);
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
