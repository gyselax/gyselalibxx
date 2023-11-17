#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>

#include <ddc/ddc.hpp>

#include <sll/bsplines_non_uniform.hpp>
#include <sll/constant_extrapolation_boundary_value.hpp>
#include <sll/greville_interpolation_points.hpp>
#include <sll/null_boundary_value.hpp>
#include <sll/spline_builder_2d.hpp>
#include <sll/spline_evaluator_2d.hpp>

#include <math.h>
#include <paraconf.h>
#include <pdi.h>
#include <stdio.h>

#include "test_utils.hpp"



namespace {
struct DimX
{
    static bool constexpr PERIODIC = false;
};
struct DimY
{
    static bool constexpr PERIODIC = false;
};
struct DimR
{
    static bool constexpr PERIODIC = false;
};

struct DimP
{
    static bool constexpr PERIODIC = true;
};

int constexpr BSDegree = 3;

// Polar dimensions
using CoordR = ddc::Coordinate<DimR>;
using CoordP = ddc::Coordinate<DimP>;
using CoordRP = ddc::Coordinate<DimR, DimP>;

using BSplinesR = NonUniformBSplines<DimR, BSDegree>;
using BSplinesP = NonUniformBSplines<DimP, BSDegree>;

using InterpPointsR
        = GrevilleInterpolationPoints<BSplinesR, BoundCond::GREVILLE, BoundCond::GREVILLE>;
using InterpPointsP
        = GrevilleInterpolationPoints<BSplinesP, BoundCond::PERIODIC, BoundCond::PERIODIC>;

using IDimR = typename InterpPointsR::interpolation_mesh_type;
using IDimP = typename InterpPointsP::interpolation_mesh_type;

using SplineRBuilder = SplineBuilder<BSplinesR, IDimR, BoundCond::GREVILLE, BoundCond::GREVILLE>;
using SplinePBuilder = SplineBuilder<BSplinesP, IDimP, BoundCond::PERIODIC, BoundCond::PERIODIC>;
using SplineRPBuilder = SplineBuilder2D<SplineRBuilder, SplinePBuilder>;

using SplineRPEvaluator = SplineEvaluator2D<BSplinesR, BSplinesP>;

using BSDomainR = ddc::DiscreteDomain<BSplinesR>;
using BSDomainP = ddc::DiscreteDomain<BSplinesP>;
using BSDomainRP = ddc::DiscreteDomain<BSplinesR, BSplinesP>;

using IDomainR = ddc::DiscreteDomain<IDimR>;
using IDomainP = ddc::DiscreteDomain<IDimP>;
using IDomainRP = ddc::DiscreteDomain<IDimR, IDimP>;

using IndexR = ddc::DiscreteElement<IDimR>;
using IndexP = ddc::DiscreteElement<IDimP>;
using IndexRP = ddc::DiscreteElement<IDimR, IDimP>;

using IVectR = ddc::DiscreteVector<IDimR>;
using IVectP = ddc::DiscreteVector<IDimP>;
using IVectRP = ddc::DiscreteVector<IDimR, IDimP>;

template <class ElementType>
using SpanR = ddc::ChunkSpan<ElementType, IDomainR>;

template <class ElementType>
using SpanP = ddc::ChunkSpan<ElementType, IDomainP>;

template <class ElementType>
using SpanRP = ddc::ChunkSpan<ElementType, IDomainRP>;

using DSpanR = SpanR<double>;
using DSpanP = SpanP<double>;
using DSpanRP = SpanRP<double>;

using IDomainRP = ddc::DiscreteDomain<IDimR, IDimP>;


template <class ElementType>
using FieldRP = ddc::Chunk<ElementType, IDomainRP>;
using DFieldRP = FieldRP<double>;



// Cartesian dimensions
using CoordX = ddc::Coordinate<DimX>;
using CoordY = ddc::Coordinate<DimY>;
using CoordXY = ddc::Coordinate<DimX, DimY>;

using BSplinesX = NonUniformBSplines<DimX, BSDegree>;
using BSplinesY = NonUniformBSplines<DimY, BSDegree>;


auto constexpr SplineXBoundary = DimX::PERIODIC ? BoundCond::PERIODIC : BoundCond::GREVILLE;
using InterpPointsX = GrevilleInterpolationPoints<BSplinesX, SplineXBoundary, SplineXBoundary>;
using IDimX = typename InterpPointsX::interpolation_mesh_type;
using SplineXBuilder = SplineBuilder<BSplinesX, IDimX, SplineXBoundary, SplineXBoundary>;

auto constexpr SplineYBoundary = DimY::PERIODIC ? BoundCond::PERIODIC : BoundCond::GREVILLE;
using InterpPointsY = GrevilleInterpolationPoints<BSplinesY, SplineYBoundary, SplineYBoundary>;
using IDimY = typename InterpPointsY::interpolation_mesh_type;
using SplineYBuilder = SplineBuilder<BSplinesY, IDimY, SplineYBoundary, SplineYBoundary>;


using SplineXYBuilder = SplineBuilder2D<SplineXBuilder, SplineYBuilder>;
using SplineXYEvaluator = SplineEvaluator2D<BSplinesX, BSplinesY>;

using BSDomainX = ddc::DiscreteDomain<BSplinesX>;
using BSDomainY = ddc::DiscreteDomain<BSplinesY>;
using BSDomainXY = ddc::DiscreteDomain<BSplinesX, BSplinesY>;

using IDomainX = ddc::DiscreteDomain<IDimX>;
using IDomainY = ddc::DiscreteDomain<IDimY>;
using IDomainXY = ddc::DiscreteDomain<IDimX, IDimY>;

using IndexX = ddc::DiscreteElement<IDimX>;
using IndexY = ddc::DiscreteElement<IDimY>;
using IndexXY = ddc::DiscreteElement<IDimX, IDimY>;

using IVectX = ddc::DiscreteVector<IDimX>;
using IVectY = ddc::DiscreteVector<IDimY>;
using IVectXY = ddc::DiscreteVector<IDimX, IDimY>;

template <class ElementType>
using SpanX = ddc::ChunkSpan<ElementType, IDomainX>;

template <class ElementType>
using SpanY = ddc::ChunkSpan<ElementType, IDomainY>;

template <class ElementType>
using SpanXY = ddc::ChunkSpan<ElementType, IDomainXY>;

using DSpanX = SpanX<double>;
using DSpanY = SpanY<double>;
using DSpanXY = SpanXY<double>;

using IDomainXY = ddc::DiscreteDomain<IDimX, IDimY>;


template <class ElementType>
using FieldXY = ddc::Chunk<ElementType, IDomainXY>;
using DFieldXY = FieldXY<double>;



/**
 * @brief Check in the polar domain if the maximum absolute value of the difference between
 * the exact value and the computed value is below the tolerance error given.
 *
 * @param[in] function_evaluated
 * 			The evaluated function on B-splines on a polar grid. .
 * @param[in] exact_function
 * 			The exact function PolarExactFunction type.
 * @param[in] coords
 * 			The polar coordinates of the grid used to define the discrete space.
 * @param[in] outside_coords
 * 			The polar coordinates where we evaluate the function. All points are supposed to be outside the domain.
 * @param[in] TOL
 *			The tolerance error on the maximum absolute value of the difference between the exact value
 * 			and the computed value.
 */
template <class Function>
void check_constant_outside_domain(
        DSpanRP function_evaluated,
        Function& exact_function,
        SpanRP<CoordRP> coords,
        SpanRP<CoordRP> outside_coords,
        double const TOL)
{
    auto r_domain = ddc::get_domain<IDimR>(coords);
    double r_max = ddc::coordinate(r_domain.back());
    double max_err(0.0);
    for_each(outside_coords.domain(), [&](IndexRP const irp) {
        CoordRP coords_edge(r_max, ddc::get<DimP>(outside_coords(irp)));
        const double err = fabs(function_evaluated(irp) - exact_function(coords_edge));
        max_err = max_err > err ? max_err : err;
    });
    EXPECT_NEAR(max_err, 0., TOL);
}


/**
 * @brief Check in the Cartesian domain if the maximum absolute value of the difference between
 * the exact value and the computed value is below the tolerance error given.
 *
 * @param[in] function_evaluated
 * 			The evaluated function on B-splines on a Cartesian grid.
 * @param[in] exact_function
 * 			The exact function CartesianExactFunction type.
 * @param[in] coords
 * 			The Cartesian coordinates of the grid used to define the discrete space.
 * @param[in] outside_coords
 * 			The Cartesian coordinates where we evaluate the function. All points are supposed to be outside the domain.
 * @param[in] TOL
 *			The tolerance error on the maximum absolute value of the difference between the exact value
 * 			and the computed value.
 */
template <class Function>
void check_constant_outside_domain(
        DSpanXY function_evaluated,
        Function& exact_function,
        SpanXY<CoordXY> coords,
        SpanXY<CoordXY> outside_coords,
        double const TOL)
{
    auto x_domain = ddc::get_domain<IDimX>(coords);
    auto y_domain = ddc::get_domain<IDimY>(coords);
    const double x_max = ddc::coordinate(x_domain.back());
    const double x_min = ddc::coordinate(x_domain.front());
    const double y_max = ddc::coordinate(y_domain.back());
    const double y_min = ddc::coordinate(y_domain.front());
    double max_err(0.0);
    for_each(outside_coords.domain(), [&](IndexXY const ixy) {
        const double x = ddc::get<DimX>(outside_coords(ixy));
        const double y = ddc::get<DimY>(outside_coords(ixy));
        const double coord_x = std::max(x_min, std::min(x_max, x));
        const double coord_y = std::max(y_min, std::min(y_max, y));
        CoordXY coords_edge(coord_x, coord_y);
        const double err = fabs(function_evaluated(ixy) - exact_function(coords_edge));
        max_err = max_err > err ? max_err : err;
    });
    EXPECT_NEAR(max_err, 0., TOL);
}

/**
 * @brief Build a polar grid with points outside the domain.
 *
 * On the @f& \theta @f$ dimension, the points are not on the mesh points.
 * On the @f$ r @f$ dimension, all of the points are outside the domain.
 *
 * @param[in] grid
 * 		A polar grid used to define the discrete space.
 * @param[out] outside_coords
 * 		The built polar grid with points outside the domain.
 */
void build_outside_grid(IDomainRP const& grid, FieldRP<CoordRP>& outside_coords)
{
    auto r_domain = ddc::get_domain<IDimR>(outside_coords);
    IndexR const ir_max(r_domain.back());
    IndexR const ir_min(r_domain.front());

    auto theta_domain = ddc::get_domain<IDimP>(outside_coords);
    IndexP const ip_min(theta_domain.front());
    IndexP const ip_max(theta_domain.back());

    CoordR const r_min(ddc::coordinate(ir_min));
    CoordR const r_max(ddc::coordinate(ir_max));

    for_each(outside_coords.domain(), [&](IndexRP const irp) {
        CoordR const coord_r(coordinate(ddc::select<IDimR>(irp)));
        CoordP const coord_p(coordinate(ddc::select<IDimP>(irp)));

        IndexR const ir(ddc::select<IDimR>(irp));
        IndexP const ip(ddc::select<IDimP>(irp));


        CoordR delta_coord_r;
        if (ir + 1 <= ir_max) {
            delta_coord_r = ddc::coordinate(ir + 1) - ddc::coordinate(ir);
        } else {
            delta_coord_r = CoordR(0.);
        }

        CoordP delta_coord_p;
        if (ip + 1 <= ip_max) {
            delta_coord_p = ddc::coordinate(ip + 1) - ddc::coordinate(ip);
        } else {
            delta_coord_p = ddc::coordinate(ip_min + 1) - ddc::coordinate(ip_min);
        }

        outside_coords(irp) = CoordRP(
                double(coord_r + delta_coord_r + (r_max - r_min)),
                fmod(double(coord_p + delta_coord_p * 1.3), 2 * M_PI));
    });
}


/**
 * @brief Build a Cartesian grid with points outside the domain.
 *
 * On the @f& x @f$ and the @f$ y @f$ dimensions, all of the points are outside the domain.
 *
 * @param[in] grid
 * 		A Cartesian grid used to define the discrete space.
 * @param[out] outside_coords
 * 		The built Cartesian grid with points outside the domain.
 */
void build_outside_grid(IDomainXY const& grid, FieldXY<CoordXY>& outside_coords)
{
    auto x_domain = ddc::get_domain<IDimX>(outside_coords);
    IndexX const ix_max(x_domain.back());
    IndexX const ix_min(x_domain.front());

    auto theta_domain = ddc::get_domain<IDimY>(outside_coords);
    IndexY const iy_min(theta_domain.front());
    IndexY const iy_max(theta_domain.back());

    CoordX const x_min(ddc::coordinate(ix_min));
    CoordX const x_max(ddc::coordinate(ix_max));

    CoordY const y_min(ddc::coordinate(iy_min));
    CoordY const y_max(ddc::coordinate(iy_max));

    for_each(outside_coords.domain(), [&](IndexXY const ixy) {
        CoordX const coord_x(coordinate(ddc::select<IDimX>(ixy)));
        CoordY const coord_y(coordinate(ddc::select<IDimY>(ixy)));

        IndexX const ix(ddc::select<IDimX>(ixy));
        IndexY const iy(ddc::select<IDimY>(ixy));


        CoordX delta_coord_x;
        if (ix + 1 <= ix_max) {
            delta_coord_x = ddc::coordinate(ix + 1) - ddc::coordinate(ix);
        } else {
            delta_coord_x = CoordX(0.);
        }

        CoordY delta_coord_y;
        if (iy + 1 <= iy_max) {
            delta_coord_y = ddc::coordinate(iy + 1) - ddc::coordinate(iy);
        } else {
            delta_coord_y = CoordY(0.);
        }

        int sign_x = (coord_x >= 0) - (coord_x < 0);
        int sign_y = (coord_y >= 0) - (coord_y < 0);

        outside_coords(ixy) = CoordXY(
                coord_x + delta_coord_x + sign_x * (x_max - x_min) / 6.,
                coord_y + delta_coord_y + sign_y * (y_max - y_min) / 6.);

        double const outside_x = ddc::get<DimX>(outside_coords(ixy));
        double const outside_y = ddc::get<DimY>(outside_coords(ixy));

        // If the outside_coord is inside the domain, it is put outside the domain:
        if ((outside_x <= x_max) and (outside_x >= x_min) and (outside_y <= y_max)
            and (outside_y >= y_min)) {
            sign_x = (outside_x >= 0) - (outside_x < 0);
            sign_y = (outside_y >= 0) - (outside_y < 0);
            int const sign_xy = (sign_x * sign_y >= 0) - (sign_x * sign_y < 0);
            outside_coords(ixy) = outside_coords(ixy)
                                  + CoordXY(
                                          -sign_xy * sign_x * (x_max - x_min) * 1 / 3.,
                                          +sign_xy * sign_y * (y_max - y_min) * 1 / 3.);
        }
    });
}

/**
 * @brief Evaluate a function on B-splines and compare the obtained values with the exact values.
 *
 * We build another grid with points outside the domain and evaluate the function on B-splines on the grid.
 * If the point is outside the domain, we expect that the value of the function will be the same
 * than the value of the exact function at the edge of the domain.
 *
 * @param[in] grid
 * 			The polar grid used to define the discrete domain.
 * @param[in] exact_function
 * 			A function PolarExactFunction type.
 * @param[in] TOL
 * 			The tolerance error on the maximum absolute value of the difference between the exact value
 * 			and the computed value.
 */
template <class Function>
void Evaluate_on_outside_coord(IDomainRP const& grid, Function& exact_function, double const TOL)
{
    // Coordinates on the grid. --------------------------------------------------------------
    FieldRP<CoordRP> coords(grid);
    for_each(grid, [&](IndexRP const irp) { coords(irp) = ddc::coordinate(irp); });


    // Build the decomposition of the function on B-splines. ---------------------------------
    SplineRPBuilder const builder(grid);
    ddc::Chunk<double, BSDomainRP> function_coefs(builder.spline_domain());
    DFieldRP function_evaluated(grid);
    for_each(grid, [&](IndexRP const irp) {
        function_evaluated(irp) = exact_function(coords(irp));
    });
    builder(function_coefs, function_evaluated);


    // Build a "outside" grid to test the evaluator. ------------------------------------------
    IDomainRP outside_grid(grid);
    FieldRP<CoordRP> outside_coords(outside_grid);
    build_outside_grid(grid, outside_coords);


    // Evaluate the function on B-splines on the "outside" grid. ------------------------------
    auto r_domain = ddc::get_domain<IDimR>(outside_coords);
    CoordR const r_min(ddc::coordinate(r_domain.front()));
    CoordR const r_max(ddc::coordinate(r_domain.back()));
    ConstantExtrapolationBoundaryValue2D<BSplinesR, BSplinesP, DimR> boundary_condition_r_left(
            r_min);
    ConstantExtrapolationBoundaryValue2D<BSplinesR, BSplinesP, DimR> boundary_condition_r_right(
            r_max);

    SplineRPEvaluator spline_evaluator(
            boundary_condition_r_left,
            boundary_condition_r_right,
            g_null_boundary_2d<BSplinesR, BSplinesP>,
            g_null_boundary_2d<BSplinesR, BSplinesP>);

    spline_evaluator(function_evaluated.span_view(), outside_coords.span_cview(), function_coefs);


    // Compare the obtained values with the exact function. ----------------------------------
    check_constant_outside_domain(function_evaluated, exact_function, coords, outside_coords, TOL);
};



/** @brief Evaluate a function on B-splines and compare the obtained values with the exact values.
*
* We build another grid with points outside the domain and evaluate the function on B-splines on the grid.
* If the point is outside the domain, we expect that the value of the function will be the same
* than the value of the exact function at the edge of the domain.
*
* @param[in] grid
* 			The Cartesian grid used to define the discrete domain.
* @param[in] exact_function
* 			A function CartesianExactFunction type.
* @param[in] TOL
* 			The tolerance error on the maximum absolute value of the difference between the exact value
* 			and the computed value.
*/
template <class Function>
void Evaluate_on_outside_coord(IDomainXY const& grid, Function& exact_function, double const TOL)
{
    // Coordinates on the grid. --------------------------------------------------------------
    FieldXY<CoordXY> coords(grid);
    for_each(grid, [&](IndexXY const ixy) {
        CoordXY coord(ddc::coordinate(ixy));
        coords(ixy) = coord;
    });


    // Build the decomposition of the function on B-splines. ---------------------------------
    SplineXYBuilder const builder(grid);
    ddc::Chunk<double, BSDomainXY> function_coefs(builder.spline_domain());
    DFieldXY function_evaluated(grid);
    for_each(grid, [&](IndexXY const ixy) {
        CoordXY coord = coords(ixy);
        function_evaluated(ixy) = exact_function(coord);
    });
    builder(function_coefs, function_evaluated);


    // Build a "outside" grid to test the evaluator. ------------------------------------------
    IDomainXY outside_grid(grid);
    FieldXY<CoordXY> outside_coords(outside_grid);
    build_outside_grid(grid, outside_coords);


    // Evaluate the function on B-splines on the "outside" grid. ------------------------------
    auto x_domain = ddc::get_domain<IDimX>(outside_coords);
    CoordX const x_min(ddc::coordinate(x_domain.front()));
    CoordX const x_max(ddc::coordinate(x_domain.back()));
    auto y_domain = ddc::get_domain<IDimY>(outside_coords);
    CoordY const y_min(ddc::coordinate(y_domain.front()));
    CoordY const y_max(ddc::coordinate(y_domain.back()));
    ConstantExtrapolationBoundaryValue2D<BSplinesX, BSplinesY, DimX>
            boundary_condition_x_left(x_min, y_min, y_max);
    ConstantExtrapolationBoundaryValue2D<BSplinesX, BSplinesY, DimX>
            boundary_condition_x_right(x_max, y_min, y_max);
    ConstantExtrapolationBoundaryValue2D<BSplinesX, BSplinesY, DimY>
            boundary_condition_y_left(y_min, x_min, x_max);
    ConstantExtrapolationBoundaryValue2D<BSplinesX, BSplinesY, DimY>
            boundary_condition_y_right(y_max, x_min, x_max);

    SplineXYEvaluator spline_evaluator(
            boundary_condition_x_left,
            boundary_condition_x_right,
            boundary_condition_y_left,
            boundary_condition_y_right);

    spline_evaluator(function_evaluated.span_view(), outside_coords.span_cview(), function_coefs);


    // Compare the obtained values with the exact function. ----------------------------------
    check_constant_outside_domain(function_evaluated, exact_function, coords, outside_coords, TOL);
};



/**
 * @brief A class to define exact function in the polar domain.
 */
class PolarExactFunction
{
public:
    /**
     * @brief Instantiate a exact function in the polar domain.
     */
    virtual ~PolarExactFunction() = default;

    /**
	 * @brief Get the value of the function at the coordinate point.
	 *
	 * @param[in] coord
	 * 			The coordinate point in the polar domain.
	 *
	 * @return The value of the function at the coordinate point which is by default 0.
	 */
    virtual double operator()(CoordRP const& coord)
    {
        return 0.0;
    };
};

/**
 * @brief A class for functions of type @f$ (r, theta) \mapsto r^{d_r}\cos^{d_c}(\theta) \sin^{d_s}(\theta)@f$.
 */
class PolarExactFunction_r_theta_cos_sin : public PolarExactFunction
{
private:
    int const m_d_r;
    int const m_d_cos;
    int const m_d_sin;

public:
    /**
     * @brief Instantiate a type @f$ (r, theta) \mapsto r^{d_r}\cos^{d_c}(\theta) \sin^{d_s}(\theta)@f$.
     *
     * @param[in] d_r
     * 			The degree of the r.
     * @param[in] d_cos
     * 			The degree of the cosine.
     * @param[in] d_sin
     * 			The degree of the sine.
     */
    PolarExactFunction_r_theta_cos_sin(int const d_r, int const d_cos, int const d_sin)
        : m_d_r(d_r)
        , m_d_cos(d_cos)
        , m_d_sin(d_sin) {};
    ~PolarExactFunction_r_theta_cos_sin() {};

    /**
     * @brief Get the value of the function at the coordinate point.
     *
     * @param[in] coord
     * 			The coordinate point in the polar domain.
     *
     * @return The value of the function at the coordinate point.
     */
    double operator()(CoordRP const& coord) override
    {
        const double r = ddc::get<DimR>(coord);
        const double t = ddc::get<DimP>(coord);
        double val_r = 1.0;
        double val_cos = 1.0;
        double val_sin = 1.0;
        for (int i(0); i < m_d_r; i++)
            val_r *= r;
        for (int i(0); i < m_d_cos; i++)
            val_cos *= std::cos(t);
        for (int i(0); i < m_d_sin; i++)
            val_sin *= std::sin(t);
        return val_r * val_cos * val_sin;
    }
};


/**
 * @brief A class for functions of type @f$ (r, theta) \mapsto r*\cos^{d_c}(n\theta)@f$.
 */
class PolarExactFunction_r_cos_theta_pulsation : public PolarExactFunction
{
private:
    int const m_d_cos;
    int const m_n;

public:
    /**
     * @brief Instantiate a type @f$ (r, theta) \mapsto r*\cos^{d_c}(n\theta)@f$.
     *
     * @param[in] d_cos
     * 			The degree of the cosine.
     * @param[in] n
     * 			The number of oscillations.
     */
    PolarExactFunction_r_cos_theta_pulsation(int const d_cos, double const n)
        : m_d_cos(d_cos)
        , m_n(n) {};
    ~PolarExactFunction_r_cos_theta_pulsation() {};

    /**
     * @brief Get the value of the function at the coordinate point.
     *
     * @param[in] coord
     * 			The coordinate point in the polar domain.
     *
     * @return The value of the function at the coordinate point.
     */
    double operator()(CoordRP const& coord) override
    {
        const double r = ddc::get<DimR>(coord);
        const double t = ddc::get<DimP>(coord);
        double val_cos = 1.0;
        for (int i(0); i < m_d_cos; i++)
            val_cos *= std::cos(t * m_n);
        return r * val_cos;
    }
};



/**
 * @brief A class to define exact function in the Cartesian domain.
 */
class CartesianExactFunction
{
public:
    /**
     * @brief Instantiate a exact function in the Cartesian domain.
     */
    virtual ~CartesianExactFunction() = default;

    /**
	 * @brief Get the value of the function at the coordinate point.
	 *
	 * @param[in] coord
	 * 			The coordinate point in the Cartesian domain.
	 *
	 * @return The value of the function at the coordinate point which is by default 0.
	 */
    virtual double operator()(CoordXY const& coord)
    {
        return 0.0;
    };
};


/**
 * @brief A class to define exact function in the cartesian domain of type
 * @f$ (x,y) \mapsto \cos^{d_x}(n_x x)\sin^{d_y}(n_y y)@f$.
 */
class CartesianExactFunction_cos_x_sin_y : public CartesianExactFunction
{
private:
    int const m_d_x;
    double const m_n_x;
    int const m_d_y;
    double const m_n_y;

public:
    /**
     * @brief Instantiate a exact function in the cartesian domain of type
     * @f$ (x,y) \mapsto \cos^{d_x}(n_x x)\sin^{d_y}(n_y y)@f$.
     *
     * @param[in] d_x
     * 			The degree of the x-cosine.
     * @param[in] n_x
     * 			The number of oscillantion in the x direction.
     * @param[in] d_y
     * 			The degree of the y-sine.
     * @param[in] n_y
     * 			The number of oscillantion in the y direction.
     */
    CartesianExactFunction_cos_x_sin_y(
            int const d_x,
            double const n_x,
            int const d_y,
            double const n_y)
        : m_d_x(d_x)
        , m_n_x(n_x)
        , m_d_y(d_y)
        , m_n_y(n_y) {};
    ~CartesianExactFunction_cos_x_sin_y() {};

    /**
	 * @brief Get the value of the function at the coordinate point.
	 *
	 * @param[in] coord
	 * 			The coordinate point in the cartesian domain.
	 *
	 * @return The value of the function at the coordinate point.
	 */
    double operator()(CoordXY const& coord)
    {
        const double x = ddc::get<DimX>(coord);
        const double y = ddc::get<DimY>(coord);
        double val_cos = 1.0;
        double val_sin = 1.0;
        for (int i(0); i < m_d_x; i++)
            val_cos *= std::cos(x * m_n_x);
        for (int i(0); i < m_d_y; i++)
            val_sin *= std::sin(y * m_n_y);
        return val_cos * val_sin;
    };
};



} // end namespace



/**
 * @brief A class for the Google tests.
 */
class ConstantExtrapolationBCEvaluator2D
    : public testing::TestWithParam<std::tuple<std::size_t, std::size_t>>
{
};


namespace fs = std::filesystem;

TEST_P(ConstantExtrapolationBCEvaluator2D, PolarDomain)
{
    // INITIALISATION OF THE DISCRETE SPACE ==================================================
    // Parameters of the grid. ---------------------------------------------------------------
    auto const [Nr, Nt] = GetParam();

    // Grid creation (uniform grid). ----------------------------------------------------------
    CoordR const r_min(0.0);
    CoordR const r_max(1.0);
    IVectR const r_size(Nr);

    CoordP const p_min(0.0);
    CoordP const p_max(2.0 * M_PI);
    IVectP const p_size(Nt);

    IndexR const r_start(0); // with the singular point.
    IndexP const p_start(0);

    double const dr((r_max - r_min) / r_size);
    double const dp((p_max - p_min) / p_size);

    std::vector<CoordR> r_knots(r_size + 1);
    std::vector<CoordP> p_knots(p_size + 1);

    for (int i(0); i < r_size; ++i) {
        r_knots[i] = CoordR(r_min + i * dr);
    }
    r_knots[r_size] = CoordR(r_max);
    for (int i(0); i < p_size + 1; ++i) {
        p_knots[i] = CoordP(p_min + i * dp);
    }


    // Creating mesh & supports
    ddc::init_discrete_space<BSplinesR>(r_knots);
    ddc::init_discrete_space<BSplinesP>(p_knots);

    ddc::init_discrete_space<IDimR>(InterpPointsR::get_sampling());
    ddc::init_discrete_space<IDimP>(InterpPointsP::get_sampling());

    IDomainR interpolation_domain_R(InterpPointsR::get_domain());
    IDomainP interpolation_domain_P(InterpPointsP::get_domain());
    IDomainRP grid(interpolation_domain_R, interpolation_domain_P);



    // TESTS ON THE DISCRETE SPACE ===========================================================
    const int spline_degree = BSplinesR::degree();
    std::cout << "Test constant extrapolation as boundary conditions for evaluator 2D with "
                 "Bsplines "
              << "of degree " << spline_degree << " on r "
              << "and degree " << BSplinesP::degree() << " on theta"
              << " on a grid of " << Nr << " x " << Nt << "." << std::endl;

    // r-polynomials functions:
    for (int degree(std::max(spline_degree - 3, 0)); degree <= spline_degree + 3; degree++) {
        PolarExactFunction_r_theta_cos_sin exact_function_Rd(degree, 0, 0);
        Evaluate_on_outside_coord(grid, exact_function_Rd, 1e-15);
    }

    // r*cosine function:
    PolarExactFunction_r_theta_cos_sin exact_function_r_theta_cos(1, 1, 0);
    Evaluate_on_outside_coord(grid, exact_function_r_theta_cos, 1e-7);

    // r*sine function:
    PolarExactFunction_r_theta_cos_sin exact_function_r_theta_sin(1, 0, 1);
    Evaluate_on_outside_coord(grid, exact_function_r_theta_sin, 1e-7);

    // r *squared cosine function:
    PolarExactFunction_r_theta_cos_sin exact_function_r_theta_cos2(1, 2, 0);
    Evaluate_on_outside_coord(grid, exact_function_r_theta_cos2, 1e-6);

    // r *squared cosine function:
    PolarExactFunction_r_theta_cos_sin exact_function_r_theta_sin2(1, 0, 2);
    Evaluate_on_outside_coord(grid, exact_function_r_theta_sin2, 1e-6);

    // very oscillating r *squared cosine function:
    PolarExactFunction_r_cos_theta_pulsation exact_function_r_cos_pul(1, 10);
    Evaluate_on_outside_coord(grid, exact_function_r_cos_pul, 1e-3);
}



TEST_P(ConstantExtrapolationBCEvaluator2D, CartesianDomain)
{
    // INITIALISATION OF THE DISCRETE SPACE ==================================================
    // Parameters of the grid. ---------------------------------------------------------------
    auto const [Nx, Ny] = GetParam();

    // Grid creation (uniform grid). ----------------------------------------------------------
    CoordX const x_min(-1.0);
    CoordX const x_max(1.0);
    IVectX const x_size(Nx);

    CoordY const y_min(-1.0);
    CoordY const y_max(1.0);
    IVectY const y_size(Ny);

    IndexX const x_start(0); // with the singular point.
    IndexX const y_start(0);

    double const dx((x_max - x_min) / x_size);
    double const dy((y_max - y_min) / y_size);

    std::vector<CoordX> x_knots(x_size + 1);
    std::vector<CoordY> y_knots(y_size + 1);

    for (int i(0); i < x_size; ++i) {
        x_knots[i] = CoordX(x_min + i * dx);
    }
    x_knots[x_size] = CoordX(x_max);

    for (int i(0); i < y_size; ++i) {
        y_knots[i] = CoordY(y_min + i * dy);
    }
    y_knots[y_size] = CoordY(y_max);


    // Creating mesh & supports
    ddc::init_discrete_space<BSplinesX>(x_knots);
    ddc::init_discrete_space<BSplinesY>(y_knots);

    ddc::init_discrete_space<IDimX>(InterpPointsX::get_sampling());
    ddc::init_discrete_space<IDimY>(InterpPointsY::get_sampling());

    IDomainX interpolation_domain_X(InterpPointsX::get_domain());
    IDomainY interpolation_domain_Y(InterpPointsY::get_domain());
    IDomainXY grid(interpolation_domain_X, interpolation_domain_Y);



    // TESTS ON THE DISCRETE SPACE ===========================================================
    const int spline_degree_x = BSplinesX::degree();
    const int spline_degree_y = BSplinesY::degree();
    std::cout << "Test constant extrapolation as boundary conditions for evaluator 2D with "
                 "Bsplines "
              << "of degree " << spline_degree_x << " on x "
              << "and degree " << spline_degree_y << " on y"
              << " on a grid of " << Nx << " x " << Ny << "." << std::endl;

    // Cosine x function:
    CartesianExactFunction_cos_x_sin_y exact_function_cos_x(1, 1, 0, 0);
    Evaluate_on_outside_coord(grid, exact_function_cos_x, 1e-7);

    // Sine y function:
    CartesianExactFunction_cos_x_sin_y exact_function_sin_y(0, 0, 1, 1);
    Evaluate_on_outside_coord(grid, exact_function_sin_y, 1e-7);


    // Cosine x sine y function:
    CartesianExactFunction_cos_x_sin_y exact_function_cos_x_sin_y(1, 1, 1, 1);
    Evaluate_on_outside_coord(grid, exact_function_cos_x_sin_y, 1e-7);

    // Very oscillating cosine x sine y function:
    CartesianExactFunction_cos_x_sin_y exact_function_cos_x_sin_y_osc(1, 10, 1, 10);
    Evaluate_on_outside_coord(grid, exact_function_cos_x_sin_y_osc, 1e-4);
}



INSTANTIATE_TEST_SUITE_P(
        MyGroup,
        ConstantExtrapolationBCEvaluator2D,
        testing::Combine(testing::Values<std::size_t>(40), testing::Values<std::size_t>(80)));
