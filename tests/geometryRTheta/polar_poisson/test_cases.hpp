#pragma once

#include <cmath>

#include <ddc/ddc.hpp>

#include <sll/mapping/circular_to_cartesian.hpp>
#include <sll/mapping/czarny_to_cartesian.hpp>
#include <sll/math_tools.hpp>

#include "geometry.hpp"

/**
 * @brief Base class for the exact solutions of the Poisson equation.
 *
 * @see PolarSplineFEMPoissonSolver
 * @see VlasovPoissonSolver
 */
template <class CurvilinearToCartesian>
class PoissonSolution
{
public:
    /**
     * @brief Type the mapping function which converts the logical (polar)
     * coordinates into the physical (Cartesian) coordinates.
     */
    using coordinate_converter_type = CurvilinearToCartesian;

private:
    using RDimX = typename CurvilinearToCartesian::cartesian_tag_x;
    using RDimY = typename CurvilinearToCartesian::cartesian_tag_y;
    using RDimR = typename CurvilinearToCartesian::curvilinear_tag_r;
    using RDimP = typename CurvilinearToCartesian::curvilinear_tag_p;

protected:
    /**
     * @brief The mapping function which converts the logical (polar)
     * coordinates into the physical (Cartesian) coordinates.
     */
    CurvilinearToCartesian const& m_coordinate_converter;

public:
    virtual ~PoissonSolution() = default;

    /**
     * @brief Instantiate a PoissonSolution.
     *
     * @param[in] coordinate_converter
     *      The mapping function which converts the logical (polar)
     *      coordinates into the physical (Cartesian) coordinates.
     */
    PoissonSolution(CurvilinearToCartesian const& coordinate_converter)
        : m_coordinate_converter(coordinate_converter)
    {
    }
    /**
     * @brief Get the value of the function at a given coordinate.
     *
     * @param[in] coord
     *      The given polar coordinate.
     *
     * @return The value of the function at a given coordinate.
     */
    virtual double operator()(ddc::Coordinate<RDimR, RDimP> const& coord) const = 0;
};

/**
 * @brief Define a curvilinear solution of the Poisson equation.
 *
 * The solution is given by
 * * @f$ \phi(x, y) = C r(x,y)^6 (r(x,y) -1)^6 \cos(m\theta) @f$,
 *
 * with @f$C = 2^{12}1e-4 @f$ and @f$ m = 11 @f$.
 */
template <class CurvilinearToCartesian>
class CurvilinearSolution : public PoissonSolution<CurvilinearToCartesian>
{
public:
    /**
     * @brief Instantiate a CurvilinearSolution.
     *
     * @param[in] coordinate_converter
     *      The mapping function which converts the logical (polar)
     *      coordinates into the physical (Cartesian) coordinates.
     */
    CurvilinearSolution(CurvilinearToCartesian const& coordinate_converter)
        : PoissonSolution<CurvilinearToCartesian>(coordinate_converter)
    {
    }

    double operator()(ddc::Coordinate<RDimR, RDimP> const& coord) const final
    {
        const double s = ddc::get<RDimR>(coord);
        const double t = ddc::get<RDimP>(coord);
        return 1e-4 * ipow(s, 6) * ipow(s - 1, 6) / ipow(0.5, 12) * std::cos(11 * t);
    }
};

/**
 * @brief Define a Cartesian solution of the Poisson equation.
 *
 * The solution is given by
 * * @f$ \phi (x,y) = C (1+r(x,y))^6  (1 - r(x,y))^6 \cos(2\pi x) \sin(2\pi y) @f$,
 *
 * with @f$C = 2^{12}1e-4 @f$.
 *
 * Its x-derivative is
 *  * @f$ \partial_x \phi(x,y) = C \cdot 12 (r(x,y)^2 -1)^5 \cdot r(x,y) \cdot \partial_x r(x,y)\cdot \cos(2\pi x) \sin(2\pi y)
 *  - C  (r(x,y)^2 -1)^6 \cdot 2\pi  \sin(2\pi x) \sin(2\pi y), @f$
 *
 * and its y-derivative,
 *  * @f$ \partial_y \phi(x,y) = C \cdot 12 (r(x,y)^2 -1)^5 \cdot r(x,y) \cdot \partial_y r(x,y) \cdot \cos(2\pi x) \sin(2\pi y)
 *  + C  (r(x,y)^2 -1)^6 \cdot 2\pi  \cos(2\pi x) \cos(2\pi y). @f$
 */
template <class CurvilinearToCartesian>
class CartesianSolution : public PoissonSolution<CurvilinearToCartesian>
{
public:
    /**
     * @brief Instantiate a CartesianSolution.
     *
     * @param[in] coordinate_converter
     *      The mapping function which converts the logical (polar)
     *      coordinates into the physical (Cartesian) coordinates.
     */
    CartesianSolution(CurvilinearToCartesian const& coordinate_converter)
        : PoissonSolution<CurvilinearToCartesian>(coordinate_converter)
    {
    }

    double operator()(ddc::Coordinate<RDimR, RDimP> const& coord) const final
    {
        const double s = ddc::get<RDimR>(coord);
        const ddc::Coordinate<RDimX, RDimY> cart_coord
                = PoissonSolution<CurvilinearToCartesian>::m_coordinate_converter(coord);
        const double x = ddc::get<RDimX>(cart_coord);
        const double y = ddc::get<RDimY>(cart_coord);
        return 1e-4 * ipow(1 + s, 6) * ipow(1 - s, 6) * std::cos(2 * M_PI * x)
               * std::sin(2 * M_PI * y) / ipow(0.5, 12);
    }

    /**
     * @brief Get the value of the x-derivative of the function at a given coordinate.
     *
     * @param[in] coord
     *      The given polar coordinate.
     *
     * @return The value of the x-derivative of the function at a given coordinate.
     */
    double derivative_x(ddc::Coordinate<RDimR, RDimP> const& coord) const
    {
        const ddc::Coordinate<RDimX, RDimY> cart_coord
                = PoissonSolution<CurvilinearToCartesian>::m_coordinate_converter(coord);
        const double x = ddc::get<RDimX>(cart_coord);
        const double y = ddc::get<RDimY>(cart_coord);

        const double r = ddc::get<RDimR>(coord);
        const double dx_r
                = PoissonSolution<CurvilinearToCartesian>::m_coordinate_converter.inv_jacobian_11(
                        coord);

        const double C = 1e-4 * ipow(2., 12);
        return C * std::sin(2 * M_PI * y)
               * (12. * r * dx_r * ipow(r * r - 1, 5) * std::cos(2 * M_PI * x)
                  - 2 * M_PI * ipow(r * r - 1, 6) * std::sin(2 * M_PI * x));
    }

    /**
     * @brief Get the value of the y-derivative of the function at a given coordinate.
     *
     * @param[in] coord
     *      The given polar coordinate.
     *
     * @return The value of the y-derivative of the function at a given coordinate.
     */
    double derivative_y(ddc::Coordinate<RDimR, RDimP> const& coord) const
    {
        const ddc::Coordinate<RDimX, RDimY> cart_coord
                = PoissonSolution<CurvilinearToCartesian>::m_coordinate_converter(coord);
        const double x = ddc::get<RDimX>(cart_coord);
        const double y = ddc::get<RDimY>(cart_coord);

        const double r = ddc::get<RDimR>(coord);
        const double dy_r
                = PoissonSolution<CurvilinearToCartesian>::m_coordinate_converter.inv_jacobian_12(
                        coord);

        const double C = 1e-4 * ipow(2., 12);
        return C * std::cos(2 * M_PI * x)
               * (12. * r * dy_r * ipow(r * r - 1, 5) * std::sin(2 * M_PI * y)
                  + 2 * M_PI * ipow(r * r - 1, 6) * std::cos(2 * M_PI * y));
    }
};

/**
 * @brief Defining the corresponding RHS of the Poisson equation
 * for a given exact solution.
 */
template <class ChosenSolution>
class ManufacturedPoissonTest
{
public:
    /**
     * @brief Type the choosen solution of the Poisson equation.
     */
    using CurvilinearToCartesian = typename ChosenSolution::coordinate_converter_type;

private:
    CurvilinearToCartesian const& m_coordinate_converter;

public:
    /**
     * @brief Instantiate a ManufacturedPoissonTest.
     *
     * @param[in] coordinate_converter
     *      The mapping function which converts the logical (polar)
     *      coordinates into the physical (Cartesian) coordinates.
     */
    ManufacturedPoissonTest(CurvilinearToCartesian const& coordinate_converter)
        : m_coordinate_converter(coordinate_converter)
    {
    }
    /**
     * @brief Get the value of the RHS at the O-point.
     *
     * @param[in] coord
     *      The given polar coordinate.
     *
     * @return The value of the RHS at the O-point.
     */
    double solution_at_pole(ddc::Coordinate<RDimR, RDimP> const& coord) const;
    /**
     * @brief Get the value of the RHS on the domain except the O-point.
     *
     * @param[in] coord
     *      The given polar coordinate.
     *
     * @return The value of the RHS on the domain except the O-point.
     */
    double non_singular_solution(ddc::Coordinate<RDimR, RDimP> const& coord) const;
    /**
     * @brief Get the value of the RHS at any point of the domain.
     *
     * @param[in] coord
     *      The given polar coordinate.
     *
     * @return The value of the RHS at any point of the domain.
     */
    double operator()(ddc::Coordinate<RDimR, RDimP> const& coord) const
    {
        if (ddc::get<RDimR>(coord) == 0.0) {
            return solution_at_pole(coord);
        } else {
            return non_singular_solution(coord);
        }
    }
};
