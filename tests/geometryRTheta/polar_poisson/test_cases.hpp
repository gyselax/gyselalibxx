// SPDX-License-Identifier: MIT
#pragma once
#include <cmath>

#include <ddc/ddc.hpp>

#include <sll/math_tools.hpp>

#include "circular_to_cartesian.hpp"
#include "czarny_to_cartesian.hpp"
#include "ddc_aliases.hpp"
#include "geometry.hpp"
#include "inverse_jacobian_matrix.hpp"

/**
 * @brief Base class for the exact solutions of the Poisson equation.
 *
 * @see PolarSplineFEMPoissonLikeSolver
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
    static_assert(is_curvilinear_2d_mapping_v<CurvilinearToCartesian>);

private:
    using X = typename CurvilinearToCartesian::cartesian_tag_x;
    using Y = typename CurvilinearToCartesian::cartesian_tag_y;
    using R = typename CurvilinearToCartesian::curvilinear_tag_r;
    using Theta = typename CurvilinearToCartesian::curvilinear_tag_theta;

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
    explicit PoissonSolution(CurvilinearToCartesian const& coordinate_converter)
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
    virtual double operator()(Coord<R, Theta> const& coord) const = 0;
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
    explicit CurvilinearSolution(CurvilinearToCartesian const& coordinate_converter)
        : PoissonSolution<CurvilinearToCartesian>(coordinate_converter)
    {
    }

    double operator()(Coord<R, Theta> const& coord) const final
    {
        const double s = ddc::get<R>(coord);
        const double t = ddc::get<Theta>(coord);
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
private:
    InverseJacobianMatrix<CurvilinearToCartesian, Coord<R, Theta>> m_inverse_jacobian;

public:
    /**
     * @brief Instantiate a CartesianSolution.
     *
     * @param[in] coordinate_converter
     *      The mapping function which converts the logical (polar)
     *      coordinates into the physical (Cartesian) coordinates.
     */
    explicit CartesianSolution(CurvilinearToCartesian const& coordinate_converter)
        : PoissonSolution<CurvilinearToCartesian>(coordinate_converter)
        , m_inverse_jacobian(coordinate_converter)
    {
    }

    double operator()(Coord<R, Theta> const& coord) const final
    {
        const double s = ddc::get<R>(coord);
        const Coord<X, Y> cart_coord
                = PoissonSolution<CurvilinearToCartesian>::m_coordinate_converter(coord);
        const double x = ddc::get<X>(cart_coord);
        const double y = ddc::get<Y>(cart_coord);
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
    double derivative_x(Coord<R, Theta> const& coord) const
    {
        const Coord<X, Y> cart_coord
                = PoissonSolution<CurvilinearToCartesian>::m_coordinate_converter(coord);
        const double x = ddc::get<X>(cart_coord);
        const double y = ddc::get<Y>(cart_coord);

        const double r = ddc::get<R>(coord);
        const double dx_r = m_inverse_jacobian.inv_jacobian_11(coord);

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
    double derivative_y(Coord<R, Theta> const& coord) const
    {
        const Coord<X, Y> cart_coord
                = PoissonSolution<CurvilinearToCartesian>::m_coordinate_converter(coord);
        const double x = ddc::get<X>(cart_coord);
        const double y = ddc::get<Y>(cart_coord);

        const double r = ddc::get<R>(coord);
        const double dy_r = m_inverse_jacobian.inv_jacobian_12(coord);

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
    explicit ManufacturedPoissonTest(CurvilinearToCartesian const& coordinate_converter)
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
    double solution_at_pole(Coord<R, Theta> const& coord) const;
    /**
     * @brief Get the value of the RHS on the index range except the O-point.
     *
     * @param[in] coord
     *      The given polar coordinate.
     *
     * @return The value of the RHS on the index range except the O-point.
     */
    double non_singular_solution(Coord<R, Theta> const& coord) const;
    /**
     * @brief Get the value of the RHS at any point of the index range.
     *
     * @param[in] coord
     *      The given polar coordinate.
     *
     * @return The value of the RHS at any point of the index range.
     */
    double operator()(Coord<R, Theta> const& coord) const
    {
        if (ddc::get<R>(coord) == 0.0) {
            return solution_at_pole(coord);
        } else {
            return non_singular_solution(coord);
        }
    }
};
