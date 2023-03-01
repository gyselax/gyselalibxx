#pragma once

#include <cmath>

#include <ddc/ddc.hpp>

#include <sll/mapping/circular_to_cartesian.hpp>
#include <sll/mapping/czarny_to_cartesian.hpp>
#include <sll/math_tools.hpp>

#include "poisson_geometry.hpp"

template <class CurvilinearToCartesian>
class PoissonSolution
{
public:
    using coordinate_converter_type = CurvilinearToCartesian;

private:
    using DimX = typename CurvilinearToCartesian::cartesian_tag_x;
    using DimY = typename CurvilinearToCartesian::cartesian_tag_y;
    using DimR = typename CurvilinearToCartesian::curvilinear_tag_r;
    using DimP = typename CurvilinearToCartesian::curvilinear_tag_p;

protected:
    CurvilinearToCartesian const& m_coordinate_converter;

public:
    PoissonSolution(CurvilinearToCartesian const& coordinate_converter)
        : m_coordinate_converter(coordinate_converter)
    {
    }
    virtual double operator()(Coordinate<DimR, DimP> const& coord) const = 0;
};

template <class CurvilinearToCartesian>
class CurvilinearSolution : public PoissonSolution<CurvilinearToCartesian>
{
public:
    CurvilinearSolution(CurvilinearToCartesian const& coordinate_converter)
        : PoissonSolution<CurvilinearToCartesian>(coordinate_converter)
    {
    }
    double operator()(Coordinate<DimR, DimP> const& coord) const final
    {
        const double s = get<DimR>(coord);
        const double t = get<DimP>(coord);
        return 1e-4 * ipow(s, 6) * ipow(s - 1, 6) / ipow(0.5, 12) * std::cos(11 * t);
    }
};

template <class CurvilinearToCartesian>
class CartesianSolution : public PoissonSolution<CurvilinearToCartesian>
{
public:
    CartesianSolution(CurvilinearToCartesian const& coordinate_converter)
        : PoissonSolution<CurvilinearToCartesian>(coordinate_converter)
    {
    }
    double operator()(Coordinate<DimR, DimP> const& coord) const final
    {
        const double s = get<DimR>(coord);
        const Coordinate<DimX, DimY> cart_coord
                = PoissonSolution<CurvilinearToCartesian>::m_coordinate_converter(coord);
        const double x = get<DimX>(cart_coord);
        const double y = get<DimY>(cart_coord);
        return 1e-4 * ipow(1 + s, 6) * ipow(1 - s, 6) * std::cos(2 * M_PI * x)
               * std::sin(2 * M_PI * y) / ipow(0.5, 12);
    }
};

template <class ChosenSolution>
class ManufacturedPoissonTest
{
public:
    using CurvilinearToCartesian = typename ChosenSolution::coordinate_converter_type;

private:
    CurvilinearToCartesian const& m_coordinate_converter;

public:
    ManufacturedPoissonTest(CurvilinearToCartesian const& coordinate_converter)
        : m_coordinate_converter(coordinate_converter)
    {
    }
    double solution_at_pole(Coordinate<DimR, DimP> const& coord) const;
    double non_singular_solution(Coordinate<DimR, DimP> const& coord) const;
    double operator()(Coordinate<DimR, DimP> const& coord) const
    {
        if (get<DimR>(coord) == 0.0) {
            return solution_at_pole(coord);
        } else {
            return non_singular_solution(coord);
        }
    }
};
