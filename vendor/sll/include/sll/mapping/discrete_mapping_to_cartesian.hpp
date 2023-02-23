#pragma once

#include <ddc/ddc.hpp>

#include <sll/mapping/curvilinear2d_to_cartesian.hpp>
#include <sll/null_boundary_value.hpp>
#include <sll/spline_builder_2d.hpp>
#include <sll/spline_evaluator_2d.hpp>

template <class DimX, class DimY, class SplineBuilder>
class DiscreteToCartesian
    : public Curvilinear2DToCartesian<
              DimX,
              DimY,
              typename SplineBuilder::bsplines_type1::tag_type,
              typename SplineBuilder::bsplines_type2::tag_type>
{
public:
    using BSplineR = typename SplineBuilder::bsplines_type1;
    using BSplineP = typename SplineBuilder::bsplines_type2;
    using cartesian_tag_x = DimX;
    using cartesian_tag_y = DimY;
    using circular_tag_r = typename BSplineR::tag_type;
    using circular_tag_p = typename BSplineP::tag_type;

private:
    using interpolation_domain = typename SplineBuilder::interpolation_domain_type;
    using spline_domain = ddc::DiscreteDomain<BSplineR, BSplineP>;

private:
    ddc::Chunk<double, spline_domain> x_spline_representation;
    ddc::Chunk<double, spline_domain> y_spline_representation;
    SplineEvaluator2D<BSplineR, BSplineP> spline_evaluator;

public:
    DiscreteToCartesian(
            ddc::Chunk<double, spline_domain>&& curvilinear_to_x,
            ddc::Chunk<double, spline_domain>&& curvilinear_to_y)
        : x_spline_representation(std::move(curvilinear_to_x))
        , y_spline_representation(std::move(curvilinear_to_y))
        , spline_evaluator(
                  g_null_boundary_2d<BSplineR, BSplineP>,
                  g_null_boundary_2d<BSplineR, BSplineP>,
                  g_null_boundary_2d<BSplineR, BSplineP>,
                  g_null_boundary_2d<BSplineR, BSplineP>)
    {
    }

    ddc::Coordinate<DimX, DimY> operator()(
            ddc::Coordinate<circular_tag_r, circular_tag_p> const& coord) const final
    {
        const double x = spline_evaluator(coord, x_spline_representation);
        const double y = spline_evaluator(coord, y_spline_representation);
        return ddc::Coordinate<DimX, DimY>(x, y);
    }

    double jacobian_11(ddc::Coordinate<circular_tag_r, circular_tag_p> const& coord) const final
    {
        return spline_evaluator.deriv_dim_1(coord, x_spline_representation);
    }

    double jacobian_12(ddc::Coordinate<circular_tag_r, circular_tag_p> const& coord) const final
    {
        return spline_evaluator.deriv_dim_2(coord, x_spline_representation);
    }

    double jacobian_21(ddc::Coordinate<circular_tag_r, circular_tag_p> const& coord) const final
    {
        return spline_evaluator.deriv_dim_1(coord, y_spline_representation);
    }

    double jacobian_22(ddc::Coordinate<circular_tag_r, circular_tag_p> const& coord) const final
    {
        return spline_evaluator.deriv_dim_2(coord, y_spline_representation);
    }

    inline const ddc::Coordinate<DimX, DimY> control_point(
            ddc::DiscreteElement<BSplineR, BSplineP> const& el) const
    {
        return ddc::
                Coordinate<DimX, DimY>(x_spline_representation(el), y_spline_representation(el));
    }

    template <class PointSamplingR, class PointSamplingP>
    static inline ddc::Coordinate<circular_tag_r, circular_tag_p> get_coord(
            ddc::DiscreteElement<PointSamplingR, PointSamplingP> const& el)
    {
        return ddc::Coordinate<circular_tag_r, circular_tag_p>(
                ddc::coordinate(ddc::select<PointSamplingR>(el)),
                ddc::coordinate(ddc::select<PointSamplingP>(el)));
    }

    template <class Mapping, class Builder2D>
    static DiscreteToCartesian analytical_to_discrete(
            Mapping const& analytical_mapping,
            Builder2D const& builder)
    {
        using Domain = typename Builder2D::interpolation_domain_type;
        ddc::Chunk<double, spline_domain> curvilinear_to_x_spline(builder.spline_domain());
        ddc::Chunk<double, spline_domain> curvilinear_to_y_spline(builder.spline_domain());
        ddc::Chunk<double, Domain> curvilinear_to_x_vals(builder.interpolation_domain());
        ddc::Chunk<double, Domain> curvilinear_to_y_vals(builder.interpolation_domain());
        for_each(
                builder.interpolation_domain(),
                [&](typename Domain::discrete_element_type const& el) {
                    ddc::Coordinate<circular_tag_r, circular_tag_p> polar_coord(get_coord(el));
                    ddc::Coordinate<DimX, DimY> cart_coord = analytical_mapping(polar_coord);
                    curvilinear_to_x_vals(el) = ddc::select<DimX>(cart_coord);
                    curvilinear_to_y_vals(el) = ddc::select<DimY>(cart_coord);
                });
        builder(curvilinear_to_x_spline, curvilinear_to_x_vals);
        builder(curvilinear_to_y_spline, curvilinear_to_y_vals);
        return DiscreteToCartesian(
                std::move(curvilinear_to_x_spline),
                std::move(curvilinear_to_y_spline));
    }
};
