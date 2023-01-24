#include <cmath>

#include <ddc/ddc.hpp>

#include <geometry.hpp>
#include <mask_tanh.hpp>
#include <quadrature.hpp>
#include <trapezoid_quadrature.hpp>

//TODO: this should be directly handled by ddc::Discretization really,
//      in the meantime, we do it ourselves
template <class IDim>
constexpr std::enable_if_t<!IDim::continuous_dimension_type::PERIODIC, double>
total_interval_length(DiscreteDomain<IDim> const& dom)
{
    return std::fabs(rlength(dom));
}

//TODO: this should be directly handled by ddc::Discretization really,
//      in the meantime, we do it ourselves
template <class RDim>
constexpr std::enable_if_t<RDim::PERIODIC, double> total_interval_length(
        DiscreteDomain<UniformPointSampling<RDim>> const& dom)
{
    return std::fabs(rlength(dom) + step<UniformPointSampling<RDim>>());
}

/**
 * Returns a mask function defined with hyperbolic tangents
 *
 * Consider the domain [xmin, xmax], and {xleft, xright} the transition coordinates
 * defined using the extent parameter.
 *
 *  If type = 'normal' the mask equals one inside the [xleft, xright] interval
 *  and zero outside.
 *  If type = 'inverted' the mask equals zero inside the [xleft, xright] interval
 *  and one outside.
 *
 *  If m_normalized = true, the mask is normalized so that its integral equals one
 */

DFieldX mask_tanh(
        IDomainX const& gridx,
        double const extent,
        double const stiffness,
        MaskType const type,
        bool const normalized)
{
    DFieldX mask(gridx);

    IVectX const Nx(gridx.size());
    CoordX const x_min(coordinate(gridx.front()));
    double const Lx = total_interval_length(gridx);
    CoordX const x_left(x_min + Lx * extent);
    CoordX const x_right(x_min + Lx - Lx * extent);

    switch (type) {
    case MaskType::Normal:
        for_each(policies::parallel_host, gridx, [&](IndexX const ix) {
            CoordX const coordx = coordinate(ix);
            mask(ix) = 0.5
                       * (std::tanh((coordx - x_left) / stiffness)
                          - std::tanh((coordx - x_right) / stiffness));
        });
        break;

    case MaskType::Inverted:
        for_each(policies::parallel_host, gridx, [&](IndexX const ix) {
            CoordX const coordx = coordinate(ix);
            mask(ix) = 1
                       - 0.5
                                 * (std::tanh((coordx - x_left) / stiffness)
                                    - std::tanh((coordx - x_right) / stiffness));
        });
        break;
    }

    if (normalized) {
        Quadrature<IDimX> const integrate_x(trapezoid_quadrature_coefficients(gridx));
        double const coeff_norm = integrate_x(mask);
        for_each(policies::parallel_host, gridx, [&](IndexX const ix) {
            mask(ix) = mask(ix) / coeff_norm;
        });
    }

    return mask;
}
