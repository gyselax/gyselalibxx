// SPDX-License-Identifier: MIT
/**
 * @file quadrature_coeffs_nd.hpp
 * File providing helper functions for defining multi-dimensional quadrature methods.
 */
#pragma once


/**
 * @brief Helper function which creates ND dimensions from N 1D quadrature coefficient functions.
 *
 * @param domain
 * 	The domain on which the coefficients will be defined.
 * @param func
 * 	The function which defines quadrature coefficients in the first dimension.
 * @param o_funcs
 * 	The functions which define quadrature coefficients in the subsequent dimensions.
 *
 * @returns The coefficients which define the quadrature method in ND.
 */
template <class IDim, class... ODims>
ddc::Chunk<double, ddc::DiscreteDomain<IDim, ODims...>> quadrature_coeffs_nd(
        ddc::DiscreteDomain<IDim, ODims...> const& domain,
        std::function<ddc::Chunk<double, ddc::DiscreteDomain<IDim>>(ddc::DiscreteDomain<IDim>)>
                func,
        std::function<ddc::Chunk<double, ddc::DiscreteDomain<ODims>>(
                ddc::DiscreteDomain<ODims>)>... o_funcs)
{
    // Get domains
    ddc::DiscreteDomain<IDim> const current_dim_domain = ddc::select<IDim>(domain);
    ddc::DiscreteDomain<ODims...> const other_dims_domain = ddc::select<ODims...>(domain);

    // Get coefficients in the first dimension
    ddc::Chunk<double, ddc::DiscreteDomain<IDim>> current_dim_coeffs = func(current_dim_domain);
    // Get coefficients over all other dimensions
    ddc::Chunk<double, ddc::DiscreteDomain<ODims...>> other_dim_coeffs
            = quadrature_coeffs_nd(other_dims_domain, o_funcs...);

    using DElemC = ddc::DiscreteElement<IDim>;
    using DElemO = ddc::DiscreteElement<ODims...>;

    // Calculate the combined coefficients over all dimensions
    ddc::Chunk<double, ddc::DiscreteDomain<IDim, ODims...>> coefficients(domain);

    ddc::for_each(current_dim_domain, [&](DElemC const idim) {
        ddc::for_each(other_dims_domain, [&](DElemO const odim) {
            coefficients(idim, odim) = current_dim_coeffs(idim) * other_dim_coeffs(odim);
        });
    });

    return coefficients;
}

/**
 * @brief Specialised 1D version of quadrature_coeffs_nd<IDim, ODims...>.
 *
 * @param domain
 * 	The domain on which the coefficients will be defined.
 * @param func
 * 	The function which defines quadrature coefficients.
 *
 * @returns The coefficients which define the quadrature method in 1D.
 */
#include <ddc/ddc.hpp>
template <class IDim>
ddc::Chunk<double, ddc::DiscreteDomain<IDim>> quadrature_coeffs_nd(
        ddc::DiscreteDomain<IDim> const& domain,
        std::function<ddc::Chunk<double, ddc::DiscreteDomain<IDim>>(ddc::DiscreteDomain<IDim>)>
                func)
{
    return func(domain);
}
