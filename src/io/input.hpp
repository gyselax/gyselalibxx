// SPDX-License-Identifier: MIT
#pragma once
#include <paraconf.h>

#include "paraconfpp.hpp"

/**
 * Initialise a domain which will serve as an interpolation domain for splines.
 *
 * The domain is initialised using information from an input yaml file.
 * If the bsplines are uniform then the information to be read is:
 * - .SplineMesh.<mesh_identifier>_min
 * - .SplineMesh.<mesh_identifier>_max
 * - .SplineMesh.<mesh_identifier>_ncells
 *
 * If the bsplines are non-uniform then the information to be read is:
 * - .SplineMesh.<mesh_identifier>_MeshFile
 *
 * This string indicates the name of a file which contains the knots of the bspline.
 *
 * This information is used to initialise the bsplines. The interpolation domain
 * is then created using the specified method.
 */
template <class IDim, class BSplines, class InterpPointInitMethod>
inline ddc::DiscreteDomain<IDim> init_spline_dependent_domain(
        PC_tree_t const& conf_voicexx,
        std::string mesh_identifier)
{
    if constexpr (BSplines::is_uniform()) {
        using RDim = typename IDim::continuous_dimension_type;
        using Coord = ddc::Coordinate<RDim>;
        Coord min(PCpp_double(conf_voicexx, ".SplineMesh." + mesh_identifier + "_min"));
        Coord max(PCpp_double(conf_voicexx, ".SplineMesh." + mesh_identifier + "_max"));
        ddc::DiscreteVector<IDim> ncells(
                PCpp_int(conf_voicexx, ".SplineMesh." + mesh_identifier + "_ncells"));
        ddc::init_discrete_space<BSplines>(min, max, ncells);
    } else {
        throw "Non-uniform bspline initialisation is not yet handled";
    }
    ddc::init_discrete_space<IDim>(InterpPointInitMethod::template get_sampling<IDim>());
    ddc::DiscreteDomain<IDim> interpolation_domain(
            InterpPointInitMethod::template get_domain<IDim>());
    return interpolation_domain;
}

/**
 * Initialise a domain which will serve as an interpolation domain for splines.
 *
 * The domain is initialised using information from an input yaml file.
 * This function should be used for non-uniform bsplines, but it initialises
 * the break points uniformly. Such splines are referred to as pseudo-uniform
 * as the cells on which the polynomials are defined are uniform. However they
 * are not strictly uniform as multiple knots will be found at the same
 * position at the boundary.
 *
 * The information to be read from the file is:
 * - .SplineMesh.<mesh_identifier>_min
 * - .SplineMesh.<mesh_identifier>_max
 * - .SplineMesh.<mesh_identifier>_ncells
 *
 * The interpolation domain is then created using the specified method.
 */
template <class IDim, class BSplines, class InterpPointInitMethod>
inline ddc::DiscreteDomain<IDim> init_pseudo_uniform_spline_dependent_domain(
        PC_tree_t const& conf_voicexx,
        std::string mesh_identifier)
{
    static_assert(!BSplines::is_uniform());
    using RDim = typename IDim::continuous_dimension_type;
    using Coord = ddc::Coordinate<RDim>;

    Coord min(PCpp_double(conf_voicexx, ".SplineMesh." + mesh_identifier + "_min"));
    Coord max(PCpp_double(conf_voicexx, ".SplineMesh." + mesh_identifier + "_max"));
    ddc::DiscreteVector<IDim> ncells(
            PCpp_int(conf_voicexx, ".SplineMesh." + mesh_identifier + "_ncells"));

    std::vector<Coord> break_points(ncells + 1);

    double const delta(double(max - min) / ncells);

    break_points[0] = min;
    for (int i(1); i < ncells; ++i) {
        break_points[i] = min + i * delta;
    }
    break_points[ncells] = max;

    ddc::init_discrete_space<BSplines>(break_points);
    ddc::init_discrete_space<IDim>(InterpPointInitMethod::template get_sampling<IDim>());
    ddc::DiscreteDomain<IDim> interpolation_domain(
            InterpPointInitMethod::template get_domain<IDim>());
    return interpolation_domain;
}
