// SPDX-License-Identifier: MIT
#pragma once
#include <sstream>

#include <ddc/ddc.hpp>

#include <paraconf.h>
#include <pdi.h>

#include "ddc_aliases.hpp"
#include "mesh_builder.hpp"
#include "non_uniform_interpolation_points.hpp"
#include "paraconfpp.hpp"
#include "pdi_helper.hpp"

/**
 * @brief Extract the paraconf configuration and the restart iteration from the executable arguments.
 *
 * @param[out] conf_voicexx The paraconf configuration describing the simulation.
 * @param[out] iter_start The index of the iteration from which the simulation should restart.
 * @param[in] argc The number of arguments passed to the executable.
 * @param[in] argv The arguments passed to the executable.
 * @param[in] params_yaml The default parameters for the yaml file.
 */
void parse_executable_arguments(
        PC_tree_t& conf_voicexx,
        long int& iter_start,
        int argc,
        char** argv,
        char const* const params_yaml);

/**
 * @brief Extract the paraconf configuration from the executable arguments.
 *
 * @param[in] argc The number of arguments passed to the executable.
 * @param[in] argv The arguments passed to the executable.
 * @param[in] params_yaml The default parameters for the yaml file.
 *
 * @returns The paraconf configuration describing the simulation.
 */
PC_tree_t parse_executable_arguments(int argc, char** argv, char const* const params_yaml);


/**
 * Initialise an index range which will serve as an interpolation index range for splines.
 *
 * The index range is initialised using information from an input yaml file.
 * If the B-splines are uniform then the information to be read is:
 * - .SplineMesh.<mesh_identifier>_min
 * - .SplineMesh.<mesh_identifier>_max
 * - .SplineMesh.<mesh_identifier>_ncells
 *
 * If the B-splines are non-uniform then the information to be read is:
 * - .SplineMesh.<mesh_identifier>_MeshFile
 *
 * This string indicates the name of a file which contains the knots of the bspline.
 *
 * This information is used to initialise the B-splines. The interpolation index range
 * is then created using the specified method.
 */
template <class Grid1D, class BSplines, class InterpPointInitMethod>
inline IdxRange<Grid1D> init_spline_dependent_idx_range(
        PC_tree_t const& conf_voicexx,
        std::string const& mesh_identifier)
{
    using Dim = typename Grid1D::continuous_dimension_type;
    using Coord1D = ddc::Coordinate<Dim>;

    std::vector<Coord1D> breakpoints;

    if constexpr (BSplines::is_uniform()) {
        // If uniform BSplines are used and interpolation points are calculated from them
        Coord1D min(PCpp_double(conf_voicexx, ".SplineMesh." + mesh_identifier + "_min"));
        Coord1D max(PCpp_double(conf_voicexx, ".SplineMesh." + mesh_identifier + "_max"));
        IdxStep<Grid1D> ncells(
                PCpp_int(conf_voicexx, ".SplineMesh." + mesh_identifier + "_ncells"));
        ddc::init_discrete_space<BSplines>(min, max, ncells);
    } else if constexpr (!ddcHelper::is_non_uniform_interpolation_points_v<InterpPointInitMethod>) {
        PDI_get_arrays("read_" + mesh_identifier, "breakpoints_" + mesh_identifier, breakpoints);
        ddc::init_discrete_space<BSplines>(breakpoints);
    }

    if constexpr (ddcHelper::is_non_uniform_interpolation_points_v<InterpPointInitMethod>) {
        // If uniform BSplines are used but the interpolation points are provided by the user
        // This may be the case if you want to test a new choice of interpolation points or
        // if you want to ensure that the interpolation points used match exactly the points
        // used to initialise values passed into the simulation.
        std::vector<Coord1D> mesh;
        std::string grid_name = "grid_" + mesh_identifier;

        if constexpr (BSplines::is_uniform()) {
            PDI_get_arrays("read_" + mesh_identifier, grid_name, mesh);
        } else {
            std::string breakpoints_name = "breakpoints_" + mesh_identifier;
            PDI_get_arrays(
                    "read_" + mesh_identifier,
                    breakpoints_name,
                    breakpoints,
                    grid_name,
                    mesh);
            ddc::init_discrete_space<BSplines>(breakpoints);
        }
        ddc::init_discrete_space<Grid1D>(
                InterpPointInitMethod::template get_sampling<Grid1D>(mesh));
    } else {
        ddc::init_discrete_space<Grid1D>(InterpPointInitMethod::template get_sampling<Grid1D>());
    }
    IdxRange<Grid1D> interpolation_idx_range(InterpPointInitMethod::template get_domain<Grid1D>());
    return interpolation_idx_range;
}

/**
 * Initialise an index range which will serve as an interpolation index range for splines.
 *
 * The index range is initialised using information from an input yaml file.
 * This function should be used for non-uniform B-splines, but it initialises
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
 * The interpolation index range is then created using the specified method.
 */
template <class Grid1D, class BSplines, class InterpPointInitMethod>
inline IdxRange<Grid1D> init_pseudo_uniform_spline_dependent_idx_range(
        PC_tree_t const& conf_voicexx,
        std::string const& mesh_identifier)
{
    static_assert(!BSplines::is_uniform());
    using Dim = typename Grid1D::continuous_dimension_type;
    using Coord1D = Coord<Dim>;

    Coord1D min(PCpp_double(conf_voicexx, ".SplineMesh." + mesh_identifier + "_min"));
    Coord1D max(PCpp_double(conf_voicexx, ".SplineMesh." + mesh_identifier + "_max"));
    IdxStep<Grid1D> ncells(PCpp_int(conf_voicexx, ".SplineMesh." + mesh_identifier + "_ncells"));

    std::vector<Coord1D> break_points = build_uniform_break_points(min, max, ncells);

    ddc::init_discrete_space<BSplines>(break_points);
    ddc::init_discrete_space<Grid1D>(InterpPointInitMethod::template get_sampling<Grid1D>());
    IdxRange<Grid1D> interpolation_idx_range(InterpPointInitMethod::template get_domain<Grid1D>());
    return interpolation_idx_range;
}
