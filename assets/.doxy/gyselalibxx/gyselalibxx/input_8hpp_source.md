

# File input.hpp

[**File List**](files.md) **>** [**io**](dir_c184e51c84f2c3f0345bbc8a0d75d3e1.md) **>** [**input.hpp**](input_8hpp.md)

[Go to the documentation of this file](input_8hpp.md)


```C++
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

void parse_executable_arguments(
        PC_tree_t& conf_gyselalibxx,
        long int& iter_start,
        int argc,
        char** argv,
        char const* const params_yaml);

PC_tree_t parse_executable_arguments(int argc, char** argv, char const* const params_yaml);


template <class Grid1D, class BSplines, class InterpPointInitMethod>
inline IdxRange<Grid1D> init_spline_dependent_idx_range(
        PC_tree_t const& conf_gyselalibxx,
        std::string const& mesh_identifier)
{
    using Dim = typename Grid1D::continuous_dimension_type;
    using Coord1D = ddc::Coordinate<Dim>;

    std::vector<Coord1D> breakpoints;

    if constexpr (BSplines::is_uniform()) {
        // If uniform BSplines are used and interpolation points are calculated from them
        Coord1D min(PCpp_double(conf_gyselalibxx, ".SplineMesh." + mesh_identifier + "_min"));
        Coord1D max(PCpp_double(conf_gyselalibxx, ".SplineMesh." + mesh_identifier + "_max"));
        IdxStep<Grid1D> ncells(
                PCpp_int(conf_gyselalibxx, ".SplineMesh." + mesh_identifier + "_ncells"));
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

        if constexpr (BSplines::is_uniform()) {
            std::string grid_name = "grid_" + mesh_identifier;
            PDI_get_arrays("read_" + mesh_identifier, grid_name, mesh);
        } else {
            std::string breakpoints_name = "breakpoints_" + mesh_identifier;
            PDI_get_arrays(
                    "read_" + mesh_identifier,
                    breakpoints_name,
                    breakpoints,
                    mesh_identifier,
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

template <class Grid1D, class BSplines, class InterpPointInitMethod>
inline IdxRange<Grid1D> init_pseudo_uniform_spline_dependent_idx_range(
        PC_tree_t const& conf_gyselalibxx,
        std::string const& mesh_identifier)
{
    static_assert(!BSplines::is_uniform());
    using Dim = typename Grid1D::continuous_dimension_type;
    using Coord1D = Coord<Dim>;

    Coord1D min(PCpp_double(conf_gyselalibxx, ".SplineMesh." + mesh_identifier + "_min"));
    Coord1D max(PCpp_double(conf_gyselalibxx, ".SplineMesh." + mesh_identifier + "_max"));
    IdxStep<Grid1D> ncells(
            PCpp_int(conf_gyselalibxx, ".SplineMesh." + mesh_identifier + "_ncells"));

    std::vector<Coord1D> break_points = build_uniform_break_points(min, max, ncells);

    ddc::init_discrete_space<BSplines>(break_points);
    ddc::init_discrete_space<Grid1D>(InterpPointInitMethod::template get_sampling<Grid1D>());
    IdxRange<Grid1D> interpolation_idx_range(InterpPointInitMethod::template get_domain<Grid1D>());
    return interpolation_idx_range;
}
```


