// SPDX-License-Identifier: MIT
#pragma once
#include <chrono>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <typeinfo>

#include <ddc/ddc.hpp>

#include "../../advection/r_theta_test_cases.hpp"

#include "bsl_advection_polar.hpp"
#include "geometry.hpp"
#include "l_norm_tools.hpp"
#include "math_tools.hpp"
#include "paraconfpp.hpp"
#include "params.yaml.hpp"
#include "polar_spline.hpp"
#include "polar_spline_evaluator.hpp"
#include "quadrature.hpp"
#include "spline_quadrature.hpp"
#include "trapezoid_quadrature.hpp"
#include "vector_field.hpp"
#include "vector_field_mem.hpp"
#include "volume_quadrature_nd.hpp"

namespace fs = std::filesystem;

std::string to_lower(std::string s)
{
    std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c) {
        if (c == ' ') {
            return '_';
        } else {
            return static_cast<char>(std::tolower(c));
        }
    });
    return s;
}

/**
 * Print the grid position to a file.
 *
 * @param[inout] out_file
 *      The stream to which the output is printed.
 * @param[in] coord_rtheta
 *      The coordinate to be printed.
 * @param[in] to_physical_mapping
 *      The mapping function from the logical domain to the physical domain.
 * @param[in] idx_range_theta
 *      An index range spanning the domain to which the poloidal coordinate should be restricted.
 */
template <class LogicalToPhysicalMapping>
void print_coordinate(
        std::ofstream& out_file,
        CoordRTheta coord_rtheta,
        LogicalToPhysicalMapping const& to_physical_mapping,
        IdxRangeTheta idx_range_theta)
{
    double const r = ddc::get<R>(coord_rtheta);
    double const th
            = ddcHelper::restrict_to_idx_range(ddc::select<Theta>(coord_rtheta), idx_range_theta);

    CoordXY coord_xy(to_physical_mapping(coord_rtheta));
    double const x = ddc::get<X>(coord_xy);
    double const y = ddc::get<Y>(coord_xy);

    out_file << std::setw(25) << r << std::setw(25) << th << std::setw(25) << x << std::setw(25)
             << y;
}

/**
 * @brief Save the characteristic feet in the logical domain
 * and the physical domain.
 *
 * @param[in] to_physical_mapping
 *      The mapping function from the logical domain to the physical
 *      domain.
 * @param[in] idx_range_rtheta
 *      The index range in the logical domain where the feet are defined.
 * @param[in] feet_coords_rtheta
 *      The characteristic feet in the logical domain.
 * @param[in] name
 *      The name of the file where the feet are saved.
 */
template <class LogicalToPhysicalMapping>
void output_feet(
        LogicalToPhysicalMapping const& to_physical_mapping,
        IdxRangeRTheta const& idx_range_rtheta,
        host_t<FieldRTheta<CoordRTheta>> const& feet_coords_rtheta,
        std::string const& name)
{
    std::ofstream file_feet(name, std::ofstream::out);
    file_feet << std::fixed << std::setprecision(16);

    IdxRangeR idx_range_r(idx_range_rtheta);
    IdxRangeTheta idx_range_theta(idx_range_rtheta);
    Idx<GridR> ir_start = idx_range_r.front();
    Idx<GridTheta> itheta_start = idx_range_theta.front();

    ddc::for_each(idx_range_rtheta, [&](IdxRTheta const irtheta) {
        IdxR ir(irtheta);
        IdxTheta itheta(irtheta);
        file_feet << std::setw(15) << (ir - ir_start).value() << std::setw(15)
                  << (itheta - itheta_start).value();
        print_coordinate(file_feet, ddc::coordinate(irtheta), to_physical_mapping, idx_range_theta);
        print_coordinate(
                file_feet,
                feet_coords_rtheta(irtheta),
                to_physical_mapping,
                idx_range_theta);
        file_feet << std::endl;
    });
    file_feet.close();
}


/**
 * @brief Save the advected function.
 *
 * @param[in] to_physical_mapping
 *      The mapping function from the logical domain to the physical
 *      domain.
 * @param[in] function
 *      The advected function.
 * @param[in] name
 *      The name of the file where the feet are saved.
 */
template <class LogicalToPhysicalMapping>
void saving_computed(
        LogicalToPhysicalMapping const& to_physical_mapping,
        host_t<DFieldRTheta> function,
        std::string const& name)
{
    IdxRangeRTheta const grid = get_idx_range(function);
    std::ofstream out_file(name, std::ofstream::out);
    out_file << std::fixed << std::setprecision(16);

    IdxRangeR idx_range_r(grid);
    IdxRangeTheta idx_range_theta(grid);
    Idx<GridR> ir_start = idx_range_r.front();
    Idx<GridTheta> itheta_start = idx_range_theta.front();

    ddc::for_each(grid, [&](IdxRTheta const irtheta) {
        IdxR const ir(irtheta);
        IdxTheta const itheta(irtheta);

        out_file << std::setw(15) << (ir - ir_start).value() << std::setw(15)
                 << (itheta - itheta_start).value();
        print_coordinate(out_file, ddc::coordinate(irtheta), to_physical_mapping, idx_range_theta);
        out_file << std::setw(25) << function(irtheta);
        out_file << std::endl;
    });
    out_file.close();
}

/**
 * @brief Get the exact characteristic feet of the simulation
 * at a given time.
 *
 * @param[in] idx_range_rtheta
 *      The logical domain where the characteristic feet are defined.
 * @param[in] to_physical_mapping
 *      The mapping function from the logical domain to the physical
 *      domain.
 * @param[in] advection_field
 *      The exact advection field in a AdvectionField object.
 * @param[in] time
 *      The time when we want to get the exact feet.
 *
 *
 * @return A FieldMem with the exact characteristic feet at the given time.
 */
template <class AdvectionField, class LogicalToPhysicalMapping, class PhysicalToLogicalMapping>
host_t<FieldMemRTheta<CoordRTheta>> compute_exact_feet_rtheta(
        IdxRangeRTheta const& idx_range_rtheta,
        LogicalToPhysicalMapping const& logical_to_physical_mapping,
        PhysicalToLogicalMapping const& physical_to_logical_mapping,
        AdvectionField const& advection_field,
        double const time)
{
    static_assert(!std::is_same_v<
                  LogicalToPhysicalMapping,
                  DiscreteToCartesian<X, Y, SplineRThetaEvaluatorConstBound_host>>);

    host_t<FieldMemRTheta<CoordRTheta>> feet_coords_rtheta(idx_range_rtheta);
    CoordXY const coord_xy_centre = CoordXY(logical_to_physical_mapping(CoordRTheta(0, 0)));
    ddc::for_each(idx_range_rtheta, [&](IdxRTheta const irtheta) {
        CoordRTheta const coord_rtheta = ddc::coordinate(irtheta);
        CoordXY const coord_xy
                = advection_field.exact_feet(logical_to_physical_mapping(coord_rtheta), time);

        CoordXY const coord_diff = coord_xy - coord_xy_centre;
        if (norm_inf(coord_diff) < 1e-15) {
            feet_coords_rtheta(irtheta) = CoordRTheta(0, 0);
        } else {
            feet_coords_rtheta(irtheta) = physical_to_logical_mapping(coord_xy);
        }
    });

    return feet_coords_rtheta;
}


/**
 * @brief Compute the L2 norm of the difference between 
 * the computed advected function and the exact
 * solution.
 *
 * @param[in] to_physical_mapping
 *      The mapping function from the logical domain to the physical
 *      domain.
 * @param[in] grid
 *      An index range spanning the logical domain where the function is defined.
 * @param[in] density_advected
 *      The computed function.
 * @param[in] function_to_be_advected
 *      The exact function.
 * @param[in] feet_coord
 *      The characteristic feet.
 *
 * @return The L2 norm of the difference between 
 * the computed function and the exact solution.
 */
template <class LogicalToPhysicalMapping, class Function>
double compute_difference_L2_norm(
        LogicalToPhysicalMapping const& to_physical_mapping,
        IdxRangeRTheta const& grid,
        host_t<DFieldRTheta> density_advected,
        Function& function_to_be_advected,
        host_t<FieldRTheta<CoordRTheta>> const& feet_coord)
{
    host_t<DFieldMemRTheta> exact_function(grid);
    host_t<DFieldMemRTheta> difference_function(grid);
    ddc::for_each(grid, [&](IdxRTheta const irtheta) {
        exact_function(irtheta) = function_to_be_advected(feet_coord(irtheta));
        difference_function(irtheta) = exact_function(irtheta) - density_advected(irtheta);
    });

    host_t<DFieldMemRTheta> const quadrature_coeffs = compute_coeffs_on_mapping(
            Kokkos::DefaultHostExecutionSpace(),
            to_physical_mapping,
            trapezoid_quadrature_coefficients<Kokkos::DefaultHostExecutionSpace>(grid));
    host_t<Quadrature<IdxRangeRTheta>> quadrature(get_const_field(quadrature_coeffs));

    double const normL2_exact_function
            = norm_L2(Kokkos::DefaultHostExecutionSpace(), quadrature, get_field(exact_function));
    double const normL2_difference_function = norm_L2(
            Kokkos::DefaultHostExecutionSpace(),
            quadrature,
            get_field(difference_function));

    return normL2_difference_function / normL2_exact_function;
}


/**
 * @brief Display the time difference.
 *
 * @param[in] start
 *      The time at the start of the period we want to display.
 * @param[in] end
 *      The time at the end of the period we want to display.
 * @param[in] title
 *      The title of the time difference for the display .
 */
void display_time(
        std::chrono::time_point<std::chrono::system_clock> const& start,
        std::chrono::time_point<std::chrono::system_clock> const& end,
        std::string const& title)
{
    double const time_h = std::chrono::duration_cast<std::chrono::hours>(end - start).count();
    double const time_min = std::chrono::duration_cast<std::chrono::minutes>(end - start).count();
    double const time_s = std::chrono::duration_cast<std::chrono::seconds>(end - start).count();
    double const time_ms
            = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << title << time_h << "h " << time_min - time_h * 60 << "min "
              << time_s - time_min * 60 << "s " << time_ms - time_s * 1000 << "ms " << std::endl;
}


/**
 * @brief Run an advection simulation.
 *
 * @param[in] to_physical_mapping
 *      The mapping function from the logical domain to the physical
 *      domain.
 * @param[in] to_logical_mapping
 *      The mapping function from the physical domain to the logical
 *      domain.
 * @param[in] analytical_to_pseudo_physical_mapping
 *      The mapping function from the logical domain to the pseudo-physical
 *      domain.
 * @param[in] analytical_to_physical_mapping
 *      The analytical version of the mapping.
 *      It can be different from the mapping if the mapping is discrete.
 * @param[in] grid
 *      An index range spanning the logical domain on which the advected function is defined.
 * @param[in] time_stepper
 *      The time integration method used to solve the characteristic
 *      equation.
 * @param[in] simulation
 *      The selected test cases.
 * @param[in] function_interpolator
 *      The B-splines interpolator used to interpolate the function at
 *      the characteristic feet.
 * @param[in] avection_builder
 *      The B-splines builder used to compute the B-splines coefficients
 *      of the advection field.
 * @param[in] advection_evaluator
 *      The B-splines evaluator used to evaluate the advection field.
 * @param[in] final_time
 *      The final time of the simulation.
 * @param[in] dt
 *      The time step.
 * @param[in] save_curves
 *      A boolean to select if the values of the function are saved in a text file
 *      for each time step. True: save in output folder; False: do not save.
 * @param[in] save_feet
 *      A boolean to select if the values of the characteristic for the last time step
 *      are saved in a text file. True: save in output folder; False: do not save.
 * @param[in] counter_function
 *      A integer referring to a test case for the name of the saved files.
 *
 * @see BslAdvection
 * @see ITimeStepper
 * @see Simulation
 */
template <
        class LogicalToPhysicalMappingHost,
        class LogicalToPhysicalMapping,
        class PhysicalToLogicalMapping,
        class LogicalToPseudoPhysicalMapping,
        class AnalyticalLogicalToPhysicalMapping,
        class PolarFootFinder,
        class AdvectionOperator,
        class Simulation>
void simulate(
        LogicalToPhysicalMappingHost const& to_physical_mapping_host,
        LogicalToPhysicalMapping const& to_physical_mapping,
        PhysicalToLogicalMapping const& to_logical_mapping,
        LogicalToPseudoPhysicalMapping const& analytical_to_pseudo_physical_mapping,
        AnalyticalLogicalToPhysicalMapping const& analytical_to_physical_mapping,
        IdxRangeRTheta const& grid,
        PolarFootFinder const& foot_finder,
        AdvectionOperator const& advection_operator,
        Simulation& simulation,
        double const final_time,
        double const dt,
        bool save_curves,
        bool save_feet,
        std::string const& output_folder)
{
    // TO CLOCK THE SIMULATION ------------------------------------------------------------------
    std::chrono::time_point<std::chrono::system_clock> start_simulation;
    std::chrono::time_point<std::chrono::system_clock> end_simulation;



    // PARAMETERS OF THE TEST -------------------------------------------------------------------
    int const iteration_number = int(final_time / dt);
    double const end_time = dt * iteration_number;


    host_t<DFieldMemRTheta> density_test(grid);

    host_t<DVectorFieldMemRTheta<X, Y>> advection_field_test_vec_host(grid);


    // START TEST -------------------------------------------------------------------------------
    start_simulation = std::chrono::system_clock::now();

    // Initialisation of the advected function:
    ddc::for_each(grid, [&](IdxRTheta const irtheta) {
        CoordRTheta coord = coordinate(irtheta);
        if (ddc::get<R>(coord) <= 1e-15) {
            ddc::get<Theta>(coord) = 0;
        }
        density_test(irtheta) = simulation.advected_function(coord);
    });



    // Definition of advection field:
    ddc::for_each(grid, [&](IdxRTheta const irtheta) {
        // Moving the coordinates in the physical domain:
        CoordXY const coord_xy = to_physical_mapping_host(ddc::coordinate(irtheta));

        // Define the advection field on the physical domain:
        ddcHelper::assign_vector_field_element(
                get_field(advection_field_test_vec_host),
                irtheta,
                simulation.advection_field(coord_xy, 0.));
    });

    auto density = ddc::
            create_mirror_view_and_copy(Kokkos::DefaultExecutionSpace(), get_field(density_test));

    auto advection_field_xy = ddcHelper::create_mirror_view_and_copy(
            Kokkos::DefaultExecutionSpace(),
            get_field(advection_field_test_vec_host));

    // SIMULATION -------------------------------------------------------------------------------
    // Advect "iteration_number" times:
    for (int i(0); i < iteration_number; ++i) {
        advection_operator(get_field(density), get_const_field(advection_field_xy), dt);

        // Save the advected function for each iteration:
        if (save_curves) {
            ddc::parallel_deepcopy(density_test, density);
            std::string const name = output_folder + "/after_" + std::to_string(i + 1) + ".txt";
            saving_computed(to_physical_mapping_host, get_field(density_test), name);
        }
    }
    ddc::parallel_deepcopy(density_test, density);



    // TREATMENT OF DATA ------------------------------------------------------------------------
    // Compute the exact characteristic feet:
    host_t<FieldMemRTheta<CoordRTheta>> feet_coords_rtheta_end_time(grid);
    host_t<FieldMemRTheta<CoordRTheta>> feet_coords_rtheta_dt(grid);
    feet_coords_rtheta_end_time = compute_exact_feet_rtheta(
            grid,
            analytical_to_physical_mapping,
            to_logical_mapping,
            simulation.advection_field,
            end_time);
    feet_coords_rtheta_dt = compute_exact_feet_rtheta(
            grid,
            analytical_to_physical_mapping,
            to_logical_mapping,
            simulation.advection_field,
            dt);


    // Compute the maximal absolute error on the space at the end of the simulation:
    double max_err = 0.;
    ddc::for_each(grid, [&](IdxRTheta const irtheta) {
        double const err
                = fabs(density_test(irtheta)
                       - simulation.advected_function(feet_coords_rtheta_end_time(irtheta)));
        max_err = max_err > err ? max_err : err;
    });


    std::cout << " - [" << iteration_number << " iterations]  Max absolute error : " << max_err;
    std::cout << std::endl;


    // Print the density difference [density conservation expected]:
    std::cout << "   ... "
              << "Relative L2 norm error: "
              << compute_difference_L2_norm(
                         to_physical_mapping_host,
                         grid,
                         get_field(density_test),
                         simulation.advected_function,
                         get_field(feet_coords_rtheta_end_time))
              << std::endl;



    end_simulation = std::chrono::system_clock::now();

    // Control simulation time
    display_time(start_simulation, end_simulation, "   ... Complete simulation time:        ");


    // SAVE DATA --------------------------------------------------------------------------------
    // Save the computed characteristic feet:
    if (save_feet) {
        FieldMemRTheta<CoordRTheta> feet_alloc(grid);
        FieldRTheta<CoordRTheta> feet = get_field(feet_alloc);
        ddc::parallel_for_each(
                grid,
                KOKKOS_LAMBDA(const IdxRTheta irtheta) {
                    feet(irtheta) = ddc::coordinate(irtheta);
                });
        auto advection_field_test_vec = ddcHelper::create_mirror_view_and_copy(
                Kokkos::DefaultExecutionSpace(),
                get_field(advection_field_test_vec_host));
        foot_finder(get_field(feet), get_const_field(advection_field_test_vec), dt);
        auto feet_host = ddc::create_mirror_view_and_copy(feet);
        std::string const name = output_folder + "/feet_computed.txt";
        output_feet(to_physical_mapping_host, grid, get_field(feet_host), name);
    }

    // Save the values of the exact function at the initial and final states:
    if (save_curves) {
        std::string const name_0 = output_folder + "/after_" + std::to_string(0) + ".txt";
        std::string const name_1
                = output_folder + "/after_" + std::to_string(iteration_number) + "_exact.txt";

        host_t<DFieldMemRTheta> initial_function(grid);
        host_t<DFieldMemRTheta> end_function(grid);
        ddc::for_each(grid, [&](const IdxRTheta irtheta) {
            initial_function(irtheta) = simulation.advected_function(ddc::coordinate(irtheta));

            // Exact final state
            end_function(irtheta)
                    = simulation.advected_function(feet_coords_rtheta_end_time(irtheta));
        });
        saving_computed(to_physical_mapping_host, get_field(initial_function), name_0);
        saving_computed(to_physical_mapping_host, get_field(end_function), name_1);
    }



    // Save the exact characteristic feet for a displacement on dt:
    if (save_feet) {
        std::string const name = output_folder + "/feet_exact.txt";
        output_feet(to_physical_mapping_host, grid, get_field(feet_coords_rtheta_dt), name);
    }


    std::cout << "--------------------------------------------------------------------------------"
              << std::endl;
}


/**
 * @brief Run three advection simulations for each test cases in the Simulation class.
 *
 * @param[in] to_physical_mapping
 *      The mapping function from the logical domain to the physical
 *      domain.
 * @param[in] to_logical_mapping
 *      The mapping function from the physical domain to the logical
 *      domain.
 * @param[in] analytical_to_pseudo_physical_mapping
 *      The mapping function from the logical domain to the pseudo-physical
 *      domain.
 * @param[in] analytical_to_physical_mapping
 *      The analytical version of the mapping.
 *      It can be different from the mapping if the mapping is discrete.
 * @param[in] grid
 *      An index range spanning the logical domain on which the advected function is defined.
 * @param[in] time_stepper
 *      The time integration method used to solve the characteristic
 *      equation.
 * @param[in] function_interpolator
 *      The B-splines interpolator used to interpolate the function at
 *      the characteristic feet.
 * @param[in] avection_builder
 *      The B-splines builder used to compute the B-splines coefficients
 *      of the advection field.
 * @param[in] advection_evaluator
 *      The B-splines evaluator used to evaluate the advection field.
 * @param[in] final_time
 *      The final time of the simulation.
 * @param[in] dt
 *      The time step.
 * @param[in] save_curves
 *      A boolean to select if the values of the function are saved in a text file
 *      for each time step. True: save in output folder; False: do not save.
 * @param[in] save_feet
 *      A boolean to select if the values of the characteristic for the last time step
 *      are saved in a text file. True: save in output folder; False: do not save.
 * @param[in] counter_function
 *      A integer referring to a test case for the name of the saved files.
 * @param[in] title
 *      The name of the test case written in the output of the console.
 *
 * @see BslAdvection
 * @see ITimeStepper
 * @see Simulation
 */
template <
        class LogicalToPhysicalMappingHost,
        class LogicalToPhysicalMapping,
        class PhysicalToLogicalMapping,
        class LogicalToPseudoPhysicalMapping,
        class AnalyticalLogicalToPhysicalMapping,
        class PolarFootFinder,
        class AdvectionOperator>
void run_simulations(
        LogicalToPhysicalMappingHost const& to_physical_mapping_host,
        LogicalToPhysicalMapping const& to_physical_mapping,
        PhysicalToLogicalMapping const& to_logical_mapping,
        LogicalToPseudoPhysicalMapping const& analytical_to_pseudo_physical_mapping,
        AnalyticalLogicalToPhysicalMapping const& analytical_to_physical_mapping,
        IdxRangeRTheta const& grid,
        PolarFootFinder const& foot_finder,
        AdvectionOperator const& advection_operator,
        double const final_time,
        double const dt,
        bool const& save_curves,
        bool const& save_feet,
        std::string const& output_stem,
        std::string const& title)
{
    IdxRangeR const r_idx_range = ddc::select<GridR>(grid);
    double const rmin = ddc::coordinate(r_idx_range.front());
    double const rmax = ddc::coordinate(r_idx_range.back());

    AdvectionSimulation simulation_translation
            = get_translation_simulation(to_physical_mapping, rmin, rmax);
    AdvectionSimulation simulation_rotation
            = get_rotation_simulation(to_physical_mapping, rmin, rmax);
    AdvectionSimulation simulation_decentred_rotation
            = get_decentred_rotation_simulation(to_physical_mapping);

    std::string const title_simu_translation = " TRANSLATION : ";
    std::string const title_simu_rotation = " ROTATION : ";
    std::string const title_simu_decentred_rotation = " DECENTRED ROTATION : ";

    std::string output_folder = output_stem + "Translation_output";
    if (save_curves or save_feet) {
        fs::create_directory(output_folder);
    }
    std::cout << title + title_simu_translation << std::endl;
    simulate(
            to_physical_mapping_host,
            to_physical_mapping,
            to_logical_mapping,
            analytical_to_pseudo_physical_mapping,
            analytical_to_physical_mapping,
            grid,
            foot_finder,
            advection_operator,
            simulation_translation,
            final_time,
            dt,
            save_curves,
            save_feet,
            output_folder);

    output_folder = output_stem + "Rotation_output";
    if (save_curves or save_feet) {
        fs::create_directory(output_folder);
    }
    std::cout << title + title_simu_rotation << std::endl;
    simulate(
            to_physical_mapping_host,
            to_physical_mapping,
            to_logical_mapping,
            analytical_to_pseudo_physical_mapping,
            analytical_to_physical_mapping,
            grid,
            foot_finder,
            advection_operator,
            simulation_rotation,
            final_time,
            dt,
            save_curves,
            save_feet,
            output_folder);

    output_folder = output_stem + "Decentred_rotation_output";
    if (save_curves or save_feet) {
        fs::create_directory(output_folder);
    }
    std::cout << title + title_simu_decentred_rotation << std::endl;
    simulate(
            to_physical_mapping_host,
            to_physical_mapping,
            to_logical_mapping,
            analytical_to_pseudo_physical_mapping,
            analytical_to_physical_mapping,
            grid,
            foot_finder,
            advection_operator,
            simulation_decentred_rotation,
            final_time,
            dt,
            save_curves,
            save_feet,
            output_folder);
}
