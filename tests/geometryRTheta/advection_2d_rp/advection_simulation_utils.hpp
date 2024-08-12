#include <chrono>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <typeinfo>

#include <ddc/ddc.hpp>

#include <sll/math_tools.hpp>
#include <sll/polar_spline.hpp>
#include <sll/polar_spline_evaluator.hpp>

#include <stdio.h>

#include "bsl_advection_rp.hpp"
#include "compute_norms.hpp"
#include "directional_tag.hpp"
#include "geometry.hpp"
#include "paraconfpp.hpp"
#include "params.yaml.hpp"
#include "quadrature.hpp"
#include "spline_interpolator_2d_rp.hpp"
#include "spline_quadrature.hpp"
#include "test_cases.hpp"
#include "trapezoid_quadrature.hpp"
#include "utils_tools.hpp"
#include "vector_field.hpp"
#include "vector_field_mem.hpp"

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
 * @param[in] coord_rp
 *      The coordinate to be printed.
 * @param[in] mapping
 *      The mapping function from the logical index range to the physical index range.
 * @param[in] p_dom
 *      The index range to which the poloidal coordinate should be restricted.
 */
template <class Mapping>
void print_coordinate(
        std::ofstream& out_file,
        CoordRTheta coord_rp,
        Mapping const& mapping,
        IdxRangeTheta p_dom)
{
    double const r = ddc::get<R>(coord_rp);
    double const th = ddcHelper::restrict_to_idx_range(ddc::select<Theta>(coord_rp), p_dom);

    CoordXY coord_xy(mapping(coord_rp));
    double const x = ddc::get<X>(coord_xy);
    double const y = ddc::get<Y>(coord_xy);

    out_file << std::setw(25) << r << std::setw(25) << th << std::setw(25) << x << std::setw(25)
             << y;
}

/**
 * @brief Save the characteristic feet in the logical index range
 * and the physical index range.
 *
 * @param[in] mapping
 *      The mapping function from the logical index range to the physical
 *      index range.
 * @param[in] rp_index range
 *      The logical index range where the feet are defined.
 * @param[in] feet_coords_rp
 *      The characteristic feet in the logical index range.
 * @param[in] name
 *      The name of the file where the feet are saved.
 */
template <class Mapping>
void save_feet(
        Mapping const& mapping,
        IdxRangeRTheta const& rp_dom,
        FieldRTheta<CoordRTheta> const& feet_coords_rp,
        std::string const& name)
{
    std::ofstream file_feet(name, std::ofstream::out);
    file_feet << std::fixed << std::setprecision(16);
    ddc::for_each(rp_dom, [&](IdxRTheta const irp) {
        IdxRangeTheta p_dom = ddc::select<GridTheta>(rp_dom);

        file_feet << std::setw(15) << ddc::select<GridR>(irp).uid() << std::setw(15)
                  << ddc::select<GridTheta>(irp).uid();
        print_coordinate(file_feet, ddc::coordinate(irp), mapping, p_dom);
        print_coordinate(file_feet, feet_coords_rp(irp), mapping, p_dom);
        file_feet << std::endl;
    });
    file_feet.close();
}


/**
 * @brief Save the advected function.
 *
 * @param[in] mapping
 *      The mapping function from the logical index range to the physical
 *      index range.
 * @param[in] function
 *      The advected function.
 * @param[in] name
 *      The name of the file where the feet are saved.
 */
template <class Mapping>
void saving_computed(Mapping const& mapping, DFieldRTheta function, std::string const& name)
{
    IdxRangeRTheta const grid = get_idx_range(function);
    std::ofstream out_file(name, std::ofstream::out);
    out_file << std::fixed << std::setprecision(16);

    ddc::for_each(grid, [&](IdxRTheta const irp) {
        IdxRangeTheta p_dom = ddc::select<GridTheta>(grid);

        IdxR const ir(ddc::select<GridR>(irp));
        IdxTheta const ip(ddc::select<GridTheta>(irp));

        out_file << std::setw(15) << ir.uid() << std::setw(15) << ip.uid();
        print_coordinate(out_file, ddc::coordinate(irp), mapping, p_dom);
        out_file << std::setw(25) << function(irp);
        out_file << std::endl;
    });
    out_file.close();
}

/**
 * @brief Get the exact characteristic feet of the simulation
 * at a given time.
 *
 * @param[in] rp_dom
 *      The logical index range where the characteristic feet are defined.
 * @param[in] mapping
 *      The mapping function from the logical index range to the physical
 *      index range.
 * @param[in] advection_field
 *      The exact advection field in a AdvectionField object.
 * @param[in] time
 *      The time when we want to get the exact feet.
 *
 *
 * @return A FieldMem with the exact characteristic feet at the given time.
 */
template <class AdvectionField, class Mapping>
FieldMemRTheta<CoordRTheta> compute_exact_feet_rp(
        IdxRangeRTheta const& rp_dom,
        Mapping const& mapping,
        AdvectionField const& advection_field,
        double const time)
{
    static_assert(!std::is_same_v<
                  Mapping,
                  DiscreteToCartesian<X, Y, SplineRThetaBuilder, SplineRThetaEvaluatorConstBound>>);

    FieldMemRTheta<CoordRTheta> feet_coords_rp(rp_dom);
    CoordXY const coord_xy_center = CoordXY(mapping(CoordRTheta(0, 0)));

    ddc::for_each(rp_dom, [&](IdxRTheta const irp) {
        CoordRTheta const coord_rp = ddc::coordinate(irp);
        CoordXY const coord_xy = advection_field.exact_feet(mapping(coord_rp), time);

        CoordXY const coord_diff = coord_xy - coord_xy_center;
        if (norm_inf(coord_diff) < 1e-15) {
            feet_coords_rp(irp) = CoordRTheta(0, 0);
        } else {
            feet_coords_rp(irp) = mapping(coord_xy);
        }
    });

    return feet_coords_rp;
}


/**
 * @brief Compute the L2 norm of the difference between 
 * the computed advected function and the exact
 * solution.
 *
 * @param[in] mapping
 *      The mapping function from the logical index range to the physical
 *      index range.
 * @param[in] grid
 *      The logical index range where the function is defined.
 * @param[in] allfdistribu_advected
 *      The computed function.
 * @param[in] function_to_be_advected
 *      The exact function.
 * @param[in] feet_coord
 *      The characteristic feet.
 *
 * @return The L2 norm of the difference between 
 * the computed function and the exact solution.
 */
template <class Mapping, class Function>
double compute_difference_L2_norm(
        Mapping const& mapping,
        IdxRangeRTheta const& grid,
        DFieldRTheta allfdistribu_advected,
        Function& function_to_be_advected,
        FieldRTheta<CoordRTheta> const& feet_coord)
{
    DFieldMemRTheta exact_function(grid);
    DFieldMemRTheta difference_function(grid);
    ddc::for_each(grid, [&](IdxRTheta const irp) {
        exact_function(irp) = function_to_be_advected(feet_coord(irp));
        difference_function(irp) = exact_function(irp) - allfdistribu_advected(irp);
    });

    DFieldMemRTheta const quadrature_coeffs = compute_coeffs_on_mapping(
            mapping,
            trapezoid_quadrature_coefficients<Kokkos::DefaultHostExecutionSpace>(grid));
    host_t<Quadrature<IdxRangeRTheta>> quadrature(get_const_field(quadrature_coeffs));

    double const normL2_exact_function = compute_L2_norm(quadrature, get_field(exact_function));
    double const normL2_difference_function
            = compute_L2_norm(quadrature, get_field(difference_function));

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
 * @param[in] mapping
 *      The mapping function from the logical index range to the physical
 *      index range.
 * @param[in] analytical_mapping
 *      The analytical version of the mapping.
 *      It can be different from the mapping if the mapping is discrete.
 * @param[in] grid
 *      The logical index range on which the advected function is defined.
 * @param[in] time_stepper
 *      The time integration method used to solve the characteristic
 *      equation.
 * @param[in] advection_index range
 *      The AdvectionIdxRange type object defining the index range where we
 *      advect the characteristic feet.
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
 * @tparam Mapping
 *      A child class of CurvilinearToCartesian.
 * @tparam AnalyticalMapping
 *      A child class of AnalyticalInvertibleCurvilinearToCartesian.
 * @tparam TimeStepper
 *      A child class of ITimeStepper.
 * @tparam AdvectionIdxRange
 *      A child class of AdvectionIdxRange.
 * @tparam Simulation
 *      A child class of AdvectionSimulation.
 *
 * @see BslAdvection
 * @see ITimeStepper
 * @see AdvectionIdxRange
 * @see Simulation
 */
template <
        class Mapping,
        class AnalyticalMapping,
        class TimeStepper,
        class AdvectionIdxRange,
        class Simulation>
void simulate(
        Mapping const& mapping,
        AnalyticalMapping const& analytical_mapping,
        IdxRangeRTheta const& grid,
        TimeStepper const& time_stepper,
        AdvectionIdxRange& advection_idx_range,
        Simulation& simulation,
        PreallocatableSplineInterpolatorRTheta<ddc::NullExtrapolationRule> const&
                function_interpolator,
        SplineRThetaBuilder const& advection_builder,
        SplineRThetaEvaluatorConstBound& advection_evaluator,
        double const final_time,
        double const dt,
        bool if_save_curves,
        bool if_save_feet,
        std::string const& output_folder)
{
    SplineFootFinder<TimeStepper, AdvectionIdxRange> const
            foot_finder(time_stepper, advection_idx_range, advection_builder, advection_evaluator);

    BslAdvectionRTheta advection_operator(function_interpolator, foot_finder, mapping);
    auto function_to_be_advected_test = simulation.get_test_function();
    auto advection_field_test = simulation.get_advection_field();



    // TO CLOCK THE SIMULATION ------------------------------------------------------------------
    std::chrono::time_point<std::chrono::system_clock> start_simulation;
    std::chrono::time_point<std::chrono::system_clock> end_simulation;



    // PARAMETERS OF THE TEST -------------------------------------------------------------------
    int const iteration_number = int(final_time / dt);
    double const end_time = dt * iteration_number;


    DFieldMemRTheta allfdistribu_test(grid);
    DFieldRTheta allfdistribu_advected_test;

    DVectorFieldMemRTheta<X, Y> advection_field_test_vec(grid);


    // START TEST -------------------------------------------------------------------------------
    start_simulation = std::chrono::system_clock::now();

    // Initialization of the advected function:
    ddc::for_each(grid, [&](IdxRTheta const irp) {
        CoordRTheta coord = coordinate(irp);
        if (ddc::get<R>(coord) <= 1e-15) {
            ddc::get<Theta>(coord) = 0;
        }
        allfdistribu_test(irp) = function_to_be_advected_test(coord);
    });



    // Definition of advection field:
    ddc::for_each(grid, [&](IdxRTheta const irp) {
        // Moving the coordinates in the physical index range:
        CoordXY const coord_xy = mapping(ddc::coordinate(irp));
        CoordXY const advection_field = advection_field_test(coord_xy, 0.);

        // Define the advection field on the physical index range:
        ddcHelper::get<X>(advection_field_test_vec)(irp) = ddc::get<X>(advection_field);
        ddcHelper::get<Y>(advection_field_test_vec)(irp) = ddc::get<Y>(advection_field);
    });


    // SIMULATION -------------------------------------------------------------------------------
    // Advect "iteration_number" times:
    for (int i(0); i < iteration_number; ++i) {
        allfdistribu_advected_test = advection_operator(
                allfdistribu_test,
                get_const_field(advection_field_test_vec),
                dt);

        // Save the advected function for each iteration:
        if (if_save_curves) {
            std::string const name = output_folder + "/after_" + std::to_string(i + 1) + ".txt";
            saving_computed(mapping, get_field(allfdistribu_advected_test), name);
        }
    }



    // TREATMENT OF DATA ------------------------------------------------------------------------
    // Compute the exact characteristic feet:
    FieldMemRTheta<CoordRTheta> feet_coords_rp_end_time(grid);
    FieldMemRTheta<CoordRTheta> feet_coords_rp_dt(grid);
    feet_coords_rp_end_time
            = compute_exact_feet_rp(grid, analytical_mapping, advection_field_test, end_time);
    feet_coords_rp_dt = compute_exact_feet_rp(grid, analytical_mapping, advection_field_test, dt);


    // Compute the maximal absolute error on the space at the end of the simulation:
    double max_err = 0.;
    ddc::for_each(grid, [&](IdxRTheta const irp) {
        double const err
                = fabs(allfdistribu_advected_test(irp)
                       - function_to_be_advected_test(feet_coords_rp_end_time(irp)));
        max_err = max_err > err ? max_err : err;
    });


    std::cout << " - [" << iteration_number << " iterations]  Max absolute error : " << max_err;
    std::cout << std::endl;


    // Print the density difference [density conservation expected]:
    std::cout << "   ... "
              << "Relative L2 norm error: "
              << compute_difference_L2_norm(
                         mapping,
                         grid,
                         allfdistribu_advected_test,
                         function_to_be_advected_test,
                         get_field(feet_coords_rp_end_time))
              << std::endl;



    end_simulation = std::chrono::system_clock::now();

    // Control simulation time
    display_time(start_simulation, end_simulation, "   ... Complete simulation time:        ");


    // SAVE DATA --------------------------------------------------------------------------------
    // Save the computed characteristic feet:
    if (if_save_feet) {
        FieldMemRTheta<CoordRTheta> feet(grid);
        ddc::for_each(grid, [&](const IdxRTheta irp) { feet(irp) = ddc::coordinate(irp); });
        foot_finder(get_field(feet), advection_field_test_vec, dt);
        std::string const name = output_folder + "/feet_computed.txt";
        save_feet(mapping, grid, get_field(feet), name);
    }

    // Save the values of the exact function at the initial and final states:
    if (if_save_curves) {
        std::string const name_0 = output_folder + "/after_" + std::to_string(0) + ".txt";
        std::string const name_1
                = output_folder + "/after_" + std::to_string(iteration_number) + "_exact.txt";

        DFieldMemRTheta initial_function(grid);
        DFieldMemRTheta end_function(grid);
        ddc::for_each(grid, [&](const IdxRTheta irp) {
            initial_function(irp) = function_to_be_advected_test(ddc::coordinate(irp));

            // Exact final state
            end_function(irp) = function_to_be_advected_test(feet_coords_rp_end_time(irp));
        });
        saving_computed(mapping, get_field(initial_function), name_0);
        saving_computed(mapping, get_field(end_function), name_1);
    }



    // Save the exact characteristic feet for a displacement on dt:
    if (if_save_feet) {
        std::string const name = output_folder + "/feet_exact.txt";
        save_feet(mapping, grid, get_field(feet_coords_rp_dt), name);
    }


    std::cout << "--------------------------------------------------------------------------------"
              << std::endl;
}


/**
 * @brief Run three advection simulations for each test cases in the Simulation class.
 *
 * @param[in] mapping
 *      The mapping function from the logical index range to the physical
 *      index range.
 * @param[in] analytical_mapping
 *      The analytical version of the mapping.
 *      It can be different from the mapping if the mapping is discrete.
 * @param[in] grid
 *      The logical index range on which the advected function is defined.
 * @param[in] time_stepper
 *      The time integration method used to solve the characteristic
 *      equation.
 * @param[in] advection_index range
 *      The AdvectionIdxRange type object defining the index range where we
 *      advect the characteristic feet.
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
 * @see AdvectionIdxRange
 * @see Simulation
 */
template <class Mapping, class AnalyticalMapping, class TimeStepper, class AdvectionIdxRange>
void simulate_the_3_simulations(
        Mapping const& mapping,
        AnalyticalMapping const& analytical_mapping,
        IdxRangeRTheta const& grid,
        TimeStepper& time_stepper,
        AdvectionIdxRange& advection_idx_range,
        PreallocatableSplineInterpolatorRTheta<ddc::NullExtrapolationRule> const&
                function_interpolator,
        SplineRThetaBuilder const& advection_builder,
        SplineRThetaEvaluatorConstBound& advection_evaluator,
        double const final_time,
        double const dt,
        bool const& save_curves,
        bool const& save_feet,
        std::string const& output_stem,
        std::string const& title)
{
    auto const r_idx_range = ddc::select<GridR>(grid);
    double const rmin = ddc::coordinate(r_idx_range.front());
    double const rmax = ddc::coordinate(r_idx_range.back());

    TranslationSimulation simulation_1(mapping, rmin, rmax);
    RotationSimulation simulation_2(mapping, rmin, rmax);
    DecentredRotationSimulation simulation_3(mapping);

    std::string const title_simu_1 = " TRANSLATION : ";
    std::string const title_simu_2 = " ROTATION : ";
    std::string const title_simu_3 = " DECENTRED ROTATION : ";

    std::string output_folder = output_stem + "Translation_output";
    if (save_curves or save_feet) {
        fs::create_directory(output_folder);
    }
    std::cout << title + title_simu_1 << std::endl;
    simulate(
            mapping,
            analytical_mapping,
            grid,
            time_stepper,
            advection_idx_range,
            simulation_1,
            function_interpolator,
            advection_builder,
            advection_evaluator,
            final_time,
            dt,
            save_curves,
            save_feet,
            output_folder);

    output_folder = output_stem + "Rotation_output";
    if (save_curves or save_feet) {
        fs::create_directory(output_folder);
    }
    std::cout << title + title_simu_2 << std::endl;
    simulate(
            mapping,
            analytical_mapping,
            grid,
            time_stepper,
            advection_idx_range,
            simulation_2,
            function_interpolator,
            advection_builder,
            advection_evaluator,
            final_time,
            dt,
            save_curves,
            save_feet,
            output_folder);

    output_folder = output_stem + "Decentred_rotation_output";
    if (save_curves or save_feet) {
        fs::create_directory(output_folder);
    }
    std::cout << title + title_simu_3 << std::endl;
    simulate(
            mapping,
            analytical_mapping,
            grid,
            time_stepper,
            advection_idx_range,
            simulation_3,
            function_interpolator,
            advection_builder,
            advection_evaluator,
            final_time,
            dt,
            save_curves,
            save_feet,
            output_folder);
}
