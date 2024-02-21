#include <chrono>
#include <typeinfo>

#include <ddc/ddc.hpp>

#include <sll/constant_extrapolation_boundary_value.hpp>
#include <sll/math_tools.hpp>
#include <sll/null_boundary_value.hpp>
#include <sll/polar_spline.hpp>
#include <sll/polar_spline_evaluator.hpp>
#include <sll/spline_builder_2d.hpp>
#include <sll/spline_evaluator_2d.hpp>

#include "geometry.hpp"
#include "paraconfpp.hpp"
#include "params.yaml.hpp"
#include "spline_interpolator_2d_rp.hpp"
#include "utils_tools.hpp"


// ...
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iostream>

#include <stdio.h>
// ...


#include <sll/bsplines_non_uniform.hpp>
#include <sll/bsplines_uniform.hpp>
#include <sll/greville_interpolation_points.hpp>
#include <sll/spline_builder.hpp>
#include <sll/spline_evaluator.hpp>

#include <directional_tag.hpp>
#include <vector_field.hpp>
#include <vector_field_span.hpp>

#include "bsl_advection_rp.hpp"
#include "compute_norms.hpp"
#include "quadrature.hpp"
#include "spline_quadrature.hpp"
#include "test_cases.hpp"
#include "trapezoid_quadrature.hpp"


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
 *      The mapping function from the logical domain to the physical domain.
 * @param[in] p_dom
 *      The domain to which the poloidal coordinate should be restricted.
 */
template <class Mapping>
void print_coordinate(
        std::ofstream& out_file,
        CoordRP coord_rp,
        Mapping const& mapping,
        IDomainP p_dom)
{
    double const r = ddc::get<RDimR>(coord_rp);
    double const th = ddcHelper::restrict_to_domain(ddc::select<RDimP>(coord_rp), p_dom);

    CoordXY coord_xy(mapping(coord_rp));
    double const x = ddc::get<RDimX>(coord_xy);
    double const y = ddc::get<RDimY>(coord_xy);

    out_file << std::setw(25) << r << std::setw(25) << th << std::setw(25) << x << std::setw(25)
             << y;
}

/**
 * @brief Save the characteristic feet in the logical domain
 * and the physical domain.
 *
 * @param[in] mapping
 *      The mapping function from the logical domain to the physical
 *      domain.
 * @param[in] rp_domain
 *      The logical domain where the feet are defined.
 * @param[in] feet_coords_rp
 *      The characteristic feet in the logical domain.
 * @param[in] name
 *      The name of the file where the feet are saved.
 */
template <class Mapping>
void save_feet(
        Mapping const& mapping,
        IDomainRP const& rp_dom,
        SpanRP<CoordRP> const& feet_coords_rp,
        std::string const& name)
{
    std::ofstream file_feet(name, std::ofstream::out);
    file_feet << std::fixed << std::setprecision(16);
    ddc::for_each(rp_dom, [&](IndexRP const irp) {
        IDomainP p_dom = ddc::select<IDimP>(rp_dom);

        IndexR const ir(ddc::select<IDimR>(irp));
        IndexP const ip(ddc::select<IDimP>(irp));

        file_feet << std::setw(15) << ddc::select<IDimR>(irp).uid() << std::setw(15)
                  << ddc::select<IDimP>(irp).uid();
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
 *      The mapping function from the logical domain to the physical
 *      domain.
 * @param[in] function
 *      The advected function.
 * @param[in] name
 *      The name of the file where the feet are saved.
 */
template <class Mapping>
void saving_computed(Mapping const& mapping, DSpanRP function, std::string const& name)
{
    IDomainRP const grid = function.domain();
    std::ofstream out_file(name, std::ofstream::out);
    out_file << std::fixed << std::setprecision(16);

    for_each(grid, [&](IndexRP const irp) {
        IDomainP p_dom = ddc::select<IDimP>(grid);

        IndexR const ir(ddc::select<IDimR>(irp));
        IndexP const ip(ddc::select<IDimP>(irp));

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
 *      The logical domain where the characteristic feet are defined.
 * @param[in] mapping
 *      The mapping function from the logical domain to the physical
 *      domain.
 * @param[in] advection_field
 *      The exact advection field in a AdvectionField object.
 * @param[in] time
 *      The time when we want to get the exact feet.
 *
 *
 * @return A Chunk with the exact characteristic feet at the given time.
 */
template <class AdvectionField, class Mapping>
FieldRP<CoordRP> compute_exact_feet_rp(
        IDomainRP const& rp_dom,
        Mapping const& mapping,
        AdvectionField const& advection_field,
        double const time)
{
    static_assert(!std::is_same_v<Mapping, DiscreteToCartesian<RDimX, RDimY, SplineRPBuilder>>);

    FieldRP<CoordRP> feet_coords_rp(rp_dom);
    CoordXY const coord_xy_center = CoordXY(mapping(CoordRP(0, 0)));

    ddc::for_each(rp_dom, [&](IndexRP const irp) {
        CoordRP const coord_rp = ddc::coordinate(irp);
        CoordXY const coord_xy = advection_field.exact_feet(mapping(coord_rp), time);

        CoordXY const coord_diff = coord_xy - coord_xy_center;
        if (norm_inf(coord_diff) < 1e-15) {
            feet_coords_rp(irp) = CoordRP(0, 0);
        } else {
            feet_coords_rp(irp) = mapping(coord_xy);
        }
    });

    return feet_coords_rp;
}


/**
 * @brief Compute the difference between the L2 norm
 * of the computed advected function and the exact
 * solution.
 *
 * @param[in] mapping
 *      The mapping function from the logical domain to the physical
 *      domain.
 * @param[in] grid
 *      The logical domain where the function is defined.
 * @param[in] allfdistribu_advected
 *      The computed function.
 * @param[in] function_to_be_advected
 *      The exact function.
 * @param[in] feet_coord
 *      The characteristic feet.
 *
 * @return The difference between the L2 norm
 * of the computed function and the exact solution.
 */
template <class Mapping, class Function>
double compute_difference_L2_norm(
        Mapping const& mapping,
        IDomainRP const& grid,
        DSpanRP allfdistribu_advected,
        Function& function_to_be_advected,
        SpanRP<CoordRP> const& feet_coord)
{
    DFieldRP exact_function(grid);
    ddc::for_each(grid, [&](IndexRP const irp) {
        exact_function(irp) = function_to_be_advected(feet_coord(irp));
    });

    DFieldRP quadrature_coeffs
            = compute_coeffs_on_mapping(mapping, trapezoid_quadrature_coefficients(grid));
    Quadrature<IDimR, IDimP> quadrature(quadrature_coeffs);

    double const normL2_exact_function = compute_L2_norm(quadrature, exact_function.span_view());
    double const normL2_computed_function = compute_L2_norm(quadrature, allfdistribu_advected);

    return (normL2_exact_function - normL2_computed_function) / normL2_exact_function;
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
 *      The mapping function from the logical domain to the physical
 *      domain.
 * @param[in] analytical_mapping
 *      The analytical version of the mapping.
 *      It can be different from the mapping if the mapping is discrete.
 * @param[in] grid
 *      The logical domain on which the advected function is defined.
 * @param[in] time_stepper
 *      The time integration method used to solve the characteristic
 *      equation.
 * @param[in] advection_domain
 *      The AdvectionDomain type object defining the domain where we
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
 * @tparam AdvectionDomain
 *      A child class of AdvectionDomain.
 * @tparam Simulation
 *      A child class of AdvectionSimulation.
 *
 * @see BslAdvection
 * @see ITimeStepper
 * @see AdvectionDomain
 * @see Simulation
 */
template <
        class Mapping,
        class AnalyticalMapping,
        class TimeStepper,
        class AdvectionDomain,
        class Simulation>
void simulate(
        Mapping const& mapping,
        AnalyticalMapping const& analytical_mapping,
        IDomainRP const& grid,
        TimeStepper const& time_stepper,
        AdvectionDomain& advection_domain,
        Simulation& simulation,
        PreallocatableSplineInterpolatorRP const& function_interpolator,
        SplineRPBuilder const& advection_builder,
        SplineRPEvaluator& advection_evaluator,
        double const final_time,
        double const dt,
        bool if_save_curves,
        bool if_save_feet,
        std::string const& output_folder)
{
    SplineFootFinder<TimeStepper, AdvectionDomain> const foot_finder(
            time_stepper,
            advection_domain,
            grid,
            advection_builder,
            advection_evaluator);

    BslAdvectionRP advection_operator(function_interpolator, foot_finder);
    auto function_to_be_advected_test = simulation.get_test_function();
    auto advection_field_test = simulation.get_advection_field();



    // TO CLOCK THE SIMULATION ------------------------------------------------------------------
    std::chrono::time_point<std::chrono::system_clock> start_simulation;
    std::chrono::time_point<std::chrono::system_clock> end_simulation;



    // PARAMETERS OF THE TEST -------------------------------------------------------------------
    int const iteration_number = int(final_time / dt);
    double const end_time = dt * iteration_number;


    DFieldRP allfdistribu_test(grid);
    DSpanRP allfdistribu_advected_test;

    VectorDFieldRP<RDimX, RDimY> advection_field_test_vec(grid);


    // START TEST -------------------------------------------------------------------------------
    start_simulation = std::chrono::system_clock::now();

    // Initialization of the advected function:
    for_each(grid, [&](IndexRP const irp) {
        CoordRP coord = coordinate(irp);
        if (ddc::get<RDimR>(coord) <= 1e-15) {
            ddc::get<RDimP>(coord) = 0;
        }
        allfdistribu_test(irp) = function_to_be_advected_test(coord);
    });



    // Definition of advection field:
    for_each(grid, [&](IndexRP const irp) {
        // Moving the coordinates in the physical domain:
        CoordXY const coord_xy = mapping(ddc::coordinate(irp));
        CoordXY const advection_field = advection_field_test(coord_xy, 0.);

        // Define the advection field on the physical domain:
        ddcHelper::get<RDimX>(advection_field_test_vec)(irp) = ddc::get<RDimX>(advection_field);
        ddcHelper::get<RDimY>(advection_field_test_vec)(irp) = ddc::get<RDimY>(advection_field);
    });


    // SIMULATION -------------------------------------------------------------------------------
    // Advect "iteration_number" times:
    for (int i(0); i < iteration_number; ++i) {
        allfdistribu_advected_test
                = advection_operator(allfdistribu_test, advection_field_test_vec.span_cview(), dt);

        // Save the advected function for each iteration:
        if (if_save_curves) {
            std::string const name = output_folder + "/after_" + std::to_string(i + 1) + ".txt";
            saving_computed(mapping, allfdistribu_advected_test.span_view(), name);
        }
    }



    // TREATMENT OF DATA ------------------------------------------------------------------------
    // Compute the exact characteristic feet:
    FieldRP<CoordRP> feet_coords_rp_end_time(grid);
    FieldRP<CoordRP> feet_coords_rp_dt(grid);
    feet_coords_rp_end_time
            = compute_exact_feet_rp(grid, analytical_mapping, advection_field_test, end_time);
    feet_coords_rp_dt = compute_exact_feet_rp(grid, analytical_mapping, advection_field_test, dt);


    // Compute the maximal absolute error on the space at the end of the simulation:
    double max_err = 0.;
    for_each(grid, [&](IndexRP const irp) {
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
                         feet_coords_rp_end_time.span_view())
              << std::endl;



    end_simulation = std::chrono::system_clock::now();

    // Control simulation time
    display_time(start_simulation, end_simulation, "   ... Complete simulation time:        ");


    // SAVE DATA --------------------------------------------------------------------------------
    // Save the computed characteristic feet:
    if (if_save_feet) {
        FieldRP<CoordRP> feet(grid);
        for_each(grid, [&](const IndexRP irp) { feet(irp) = ddc::coordinate(irp); });
        foot_finder(feet.span_view(), advection_field_test_vec, dt);
        std::string const name = output_folder + "/feet_computed.txt";
        save_feet(mapping, grid, feet.span_view(), name);
    }

    // Save the values of the exact function at the initial and final states:
    if (if_save_curves) {
        std::string const name_0 = output_folder + "/after_" + std::to_string(0) + ".txt";
        std::string const name_1
                = output_folder + "/after_" + std::to_string(iteration_number) + "_exact.txt";

        DFieldRP initial_function(grid);
        DFieldRP end_function(grid);
        for_each(grid, [&](const IndexRP irp) {
            initial_function(irp) = function_to_be_advected_test(ddc::coordinate(irp));

            // Exact final state
            end_function(irp) = function_to_be_advected_test(feet_coords_rp_end_time(irp));
        });
        saving_computed(mapping, initial_function.span_view(), name_0);
        saving_computed(mapping, end_function.span_view(), name_1);
    }



    // Save the exact characteristic feet for a displacement on dt:
    if (if_save_feet) {
        std::string const name = output_folder + "/feet_exact.txt";
        save_feet(mapping, grid, feet_coords_rp_dt.span_view(), name);
    }


    std::cout << "--------------------------------------------------------------------------------"
              << std::endl;
}


/**
 * @brief Run three advection simulations for each test cases in the Simulation class.
 *
 * @param[in] mapping
 *      The mapping function from the logical domain to the physical
 *      domain.
 * @param[in] analytical_mapping
 *      The analytical version of the mapping.
 *      It can be different from the mapping if the mapping is discrete.
 * @param[in] grid
 *      The logical domain on which the advected function is defined.
 * @param[in] time_stepper
 *      The time integration method used to solve the characteristic
 *      equation.
 * @param[in] advection_domain
 *      The AdvectionDomain type object defining the domain where we
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
 * @see AdvectionDomain
 * @see Simulation
 */
template <class Mapping, class AnalyticalMapping, class TimeStepper, class AdvectionDomain>
void simulate_the_3_simulations(
        Mapping const& mapping,
        AnalyticalMapping const& analytical_mapping,
        IDomainRP const& grid,
        TimeStepper& time_stepper,
        AdvectionDomain& advection_domain,
        PreallocatableSplineInterpolatorRP const& function_interpolator,
        SplineRPBuilder const& advection_builder,
        SplineRPEvaluator& advection_evaluator,
        double const final_time,
        double const dt,
        bool const& save_curves,
        bool const& save_feet,
        std::string const& output_stem,
        std::string const& title)
{
    auto const r_domain = ddc::select<IDimR>(grid);
    double const rmin = ddc::coordinate(r_domain.front());
    double const rmax = ddc::coordinate(r_domain.back());

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
            advection_domain,
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
            advection_domain,
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
            advection_domain,
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
