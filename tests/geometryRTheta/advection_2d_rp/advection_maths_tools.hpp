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
    FILE* file_feet_ptr = fopen(name.c_str(), "w");
    ddc::for_each(rp_dom, [&](IndexRP const irp) {
        double const r = coordinate(ddc::select<IDimR>(irp));
        double const th = coordinate(ddc::select<IDimP>(irp));

        CoordXY coord_xy(mapping(CoordRP(r, th)));
        double const x = ddc::get<RDimX>(coord_xy);
        double const y = ddc::get<RDimY>(coord_xy);

        double const feet_r = ddc::get<RDimR>(feet_coords_rp(irp));
        double const feet_th = ddcHelper::
                restrict_to_domain(ddc::select<RDimP>(feet_coords_rp(irp)), IDomainP(rp_dom));

        CoordXY feet_xy(mapping(CoordRP(feet_r, feet_th)));
        double const feet_x = ddc::get<RDimX>(feet_xy);
        double const feet_y = ddc::get<RDimY>(feet_xy);


        fprintf(file_feet_ptr, "%10.0f", double(ddc::select<IDimR>(irp).uid()));
        fprintf(file_feet_ptr, "   %10.0f", double(ddc::select<IDimP>(irp).uid()));
        fprintf(file_feet_ptr, "   %10.16f", r);
        fprintf(file_feet_ptr, "   %10.16f", th);
        fprintf(file_feet_ptr, "   %10.16f", x);
        fprintf(file_feet_ptr, "   %10.16f", y);
        fprintf(file_feet_ptr, "   %10.16f", feet_r);
        fprintf(file_feet_ptr, "   %10.16f", feet_th);
        fprintf(file_feet_ptr, "   %10.16f", feet_x);
        fprintf(file_feet_ptr, "   %10.16f\n", feet_y);
    });
    fclose(file_feet_ptr);
};


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
    auto const grid = ddc::get_domain<IDimR, IDimP>(function);


    FILE* file_ptr = fopen(name.c_str(), "w");
    for_each(grid, [&](IndexRP const irp) {
        double const r = ddc::coordinate(ddc::select<IDimR>(irp));
        double const th = ddc::coordinate(ddc::select<IDimP>(irp));

        CoordXY coord_xy(mapping(CoordRP(r, th)));
        double const x = ddc::get<RDimX>(coord_xy);
        double const y = ddc::get<RDimY>(coord_xy);

        IndexR const ir(ddc::select<IDimR>(irp));
        IndexP const ip(ddc::select<IDimP>(irp));

        fprintf(file_ptr, "%10.0f", double(ir.uid()));
        fprintf(file_ptr, "   %10.0f", double(ip.uid()));
        fprintf(file_ptr, "   %10.0f", 0.);
        fprintf(file_ptr, "   %10.0f", 0.);
        fprintf(file_ptr, "   %10.16f", r);
        fprintf(file_ptr, "   %10.16f", th);
        fprintf(file_ptr, "   %10.16f", x);
        fprintf(file_ptr, "   %10.16f", y);
        fprintf(file_ptr, "   %10.16f \n", function(irp));
    });
    fclose(file_ptr);
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
    assert(typeid(mapping) != typeid(DiscreteToCartesian<RDimX, RDimY, SplineRPBuilder>));

    FieldRP<CoordRP> feet_coords_rp(rp_dom);
    CoordXY const coord_xy_center = CoordXY(mapping(CoordRP(0, 0)));

    ddc::for_each(rp_dom, [&](IndexRP const irp) {
        CoordRP const coord_rp = ddc::coordinate(irp);
        CoordXY coord_xy = mapping(coord_rp);
        double const x = ddc::get<RDimX>(coord_xy);
        double const y = ddc::get<RDimY>(coord_xy);

        double const feet_x = advection_field.exact_feet_x(x, y, time);
        double const feet_y = advection_field.exact_feet_y(x, y, time);
        coord_xy = CoordXY(feet_x, feet_y);

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
        double const feet_r = ddc::get<RDimR>(feet_coord(irp));
        double const feet_th = ddc::get<RDimP>(feet_coord(irp));
        exact_function(irp) = function_to_be_advected(feet_r, feet_th);
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
 * @param[in] if_save_curves
 *      A boolean to select if the values of the function are saved in a text file
 *      for each time step. True: save in output folder; False: do not save.
 * @param[in] if_save_feet
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
        bool const& if_save_curves,
        bool const& if_save_feet,
        int const& counter_function)
{
    if (if_save_curves or if_save_feet) {
        fs::create_directory("output");
    }
    if (if_save_curves) {
        fs::create_directory("output/curves_" + std::to_string(counter_function));
    }

    SplineFootFinder<TimeStepper, AdvectionDomain> const foot_finder(
            time_stepper,
            advection_domain,
            grid,
            advection_builder,
            advection_evaluator);

    BslAdvectionRP advection_operator(function_interpolator, foot_finder);
    auto function_to_be_advected_test = *(simulation.get_test_function());
    auto advection_field_test = *(simulation.get_advection_field());



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
        double const r = coordinate(ddc::select<IDimR>(irp));
        double th = coordinate(ddc::select<IDimP>(irp));
        if (r <= 1e-15) {
            th = 0;
        }
        allfdistribu_test(irp) = function_to_be_advected_test(r, th);
    });



    // Definition of advection field:
    for_each(grid, [&](IndexRP const irp) {
        // Moving the coordinates in the physical domain:
        CoordXY const coord_xy = mapping(ddc::coordinate(irp));
        double const x = ddc::get<RDimX>(coord_xy);
        double const y = ddc::get<RDimY>(coord_xy);

        // Define the advection field on the physical domain:
        ddcHelper::get<RDimX>(advection_field_test_vec)(irp)
                = advection_field_test.x_value(x, y, 0.);
        ddcHelper::get<RDimY>(advection_field_test_vec)(irp)
                = advection_field_test.y_value(x, y, 0.);
    });


    // SIMULATION -------------------------------------------------------------------------------
    // Advect "iteration_number" times:
    for (int i(0); i < iteration_number; ++i) {
        allfdistribu_advected_test
                = advection_operator(allfdistribu_test, advection_field_test_vec.span_cview(), dt);

        // Save the advected function for each iteration:
        if (if_save_curves) {
            std::string const name = "output/curves_" + std::to_string(counter_function) + "/after_"
                                     + std::to_string(i + 1) + ".txt";
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
        double const feet_r = ddc::get<RDimR>(feet_coords_rp_end_time(irp));
        double const feet_th = ddcHelper::restrict_to_domain(
                ddc::select<RDimP>(feet_coords_rp_end_time(irp)),
                ddc::get_domain<IDimP>(feet_coords_rp_end_time));

        double const err = fabs(
                allfdistribu_advected_test(irp) - function_to_be_advected_test(feet_r, feet_th));
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
        foot_finder(feet.span_view(), advection_field_test_vec, dt);
        std::string const name
                = "output/feet_computed_" + std::to_string(counter_function) + ".txt";
        save_feet(mapping, grid, feet.span_view(), name);
    }

    // Save the values of the exact function at the initial and final states:
    if (if_save_curves) {
        std::string const name_0 = "output/curves_" + std::to_string(counter_function) + "/after_"
                                   + std::to_string(0) + ".txt";
        std::string const name_1 = "output/curves_" + std::to_string(counter_function) + "/after_"
                                   + std::to_string(iteration_number) + "_exact.txt";

        DFieldRP initial_function(grid);
        DFieldRP end_function(grid);
        for_each(grid, [&](const IndexRP irp) {
            double const r = coordinate(ddc::select<IDimR>(irp));
            double const th = coordinate(ddc::select<IDimP>(irp));
            initial_function(irp) = function_to_be_advected_test(r, th);

            // Exact final state
            double const feet_r = ddc::get<RDimR>(feet_coords_rp_end_time(irp));
            double const feet_th = ddc::get<RDimP>(feet_coords_rp_end_time(irp));
            end_function(irp) = function_to_be_advected_test(feet_r, feet_th);
        });
        saving_computed(mapping, initial_function.span_view(), name_0);
        saving_computed(mapping, end_function.span_view(), name_1);
    }



    // Save the exact characteristic feet for a displacement on dt:
    if (if_save_feet) {
        std::string const name = "output/feet_exact_" + std::to_string(counter_function) + ".txt";
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
 * @param[in] if_save_curves
 *      A boolean to select if the values of the function are saved in a text file
 *      for each time step. True: save in output folder; False: do not save.
 * @param[in] if_save_feet
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
        bool const& if_save_curves,
        bool const& if_save_feet,
        int& counter_function,
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
            if_save_curves,
            if_save_feet,
            counter_function++);

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
            if_save_curves,
            if_save_feet,
            counter_function++);

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
            if_save_curves,
            if_save_feet,
            counter_function++);
}
