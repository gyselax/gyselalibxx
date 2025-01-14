/*
 @licstart  The following is the entire license notice for the JavaScript code in this file.

 The MIT License (MIT)

 Copyright (C) 1997-2020 by Dimitri van Heesch

 Permission is hereby granted, free of charge, to any person obtaining a copy of this software
 and associated documentation files (the "Software"), to deal in the Software without restriction,
 including without limitation the rights to use, copy, modify, merge, publish, distribute,
 sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
 furnished to do so, subject to the following conditions:

 The above copyright notice and this permission notice shall be included in all copies or
 substantial portions of the Software.

 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING
 BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
 DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

 @licend  The above is the entire license notice for the JavaScript code in this file
*/
var NAVTREE =
[
  [ "Gyselalib++", "index.html", [
    [ "Set-up", "index.html#__Set-up", null ],
    [ "Compilation", "index.html#__Compilation", null ],
    [ "Execution", "index.html#__Execution", null ],
    [ "Dependencies", "index.html#__Dependencies", null ],
    [ "Pre-made build settings", "toolchains.html", [
      [ "General usage", "toolchains.html#toolchains__General_usage", null ],
      [ "Toolchains", "toolchains.html#toolchains__Toolchains", null ],
      [ "Environments", "toolchains.html#toolchains__Environments", null ],
      [ "Preparing environments", "toolchains.html#toolchains__Preparing_environments", null ],
      [ "General notes", "toolchains.html#toolchains__General_notes", null ]
    ] ],
    [ "Adding Documentation", "docs_Adding_docs.html", [
      [ "Building documentation locally", "docs_Adding_docs.html#docs_Adding_docs__Building_documentation_locally", null ],
      [ "Documentation describing code structures", "docs_Adding_docs.html#docs_Adding_docs__Documentation_describing_code_structures", null ],
      [ "Documentation describing general methods", "docs_Adding_docs.html#docs_Adding_docs__Documentation_describing_general_methods", null ],
      [ "Documenting functions", "docs_Adding_docs.html#docs_Adding_docs__Documenting_functions", null ],
      [ "Mathematical notation in documentation", "docs_Adding_docs.html#docs_Adding_docs__Mathematical_notation_in_documentation", null ]
    ] ],
    [ "Coding Standards", "docs_CODING_STANDARD.html", [
      [ "C++ Features", "docs_CODING_STANDARD.html#docs_CODING_STANDARD__Cxx_Features", null ],
      [ "Parameter passing", "docs_CODING_STANDARD.html#docs_CODING_STANDARD__Parameter_passing", null ],
      [ "Naming", "docs_CODING_STANDARD.html#docs_CODING_STANDARD__Naming", null ],
      [ "Style", "docs_CODING_STANDARD.html#docs_CODING_STANDARD__Style", null ],
      [ "Operators", "docs_CODING_STANDARD.html#docs_CODING_STANDARD__Operators", null ],
      [ "Code Organisation", "docs_CODING_STANDARD.html#docs_CODING_STANDARD__Code_Organisation", null ]
    ] ],
    [ "Common compilation problems", "docs_Common_compilation_problems.html", [
      [ "The closure type for a lambda cannot be used in the template argument type of a '__global__' function", "docs_Common_compilation_problems.html#docs_Common_compilation_problems__The_closure_type_for_a_lambda_cannot_be_used_in_the_template_argument_type_of_a____global____function", null ],
      [ "Implicit capture of 'this' in extended lambda expression", "docs_Common_compilation_problems.html#docs_Common_compilation_problems__Implicit_capture_of__this__in_extended_lambda_expression", null ],
      [ "Accessing allocated data", "docs_Common_compilation_problems.html#docs_Common_compilation_problems__Accessing_allocated_data", null ],
      [ "The enclosing parent function for an extended '__host__' '__device__' lambda cannot have private or protected access within its class", "docs_Common_compilation_problems.html#docs_Common_compilation_problems__The_enclosing_parent_function_for_an_extended____host_______device____lambda_cannot_have_private_or_protected_access_within_its_class", null ],
      [ "The enclosing parent function for an extended '__host__' '__device__' lambda must allow its address to be taken", "docs_Common_compilation_problems.html#docs_Common_compilation_problems__The_enclosing_parent_function_for_an_extended____host_______device____lambda_must_allow_its_address_to_be_taken", null ],
      [ "A nonstatic member reference must be relative to a specific object", "docs_Common_compilation_problems.html#docs_Common_compilation_problems__A_nonstatic_member_reference_must_be_relative_to_a_specific_object", null ],
      [ "X is not defined", "docs_Common_compilation_problems.html#docs_Common_compilation_problems__X_is_not_defined", null ]
    ] ],
    [ "Using DDC in Gyselalibxx", "docs_DDC_in_gyselalibxx.html", [
      [ "Contents", "docs_DDC_in_gyselalibxx.html#docs_DDC_in_gyselalibxx__Contents", null ],
      [ "Coordinates", "docs_DDC_in_gyselalibxx.html#docs_DDC_in_gyselalibxx__Coordinates", null ],
      [ "Indexing and associated concepts", "docs_DDC_in_gyselalibxx.html#docs_DDC_in_gyselalibxx__Indexing_and_associated_concepts", null ],
      [ "Grid", "docs_DDC_in_gyselalibxx.html#docs_DDC_in_gyselalibxx__Grid", null ],
      [ "Index", "docs_DDC_in_gyselalibxx.html#docs_DDC_in_gyselalibxx__Index", null ],
      [ "Index Step", "docs_DDC_in_gyselalibxx.html#docs_DDC_in_gyselalibxx__Index_Step", null ],
      [ "Index Range", "docs_DDC_in_gyselalibxx.html#docs_DDC_in_gyselalibxx__Index_Range", null ],
      [ "Data Storage", "docs_DDC_in_gyselalibxx.html#docs_DDC_in_gyselalibxx__Data_Storage", null ],
      [ "Example", "docs_DDC_in_gyselalibxx.html#docs_DDC_in_gyselalibxx__Example", null ],
      [ "Pitfalls", "docs_DDC_in_gyselalibxx.html#docs_DDC_in_gyselalibxx__Pitfalls", null ]
    ] ],
    [ "Using Git", "docs_Using_git.html", [
      [ "Branches", "docs_Using_git.html#docs_Using_git__Branches", null ],
      [ "Submodules", "docs_Using_git.html#docs_Using_git__Submodules", null ]
    ] ],
    [ "Developer's FAQ", "docs_developer_FAQ.html", [
      [ "What is the difference between Debug and Release mode?", "docs_developer_FAQ.html#docs_developer_FAQ__What_is_the_difference_between_Debug_and_Release_mode", null ],
      [ "Should I use abort, assert, or static_assert to raise an error?", "docs_developer_FAQ.html#docs_developer_FAQ__Should_I_use_abort_assert_or_static_assert_to_raise_an_error", null ]
    ] ],
    [ "Getting Started with Gyselalib++", "docs_getting_started.html", [
      [ "Understanding Functional Programming in Gyselalib++", "docs_getting_started.html#docs_getting_started__Understanding_Functional_Programming_in_Gyselalibxx", null ],
      [ "Navigating the Gyselalib++ Codebase", "docs_getting_started.html#docs_getting_started__Navigating_the_Gyselalibxx_Codebase", null ],
      [ "Recommended Steps for Getting Started", "docs_getting_started.html#docs_getting_started__Recommended_Steps_for_Getting_Started", null ]
    ] ],
    [ "Gyselalib++ simulations", "simulations.html", [
      [ "Simulations in (r, theta) geometry", "simulations_geometryRTheta.html", [
        [ "Diocotron instability", "simulations_geometryRTheta_diocotron.html", [
          [ "Studied problem", "simulations_geometryRTheta_diocotron.html#simulations_geometryRTheta_diocotron__Studied_problem", null ],
          [ "Test case - diocotron instability", "simulations_geometryRTheta_diocotron.html#simulations_geometryRTheta_diocotron__Test_case_-_diocotron_instability", null ],
          [ "References", "simulations_geometryRTheta_diocotron.html#simulations_geometryRTheta_diocotron__References", null ],
          [ "Contents", "simulations_geometryRTheta_diocotron.html#simulations_geometryRTheta_diocotron__Contents", null ]
        ] ],
        [ "Vortex merger", "simulations_geometryRTheta_vortex_merger.html", [
          [ "Studied problem", "simulations_geometryRTheta_vortex_merger.html#simulations_geometryRTheta_vortex_merger__Studied_problem", null ],
          [ "Test case - vortex merger", "simulations_geometryRTheta_vortex_merger.html#simulations_geometryRTheta_vortex_merger__Test_case_-_vortex_merger", null ],
          [ "References", "simulations_geometryRTheta_vortex_merger.html#simulations_geometryRTheta_vortex_merger__References", null ],
          [ "Useful references", "simulations_geometryRTheta_vortex_merger.html#simulations_geometryRTheta_vortex_merger__Useful_references", null ],
          [ "Contents", "simulations_geometryRTheta_vortex_merger.html#simulations_geometryRTheta_vortex_merger__Contents", null ]
        ] ]
      ] ],
      [ "Simulations in (x, vx) geometry", "simulations_geometryXVx.html", null ],
      [ "Simulations in (x, y) geometry", "simulations_geometryXY.html", [
        [ "Guiding center (X,Y) simulation", "simulations_geometryXY_guiding_center.html", [
          [ "Equations", "simulations_geometryXY_guiding_center.html#simulations_geometryXY_guiding_center__Equations", null ],
          [ "Simulation", "simulations_geometryXY_guiding_center.html#simulations_geometryXY_guiding_center__Simulation", null ],
          [ "Contents", "simulations_geometryXY_guiding_center.html#simulations_geometryXY_guiding_center__Contents", null ],
          [ "References", "simulations_geometryXY_guiding_center.html#simulations_geometryXY_guiding_center__References", null ]
        ] ]
      ] ]
    ] ],
    [ "Gyselalib++ contents", "src.html", [
      [ "Advection methods", "src_advection.html", [
        [ "Spatial advection", "src_advection.html#src_advection__Spatial_advection", null ],
        [ "Velocity advection", "src_advection.html#src_advection__Velocity_advection", null ],
        [ "1D advection with a given advection field", "src_advection.html#src_advection__1D_advection_with_a_given_advection_field", null ]
      ] ],
      [ "Collisions", "src_collisions.html", null ],
      [ "Data Storage Types", "src_data_types.html", [
        [ "VectorField", "src_data_types.html#src_data_types__VectorField", null ],
        [ "DerivField", "src_data_types.html#src_data_types__DerivField", null ]
      ] ],
      [ "Geometry (r, theta)", "src_geometryRTheta.html", [
        [ "Advection operator", "src_geometryRTheta_advection.html", [
          [ "Studied equation", "src_geometryRTheta_advection.html#src_geometryRTheta_advection__Studied_equation", null ],
          [ "Backward Semi-Lagrangian method", "src_geometryRTheta_advection.html#src_geometryRTheta_advection__Backward_Semi-Lagrangian_method", null ],
          [ "Time integration methods", "src_geometryRTheta_advection.html#src_geometryRTheta_advection__Time_integration_methods", null ],
          [ "Advection domain", "src_geometryRTheta_advection.html#src_geometryRTheta_advection__Advection_domain", null ],
          [ "Advection Field", "src_geometryRTheta_advection.html#src_geometryRTheta_advection__Advection_Field", null ],
          [ "Unit tests", "src_geometryRTheta_advection.html#autotoc_md51", null ],
          [ "References", "src_geometryRTheta_advection.html#autotoc_md52", null ],
          [ "Contents", "src_geometryRTheta_advection.html#autotoc_md53", null ]
        ] ],
        [ "Advection Field finder", "src_geometryRTheta_advection_field.html", [
          [ "Guiding center case", "src_geometryRTheta_advection_field.html#src_geometryRTheta_advection_field__Guiding_center_case", null ],
          [ "References", "src_geometryRTheta_advection_field.html#autotoc_md58", null ],
          [ "Contents", "src_geometryRTheta_advection_field.html#autotoc_md59", null ]
        ] ],
        [ "Geometry RTheta", "src_geometryRTheta_geometry.html", [
          [ "Shortcuts", "src_geometryRTheta_geometry.html#src_geometryRTheta_geometry__Shortcuts", null ]
        ] ],
        [ "Initialization", "src_geometryRTheta_initialization.html", [
          [ "Diocotron instability", "src_geometryRTheta_initialization.html#src_geometryRTheta_initialization__Diocotron_instability", null ],
          [ "Vortex merger", "src_geometryRTheta_initialization.html#src_geometryRTheta_initialization__Vortex_merger", null ],
          [ "Contents", "src_geometryRTheta_initialization.html#src_geometryRTheta_initialization__Contents", null ]
        ] ],
        [ "Spline interpolator in polar coordinates", "src_geometryRTheta_interpolation.html", null ],
        [ "Polar Poisson solver", "src_geometryRTheta_poisson.html", [
          [ "The Poisson-like equation", "src_geometryRTheta_poisson.html#src_geometryRTheta_poisson__The_Poisson-like_equation", null ],
          [ "Unit tests", "src_geometryRTheta_poisson.html#src_geometryRTheta_poisson__Unit_tests", null ],
          [ "References", "src_geometryRTheta_poisson.html#src_geometryRTheta_poisson__References", null ],
          [ "Contents", "src_geometryRTheta_poisson.html#src_geometryRTheta_poisson__Contents", null ]
        ] ],
        [ "Predictor-corrector methods", "src_geometryRTheta_time_solver.html", [
          [ "Predictor-corrector", "src_geometryRTheta_time_solver.html#src_geometryRTheta_time_solver__Predictor-corrector", null ],
          [ "Explicit predictor-corrector", "src_geometryRTheta_time_solver.html#src_geometryRTheta_time_solver__Explicit_predictor-corrector", null ],
          [ "Implicit predictor-corrector", "src_geometryRTheta_time_solver.html#src_geometryRTheta_time_solver__Implicit_predictor-corrector", null ],
          [ "References", "src_geometryRTheta_time_solver.html#src_geometryRTheta_time_solver__References", null ],
          [ "Contents", "src_geometryRTheta_time_solver.html#src_geometryRTheta_time_solver__Contents", null ]
        ] ]
      ] ],
      [ "Geometry (vpar, mu)", "src_geometryVparMu.html", [
        [ "GeometryVparMu :", "src_geometryVparMu_geometry.html", null ],
        [ "Initialization methods", "src_geometryVparMu_initialization.html", null ]
      ] ],
      [ "Geometry (x, v_x)", "src_geometryXVx.html", [
        [ "Boltzmann solver", "src_geometryXVx_boltzmann.html", null ],
        [ "Geometry X-Vx", "src_geometryXVx_geometry.html", null ],
        [ "Initialization methods", "src_geometryXVx_initialization.html", null ],
        [ "Quasi-Neutrality Solver", "src_geometryXVx_poisson.html", [
          [ "Charge Density", "src_geometryXVx_poisson.html#src_geometryXVx_poisson__Charge_Density", null ],
          [ "Poisson Solver", "src_geometryXVx_poisson.html#src_geometryXVx_poisson__Poisson_Solver", null ]
        ] ],
        [ "RHS", "src_geometryXVx_rhs.html", null ],
        [ "Time integration", "src_geometryXVx_time_integration.html", null ],
        [ "Utils", "src_geometryXVx_utils.html", null ]
      ] ],
      [ "Geometry (x, y)", "src_geometryXY.html", [
        [ "Geometry XY", "src_geometryXY_geometry.html", null ],
        [ "Initialization on (x,y) geometry", "src_geometryXY_initialization.html", [
          [ "Kelvin-Helmholtz instability test case", "src_geometryXY_initialization.html#src_geometryXY_initialization__Kelvin-Helmholtz_instability_test_case", null ]
        ] ],
        [ "Predictor-corrector methods", "src_geometryXY_time_integration.html", [
          [ "Predictor-corrector based on RK2", "src_geometryXY_time_integration.html#src_geometryXY_time_integration__Predictor-corrector_based_on_RK2", null ]
        ] ]
      ] ],
      [ "Geometry (x, y, v_x, v_y)", "src_geometryXYVxVy.html", [
        [ "Geometry X Y-Vx Vy", "src_geometryXYVxVy_geometry.html", null ],
        [ "Quasi-Neutrality Solver", "src_geometryXYVxVy_poisson.html", [
          [ "Charge Density", "src_geometryXYVxVy_poisson.html#src_geometryXYVxVy_poisson__Charge_Density", null ],
          [ "Quasi-Neutrality Solver", "src_geometryXYVxVy_poisson.html#src_geometryXYVxVy_poisson__Quasi-Neutrality_Solver", null ]
        ] ]
      ] ],
      [ "Interpolation Methods", "src_interpolation.html", [
        [ "Spline Interpolation", "src_interpolation.html#src_interpolation__Spline_Interpolation", null ],
        [ "Memory concerns", "src_interpolation.html#src_interpolation__Memory_concerns", null ],
        [ "Polar Splines", "src_interpolation_polar_splines.html", null ]
      ] ],
      [ "Functions used for input and output.", "src_io.html", null ],
      [ "Mappings", "src_mapping.html", null ],
      [ "Utility Functions", "src_math_tools.html", [
        [ "Derivative tools", "src_math_tools.html#src_math_tools__Derivative_tools", null ],
        [ "Utility tools", "src_math_tools.html#src_math_tools__Utility_tools", null ]
      ] ],
      [ "Matrix tools", "src_matrix_tools.html", null ],
      [ "Multipatch", "src_multipatch.html", [
        [ "Multipatch connectivity", "src_multipatch_connectivity.html", [
          [ "Patch", "src_multipatch_connectivity.html#src_multipatch_connectivity__Patch", null ],
          [ "Interfaces", "src_multipatch_connectivity.html#src_multipatch_connectivity__Interfaces", null ],
          [ "Patch locator", "src_multipatch_connectivity.html#src_multipatch_connectivity__Patch_locator", null ],
          [ "References", "src_multipatch_connectivity.html#src_multipatch_connectivity__References", null ],
          [ "Contents", "src_multipatch_connectivity.html#src_multipatch_connectivity__Contents", null ]
        ] ],
        [ "Data Types for Multipatch Geometry", "src_multipatch_data_types.html", [
          [ "MultipatchType", "src_multipatch_data_types.html#src_multipatch_data_types__MultipatchType", null ]
        ] ],
        [ "Spline on multipatch geometry", "src_multipatch_spline.html", [
          [ "Multipatch spline builder", "src_multipatch_spline.html#src_multipatch_spline__Multipatch_spline_builder", null ],
          [ "Multipatch spline evaluator", "src_multipatch_spline.html#src_multipatch_spline__Multipatch_spline_evaluator", null ],
          [ "Multipatch extrapolation rules", "src_multipatch_spline.html#src_multipatch_spline__Multipatch_extrapolation_rules", null ],
          [ "Contents", "src_multipatch_spline.html#src_multipatch_spline__Contents", null ]
        ] ],
        [ "Multipatch utilitary functions", "src_multipatch_utils.html", null ]
      ] ],
      [ "Parallelisation", "src_mpi_parallelisation.html", [
        [ "Layout", "src_mpi_parallelisation.html#src_mpi_parallelisation__Layout", null ],
        [ "TransposeOperator", "src_mpi_parallelisation.html#src_mpi_parallelisation__TransposeOperator", null ],
        [ "Alltoall Transpose Operator", "src_mpi_parallelisation.html#src_mpi_parallelisation__Alltoall_Transpose_Operator", null ]
      ] ],
      [ "PDE Solvers", "src_pde_solvers.html", [
        [ "Poisson's equation", "src_pde_solvers.html#src_pde_solvers__Poisson_s_equation", null ]
      ] ],
      [ "Quadrature Methods", "src_quadrature.html", null ],
      [ "SpeciesInfo (x, v_x)", "src_speciesinfo.html", null ],
      [ "Time Stepping Methods", "src_timestepper.html", null ],
      [ "Utility Functions", "src_utils.html", null ]
    ] ],
    [ "Gyselalib++ tests", "tests.html", [
      [ "Tests on the templated advection operators", "tests_advection.html", [
        [ "Contents", "tests_advection.html#tests_advection__Contents", null ]
      ] ],
      [ "Tests : Geometry (r, theta)", "tests_geometryRTheta.html", [
        [ "Tests on the 2D polar advection operator", "tests_geometryRTheta_advection_2d_rp.html", [
          [ "Tests on the 2D polar advection operator", "tests_geometryRTheta_advection_2d_rp.html#tests_geometryRTheta_advection_2d_rp__Tests_on_the_2D_polar_advection_operator", null ],
          [ "Python tests", "tests_geometryRTheta_advection_2d_rp.html#tests_geometryRTheta_advection_2d_rp__Python_tests", null ],
          [ "References", "tests_geometryRTheta_advection_2d_rp.html#tests_geometryRTheta_advection_2d_rp__References", null ],
          [ "Contents", "tests_geometryRTheta_advection_2d_rp.html#tests_geometryRTheta_advection_2d_rp__Contents", null ]
        ] ],
        [ "Tests on spline interpolator in polar coordinates", "tests_geometryRTheta_2d_spline_interpolator.html", null ],
        [ "Tests on the 2D polar poisson solver", "tests_geometryRTheta_polar_poisson.html", [
          [ "Polar Poisson solver", "tests_geometryRTheta_polar_poisson.html#tests_geometryRTheta_polar_poisson__Polar_Poisson_solver", null ],
          [ "References", "tests_geometryRTheta_polar_poisson.html#tests_geometryRTheta_polar_poisson__References", null ],
          [ "Contents", "tests_geometryRTheta_polar_poisson.html#tests_geometryRTheta_polar_poisson__Contents", null ]
        ] ]
      ] ],
      [ "Multipatch geometry tests", "tests_multipatch.html", [
        [ "Multipatch geometries", "tests_multipatch_geometries.html", null ]
      ] ]
    ] ],
    [ "Selalib++", "vendor_sll.html", null ],
    [ "API reference", "annotated.html", "annotated" ],
    [ "Files", "files.html", "files" ]
  ] ]
];

var NAVTREEINDEX =
[
"2d__spline__interpolator_2params_8yaml_8hpp_source.html",
"classCollisionsIntra.html#a164c6c249095a533900a12e2f1f48e25",
"classFluidMoments.html",
"classMatrix__Banded.html#a5b048b1936941d7188b996cd531cd630",
"classPolarBSplines_1_1Impl.html#a76ae9aee443a828ec5f62ddabbda37fb",
"classddcHelper_1_1NonUniformInterpolationPoints.html#afe391dfd5787f205a4096fef58e7beab",
"neumann__spline__quadrature_8hpp_source.html",
"structPatch_3_01grid1_00_01grid2_00_01bsplines__dim1_00_01bsplines__dim2_01_4.html#aed38b340e11a628505083d8ccace02a2"
];

var SYNCONMSG = 'click to disable panel synchronisation';
var SYNCOFFMSG = 'click to enable panel synchronisation';