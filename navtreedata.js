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
    [ "Gyselalib++ simulations", "simulations.html", "simulations" ],
    [ "Gyselalib++ contents", "src.html", "src" ],
    [ "Gyselalib++ tests", "tests.html", "tests" ],
    [ "Selalib++", "vendor_sll.html", "vendor_sll" ],
    [ "API reference", "annotated.html", "annotated" ],
    [ "Files", "files.html", "files" ]
  ] ]
];

var NAVTREEINDEX =
[
"2d__spline__interpolator_2params_8yaml_8hpp_source.html",
"classCollisionsInter.html",
"classFFTPoissonSolver_3_01IdxRange_3_01GridPDEDim1D_8_8_8_01_4_00_01IdxRangeFull_00_01ExecSpace_00_01LayoutSpace_01_4.html#a6e31dd601d961ee711979be414fd2036",
"classMatrix.html#a0e72216b991d1f4af19982fe936067ab",
"classMultipatchType.html",
"classSplitRightHandSideSolver.html#a350971eed31753861f6e2917b9b6b968",
"geometryXVx_2poisson_2nullqnsolver_8hpp_source.html",
"src_geometryXYVxVy_poisson.html#src_geometryXYVxVy_poisson__Quasi-Neutrality_Solver",
"tests_geometryRTheta_2d_spline_interpolator.html#autotoc_md79"
];

var SYNCONMSG = 'click to disable panel synchronisation';
var SYNCOFFMSG = 'click to enable panel synchronisation';