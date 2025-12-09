# SingleInterfaceDerivativeCalculator test =======================================================
function (multipatch_single_interface_derivatives_executable MESH_TYPE)
    set(test_name_gtest "multipatch_single_interface_derivatives_${MESH_TYPE}_test")
    add_executable("${test_name_gtest}"
        single_interface_derivatives_calculator_test.cpp
        ../../main.cpp
    )
    target_link_libraries("${test_name_gtest}"
    PUBLIC
        DDC::core
        GTest::gtest
        GTest::gmock

        gslx::multipatch_geometries
        gslx::multipatch_spline
        gslx::multipatch_interface_derivatives
    )
    target_compile_definitions("${test_name_gtest}" PUBLIC -D${MESH_TYPE})
    gtest_discover_tests("${test_name_gtest}" DISCOVERY_MODE PRE_TEST)
endfunction()

foreach(MESH_TYPE "UNIFORM_MESH" "NON_UNIFORM_MESH")
    multipatch_single_interface_derivatives_executable(${MESH_TYPE}) 
endforeach()


# SingleInterfaceDerivativeCalculatorCollection test =============================================
set(test_name_gtest "interface_derivatives_test")
add_executable("${test_name_gtest}"
    single_interface_derivatives_calculator_collection_test.cpp
    ../../main.cpp
)
target_link_libraries("${test_name_gtest}"
PUBLIC
    DDC::core
    GTest::gtest
    GTest::gmock

    gslx::multipatch_geometries
    gslx::multipatch_spline
    gslx::multipatch_interface_derivatives
)
target_compile_definitions("${test_name_gtest}" PUBLIC)
gtest_discover_tests("${test_name_gtest}" DISCOVERY_MODE PRE_TEST)


# InterfaceExactDerivativeMatrix test =================================================================
function (interface_derivative_matrix_executable REVERSED_PATCH)
    # The test with Greville do not depend on the reversed patch. Do not need to duplicate test.
     if (REVERSED_PATCH STREQUAL "REVERSE_PATCH1")
        set(test_name_gtest "interface_derivative_matrix_start_patch_reversed_test")
        add_executable("${test_name_gtest}"
            interface_derivative_matrix_Greville_periodic_test.cpp
            interface_derivative_matrix_Hermite_test.cpp
            ../../main.cpp
        )
    elseif (REVERSED_PATCH STREQUAL "REVERSE_PATCH2")
        set(test_name_gtest "interface_derivative_matrix_middle_patch_reversed_test")
        add_executable("${test_name_gtest}"
            interface_derivative_matrix_Hermite_test.cpp
            ../../main.cpp
        )
    elseif (REVERSED_PATCH STREQUAL "REVERSE_PATCH3")
        set(test_name_gtest "interface_derivative_matrix_end_patch_reversed_test")
        add_executable("${test_name_gtest}"
            interface_derivative_matrix_Hermite_test.cpp
            ../../main.cpp
        )
    else()
        set(test_name_gtest "interface_derivative_matrix_${REVERSED_PATCH}_test")
        add_executable("${test_name_gtest}"
            interface_derivative_matrix_Hermite_test.cpp
            ../../main.cpp
        )
    endif()
    target_link_libraries("${test_name_gtest}"
    PUBLIC
        DDC::core
        GTest::gtest
        GTest::gmock

        gslx::multipatch_geometries
        gslx::multipatch_spline
        gslx::multipatch_interface_derivatives
    )
    target_compile_definitions("${test_name_gtest}" PUBLIC -D${REVERSED_PATCH})
    gtest_discover_tests("${test_name_gtest}" DISCOVERY_MODE PRE_TEST)
endfunction()

interface_derivative_matrix_executable("REVERSE_PATCH1")
interface_derivative_matrix_executable("REVERSE_PATCH2")
interface_derivative_matrix_executable("REVERSE_PATCH3")
interface_derivative_matrix_executable("CHANGE_BOUND1")
interface_derivative_matrix_executable("CHANGE_BOUND3")
