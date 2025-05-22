# To remove ======================================================================================
set(test_name_gtest "Ginkgo_test")
add_executable("${test_name_gtest}"
    test_with_ginkgo.cpp
    #../../main.cpp
)
target_link_libraries("${test_name_gtest}"
PUBLIC
    DDC::core
    gslx::utils
    # GTest::gtest
    # GTest::gmock
)
target_compile_definitions("${test_name_gtest}" PUBLIC)
gtest_discover_tests("${test_name_gtest}" DISCOVERY_MODE PRE_TEST)
# ================================================================================================


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
