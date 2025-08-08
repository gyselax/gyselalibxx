# SPDX-License-Identifier: MIT

find_package(Python3 REQUIRED COMPONENTS Interpreter)

function (diocotron_test_executable PREDCORR_METHOD TIME_METHOD )
    set(test_name "TestSimulationDiocotron_${PREDCORR_METHOD}_${TIME_METHOD}")

    add_test(NAME "${test_name}"
        COMMAND bash "${CMAKE_CURRENT_SOURCE_DIR}/test_diocotron.sh"
            "${PROJECT_SOURCE_DIR}"
            "$<TARGET_FILE:diocotron_${PREDCORR_METHOD}_${TIME_METHOD}>"
            "$<TARGET_FILE:Python3::Interpreter>")
            
    # The test should take less than 15 min.         
    set_property(TEST "${test_name}" PROPERTY TIMEOUT 900000)
    set_property(TEST "${test_name}" PROPERTY COST 100)
    set_property(TEST "${test_name}" PROPERY LABEL Release_only)
endfunction()

# Select the parameters of the simulation:
if(defined ENV{BUILD_END_TO_END_TESTS})
    diocotron_test_executable(PREDCORR EULER_METHOD)
    diocotron_test_executable(EXPLICIT_PREDCORR EULER_METHOD)
    diocotron_test_executable(IMPLICIT_PREDCORR EULER_METHOD)
endif()