# SPDX-License-Identifier: MIT

find_package(Python3 REQUIRED COMPONENTS Interpreter)

function (diocotron_executable PREDCORR_METHOD TIME_METHOD )
    set(test_name "TestSimulationDiocotron_${PREDCORR_METHOD}_${TIME_METHOD}")

    add_test(NAME "${test_name}"
        COMMAND bash "${CMAKE_CURRENT_SOURCE_DIR}/test_diocotron.sh"
            "${PROJECT_SOURCE_DIR}"
            "$<TARGET_FILE:diocotron_${PREDCORR_METHOD}_${TIME_METHOD}>"
            "$<TARGET_FILE:Python3::Interpreter>")
            
    set_property(TEST "${test_name}" PROPERTY TIMEOUT 10000)
    set_property(TEST "${test_name}" PROPERTY COST 100)
endfunction()

# Select the parameters of the simulation:
diocotron_executable(PREDCORR RK3_METHOD)
diocotron_executable(EXPLICIT_PREDCORR EULER_METHOD)
diocotron_executable(IMPLICIT_PREDCORR EULER_METHOD)
