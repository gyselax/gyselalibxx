# SPDX-License-Identifier: MIT

find_package(Python3 REQUIRED COMPONENTS Interpreter)

add_test(NAME TestCollisionsVparMuTwoSpeciesDeltatZero
    COMMAND bash "${CMAKE_CURRENT_SOURCE_DIR}/test_collisions_vparmu_twospecies_deltat0.sh"
        "${PROJECT_SOURCE_DIR}"
        "$<TARGET_FILE:test_collSpVparMu>"
        "$<TARGET_FILE:Python3::Interpreter>")
set_property(TEST TestCollisionsVparMuTwoSpeciesDeltatZero PROPERTY TIMEOUT 20)
set_property(TEST TestCollisionsVparMuTwoSpeciesDeltatZero PROPERTY COST 10)

