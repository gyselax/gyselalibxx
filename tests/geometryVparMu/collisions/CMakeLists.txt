# SPDX-License-Identifier: MIT

find_package(Python3 REQUIRED COMPONENTS Interpreter)

add_test(NAME TestCollisionsVparMuTwoSpeciesDeltatZero
    COMMAND bash "${CMAKE_CURRENT_SOURCE_DIR}/test_collisions_vparmu_twospecies_deltat0.sh"
        "$<TARGET_FILE:test_collSpVparMu>")
set_property(TEST TestCollisionsVparMuTwoSpeciesDeltatZero PROPERTY TIMEOUT 200)
set_property(TEST TestCollisionsVparMuTwoSpeciesDeltatZero PROPERTY COST 10)

