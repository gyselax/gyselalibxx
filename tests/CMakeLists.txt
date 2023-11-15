# SPDX-License-Identifier: MIT

include(GoogleTest)

add_executable(unit_tests_common
    main.cpp
    species_info.cpp
)
target_compile_features(unit_tests_common PUBLIC cxx_std_17)
target_link_libraries(unit_tests_common
    PUBLIC
        GTest::gmock
        gslx::speciesinfo
)

gtest_discover_tests(unit_tests_common)

add_subdirectory(geometryXVx)
add_subdirectory(geometryXYVxVy)
add_subdirectory(geometryRTheta)
add_subdirectory(timestepper)
add_subdirectory(utils)
