# SPDX-License-Identifier: MIT

include(GoogleTest)

add_executable(unit_tests_parallelisation_layout
    layout.cpp
    ../main.cpp
)

target_link_libraries(unit_tests_parallelisation_layout
    PUBLIC
        DDC::DDC
        GTest::gtest
        GTest::gmock
        gslx::mpi_parallelisation
)

gtest_discover_tests(unit_tests_parallelisation_layout)

add_executable(unit_tests_parallelisation
    alltoall.cpp
    layout.cpp
    main.cpp
)
target_link_libraries(unit_tests_parallelisation
    PUBLIC
        DDC::DDC
        GTest::gtest
        GTest::gmock
        gslx::mpi_parallelisation
)

set_property(TARGET unit_tests_parallelisation PROPERTY CROSSCOMPILING_EMULATOR '${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2')

gtest_discover_tests(unit_tests_parallelisation)

