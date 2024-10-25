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
        gslx::utils

)

gtest_discover_tests(unit_tests_parallelisation_layout DISCOVERY_MODE PRE_TEST)

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
        gslx::utils

)

function(make_mpi_test test_name)
    add_test(NAME ${test_name}
        COMMAND
        "${MPIEXEC_EXECUTABLE}"
        "-n"
        "2"
        "$<TARGET_FILE:unit_tests_parallelisation>"
        "--gtest_filter=${test_name}"
    )
endfunction()

make_mpi_test(MPIParallelisation.AllToAll2D_CPU)
make_mpi_test(MPIParallelisation.AllToAll2D_GPU)
make_mpi_test(MPIParallelisation.AllToAll3D_CPU)
make_mpi_test(MPIParallelisation.AllToAll4D_CPU)
make_mpi_test(Layout.MinimalDomainDistribution)
make_mpi_test(Layout.SpreadDomainDistribution)
make_mpi_test(Layout.DomainSelection)
