cmake_minimum_required(VERSION 3.15)
project(SLL CXX)

# Our project
option(SLL_BUILD_TESTING "Build sll tests" ON)

## The library itself
add_library(SLL INTERFACE)
target_include_directories(SLL
    INTERFACE
        "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
)
target_link_libraries(SLL
    INTERFACE
        DDC::core
        DDC::splines
        Kokkos::kokkos
)
add_library(sll::SLL ALIAS SLL)

## if tests are enabled, build the tests in `tests/`
if("${GYSELALIBXX_BUILD_TESTING}" AND "${SLL_BUILD_TESTING}")
    add_subdirectory(tests/)
endif()
