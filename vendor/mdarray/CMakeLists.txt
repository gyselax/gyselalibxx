
cmake_minimum_required(VERSION 3.12)
project(MDArray
  VERSION 0.0.1
  LANGUAGES CXX
)

################################################################################

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

################################################################################

option(MDARRAY_ENABLE_TESTS "Enable tests." Off)
option(MDARRAY_ENABLE_EXAMPLES "Build examples." Off)
#option(MDARRAY_ENABLE_BENCHMARKS "Enable benchmarks." Off)
#option(MDARRAY_ENABLE_COMP_BENCH "Enable compilation benchmarks." Off)

# Option to override which C++ standard to use
set(MDARRAY_CXX_STANDARD DETECT CACHE STRING "Override the default CXX_STANDARD to compile with.")
set_property(CACHE MDARRAY_CXX_STANDARD PROPERTY STRINGS DETECT 11 14 17 20)

option(MDARRAY_ENABLE_CONCEPTS "Try to enable concepts support by giving extra flags." On)

################################################################################

# Decide on the standard to use
if(MDARRAY_CXX_STANDARD STREQUAL "17")
  if("cxx_std_17" IN_LIST CMAKE_CXX_COMPILE_FEATURES)
    message(STATUS "Using C++17 standard")
    set(CMAKE_CXX_STANDARD 17)
  else()
    message(FATAL_ERROR "Requested MDARRAY_CXX_STANDARD \"17\" not supported by provided C++ compiler")
  endif()
elseif(MDARRAY_CXX_STANDARD STREQUAL "14")
  if("cxx_std_14" IN_LIST CMAKE_CXX_COMPILE_FEATURES)
    message(STATUS "Using C++14 standard")
    set(CMAKE_CXX_STANDARD 14)
  else()
    message(FATAL_ERROR "Requested MDARRAY_CXX_STANDARD \"14\" not supported by provided C++ compiler")
  endif()
elseif(MDARRAY_CXX_STANDARD STREQUAL "11")
  message(STATUS "Using C++11 standard")
  set(CMAKE_CXX_STANDARD 11)
elseif(MDARRAY_CXX_STANDARD STREQUAL "20")
  if("cxx_std_20" IN_LIST CMAKE_CXX_COMPILE_FEATURES)
    message(STATUS "Using C++20 standard")
    set(CMAKE_CXX_STANDARD 20)
  else()
    message(FATAL_ERROR "Requested MDARRAY_CXX_STANDARD \"20\" not supported by provided C++ compiler")
  endif()
else()
  if("cxx_std_20" IN_LIST CMAKE_CXX_COMPILE_FEATURES)
    set(CMAKE_CXX_STANDARD 20)
    message(STATUS "Detected support for C++20 standard")
  elseif("cxx_std_17" IN_LIST CMAKE_CXX_COMPILE_FEATURES)
    set(CMAKE_CXX_STANDARD 17)
    message(STATUS "Detected support for C++17 standard")
  elseif("cxx_std_14" IN_LIST CMAKE_CXX_COMPILE_FEATURES)
    set(CMAKE_CXX_STANDARD 14)
    message(STATUS "Detected support for C++14 standard")
  elseif("cxx_std_11" IN_LIST CMAKE_CXX_COMPILE_FEATURES)
    set(CMAKE_CXX_STANDARD 11)
    message(STATUS "Newest supported C++ standard detected is C++11.")
  else()
    message(FATAL_ERROR "Cannot detect CXX_STANDARD of C++11 or newer.")
  endif()
endif()

################################################################################

if(MDARRAY_ENABLE_CONCEPTS)
  if(CMAKE_CXX_STANDARD STREQUAL "20")
    include(CheckCXXCompilerFlag)
    CHECK_CXX_COMPILER_FLAG("-fconcepts" COMPILER_SUPPORTS_FCONCEPTS)
    if(COMPILER_SUPPORTS_FCONCEPTS)
      message(STATUS "Using \"-fconcepts\" to enable concepts support")
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fconcepts")
    else()
      CHECK_CXX_COMPILER_FLAG("-fconcepts-ts" COMPILER_SUPPORTS_FCONCEPTS_TS)
      if(COMPILER_SUPPORTS_FCONCEPTS)
        message(STATUS "Using \"-fconcepts-ts\" to enable concepts support")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fconcepts-ts")
      endif()
    endif()
    # Otherwise, it's possible that the compiler supports concepts without flags,
    # but if it doesn't, they just won't be used, which is fine
  endif()
endif()

################################################################################

if(NOT TARGET std::mdspan)
  find_package(MDSpan REQUIRED)
  message(STATUS "Found MDSpan in ${MDSpan_DIR}")
endif()


################################################################################

add_library(mdarray INTERFACE)
add_library(std::mdarray ALIAS mdarray)

target_link_libraries(mdarray INTERFACE std::mdspan)

target_include_directories(mdarray INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
)

################################################################################

install(TARGETS mdarray EXPORT mdarrayTargets
  INCLUDES DESTINATION include
)

install(EXPORT mdarrayTargets
    FILE mdarrayTargets.cmake
    NAMESPACE std::
    DESTINATION cmake
)

export(TARGETS mdarray
    NAMESPACE std::
    FILE mdarrayTargets.cmake
)

install(DIRECTORY include/experimental DESTINATION include)

include(CMakePackageConfigHelpers)
configure_package_config_file(cmake/MDArrayConfig.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/MDArrayConfig.cmake
  INSTALL_DESTINATION cmake
)
write_basic_package_version_file(${CMAKE_CURRENT_BINARY_DIR}/MDArrayConfigVersion.cmake
  COMPATIBILITY SameMajorVersion
)

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/MDArrayConfig.cmake ${CMAKE_CURRENT_BINARY_DIR}/MDArrayConfigVersion.cmake
    DESTINATION cmake
)

################################################################################

if(MDARRAY_ENABLE_TESTS)
  enable_testing()
  add_subdirectory(tests)
  add_subdirectory(compilation_tests)
endif()

if(MDARRAY_ENABLE_EXAMPLES)
 add_subdirectory(examples)
endif()

#if(MDARRAY_ENABLE_BENCHMARKS)
# add_subdirectory(benchmarks)
#endif()

#if(MDARRAY_ENABLE_COMP_BENCH)
#  add_subdirectory(comp_bench)
#endif()
