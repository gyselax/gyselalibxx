# SPDX-License-Identifier: MIT
add_library(gyroaverage INTERFACE)

target_compile_features("gyroaverage" INTERFACE cxx_std_20)

target_sources("gyroaverage"
    INTERFACE
        FILE_SET HEADERS
            FILES
                gyroaverage_operator.hpp
)

target_link_libraries(gyroaverage
    INTERFACE
        DDC::core
        gslx::utils
        gslx::coord_transformations
)

add_library(gslx::gyroaverage ALIAS gyroaverage)

install(TARGETS gyroaverage EXPORT GyselalibxxTargets FILE_SET HEADERS)
