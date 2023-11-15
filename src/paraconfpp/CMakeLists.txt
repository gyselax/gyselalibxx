# SPDX-License-Identifier: MIT

add_library(paraconfpp INTERFACE)
target_include_directories(paraconfpp
    INTERFACE
        ${CMAKE_CURRENT_SOURCE_DIR}
)
target_link_libraries(paraconfpp INTERFACE
        paraconf::paraconf
)
add_library(gslx::paraconfpp ALIAS paraconfpp)
