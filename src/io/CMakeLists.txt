# SPDX-License-Identifier: MIT

add_library("io"
  STATIC
    input.cpp
)

target_include_directories("io"
    INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}"
)

target_link_libraries("io"
    PUBLIC
        DDC::core
        PDI::pdi
        gslx::paraconfpp
        gslx::utils

)

add_library("gslx::io" ALIAS "io")
