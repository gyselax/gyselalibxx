# SPDX-License-Identifier: MIT

add_library("io" INTERFACE)

target_compile_features("io"
    INTERFACE
        cxx_std_17
)

target_include_directories("io"
    INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}"
)

target_link_libraries("io"
    INTERFACE
        DDC::DDC
)

add_library("gslx::io" ALIAS "io")