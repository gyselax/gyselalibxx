# SPDX-License-Identifier: MIT

add_library("utils" STATIC
            assert.cpp)

target_compile_features("utils"
    INTERFACE
        cxx_std_17
)

target_include_directories("utils"
    INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}"
)

target_link_libraries("utils"
    INTERFACE
        DDC::DDC
        gslx::data_types
)

add_library("gslx::utils" ALIAS "utils")

