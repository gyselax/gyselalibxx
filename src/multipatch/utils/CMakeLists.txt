# SPDX-License-Identifier: MIT

add_library("multipatch_utils" INTERFACE)

target_sources("multipatch_utils"
    INTERFACE
        FILE_SET HEADERS
            FILES
                multipatch_math_tools.hpp
)

target_link_libraries("multipatch_utils"
    INTERFACE
        gslx::multipatch_data_types
        gslx::utils
)

add_library("gslx::multipatch_utils" ALIAS "multipatch_utils")

install(TARGETS multipatch_utils EXPORT GyselalibxxTargets FILE_SET HEADERS)
