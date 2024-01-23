# SPDX-License-Identifier: MIT
    
add_library("simulation_utils" INTERFACE)

target_compile_features("simulation_utils"
    INTERFACE
        cxx_std_17
)

target_include_directories("simulation_utils"
    INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}"
)

target_link_libraries("simulation_utils"
    INTERFACE
        DDC::DDC
)

add_library("gslx::simulation_utils" ALIAS "simulation_utils")
