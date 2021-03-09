add_library(math_tools INTERFACE)

target_include_directories( math_tools
    INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}"
)
add_library( selalib::math_tools ALIAS math_tools )
