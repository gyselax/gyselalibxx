add_library(geometry INTERFACE)

target_include_directories( geometry
    INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}"
)

add_library(vcx::geometry ALIAS geometry)
