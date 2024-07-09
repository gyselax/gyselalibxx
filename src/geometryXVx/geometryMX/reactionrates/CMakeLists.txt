# SPDX-License-Identifier: MIT

foreach(GEOMETRY_VARIANT IN LISTS GEOMETRY_XVx_VARIANTS_LIST)

add_library("reactionrates_${GEOMETRY_VARIANT}" STATIC
    constantrate.cpp
)

target_include_directories("reactionrates_${GEOMETRY_VARIANT}"
    PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}"
)

target_link_libraries("reactionrates_${GEOMETRY_VARIANT}"
    PUBLIC
        DDC::DDC
        gslx::geometry_${GEOMETRY_VARIANT}
        gslx::speciesinfo
)

add_library("gslx::reactionrates_${GEOMETRY_VARIANT}" ALIAS "reactionrates_${GEOMETRY_VARIANT}")

endforeach()
