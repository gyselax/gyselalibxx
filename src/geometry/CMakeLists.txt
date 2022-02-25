# SPDX-License-Identifier: MIT

cmake_minimum_required(VERSION 3.15)

add_library(geometry INTERFACE)
target_include_directories(geometry
    INTERFACE
        "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>"
)
target_link_libraries(geometry INTERFACE
    DDC::DDC
    sll::splines
)

if(VOICEXX_ENABLE_PERIODIC_RDIMX)
    target_compile_definitions(geometry INTERFACE -DENABLE_PERIODIC_RDIMX=1)
else()
    target_compile_definitions(geometry INTERFACE -DENABLE_PERIODIC_RDIMX=0)
endif()

add_library(vcx::geometry ALIAS geometry)
