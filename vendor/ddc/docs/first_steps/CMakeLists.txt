# SPDX-License-Identifier: MIT

if("${DDC_BUILD_PDI_WRAPPER}")

add_executable(heat_equation heat_equation.cpp)
target_link_libraries(heat_equation PUBLIC DDC::DDC DDC::PDI)
endif()
