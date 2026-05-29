# SPDX-License-Identifier: MIT

foreach(PRECISION_VARIANT IN LISTS PRECISION_VARIANTS_LIST)

add_library("time_integration_xyvxvy_${PRECISION_VARIANT}" STATIC
    predcorr.cpp
)

target_include_directories("time_integration_xyvxvy_${PRECISION_VARIANT}"
    PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}"
)

target_link_libraries("time_integration_xyvxvy_${PRECISION_VARIANT}"
    PUBLIC
        DDC::core
        DDC::pdi
        gslx::poisson_xy_${PRECISION_VARIANT}
        gslx::speciesinfo
        gslx::vlasov_xyvxvy_${PRECISION_VARIANT}
        gslx::utils

)

add_library("gslx::time_integration_xyvxvy_${PRECISION_VARIANT}" ALIAS "time_integration_xyvxvy_${PRECISION_VARIANT}")

endforeach()

