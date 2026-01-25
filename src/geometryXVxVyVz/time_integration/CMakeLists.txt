# SPDX-License-Identifier: MIT

add_library("time_integration_xvxvyvz" STATIC
    hybridsplitting.cpp
)

target_include_directories("time_integration_xvxvyvz"
    PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}"
)

target_link_libraries("time_integration_xvxvyvz"
    PUBLIC
        DDC::core
        DDC::pdi
        gslx::hybrid_field_solver_xvxvyvz
        gslx::speciesinfo
        gslx::vlasov_xvxvyvz
        gslx::utils
        gslx::geometry_xvxvyvz

)

add_library("gslx::time_integration_xvxvyvz" ALIAS "time_integration_xvxvyvz")

