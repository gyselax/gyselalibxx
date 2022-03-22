# SPDX-License-Identifier: MIT

add_library(time STATIC
    predcorr.cpp
)

target_compile_features(time
    PUBLIC
        cxx_std_17
)

target_include_directories(time
    PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}"
)

target_link_libraries(time
    PUBLIC
        DDC::DDC
        DDC::PDI_Wrapper
        vcx::speciesinfo
        vcx::poisson
        vcx::vlasov
)

add_library(vcx::time ALIAS time)
