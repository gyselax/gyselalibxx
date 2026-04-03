# SPDX-License-Identifier: MIT

add_library("io"
  STATIC
    input.cpp
)

target_compile_features("io" PUBLIC cxx_std_20)

target_sources("io"
    INTERFACE
        FILE_SET HEADERS
            FILES
                input.hpp
                output.hpp
                pdi_helper.hpp
)

target_link_libraries("io"
    PUBLIC
        DDC::core
        PDI::PDI_C
        gslx::data_types
        gslx::paraconfpp
        gslx::utils
)

add_library("gslx::io" ALIAS "io")

install(TARGETS io EXPORT GyselalibxxTargets FILE_SET HEADERS)
