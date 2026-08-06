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
        gyselalibxx::data_types
        gyselalibxx::paraconfpp
        gyselalibxx::utils
)

add_library("gyselalibxx::io" ALIAS "io")

install(TARGETS io EXPORT GyselalibxxTargets FILE_SET HEADERS)
