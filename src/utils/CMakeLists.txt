# SPDX-License-Identifier: MIT

add_library("utils" STATIC
            assert.cpp)

target_sources("utils"
    INTERFACE
        FILE_SET HEADERS
            FILES
                assert.hpp
                ddc_aliases.hpp
                ddc_alias_inline_functions.hpp
                ddc_helper.hpp
                geometry_descriptors.hpp
                mesh_builder.hpp
                non_uniform_interpolation_points.hpp
                preprocessor.hpp
                spline_builder_deriv_field_2d.hpp
                transpose.hpp
                type_seq_tools.hpp
                view.hpp
)

target_link_libraries("utils"
    INTERFACE
        DDC::core
        DDC::splines
)

add_library("gslx::utils" ALIAS "utils")

install(TARGETS utils EXPORT GyselalibxxTargets FILE_SET HEADERS)
