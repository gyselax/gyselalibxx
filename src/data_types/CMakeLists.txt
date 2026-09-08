# SPDX-License-Identifier: MIT

add_library("data_types" INTERFACE)

target_compile_features("data_types" INTERFACE cxx_std_20)

target_sources("data_types"
    INTERFACE
        FILE_SET HEADERS
            FILES
                deriv_details.hpp
                derivative_field.hpp
                derivative_field_common.hpp
                derivative_field_mem.hpp
                idx_range_slice.hpp
                indexed_tensor.hpp
                static_tensors.hpp
                tensor.hpp
                tensor_common.hpp
                tensor_index_tools.hpp
                vector_field.hpp
                vector_field_common.hpp
                vector_field_mem.hpp
                vector_index_tools.hpp
)

target_link_libraries("data_types"
    INTERFACE
        DDC::core
        DDC::splines
        gyselalibxx::utils
)

add_library("gyselalibxx::data_types" ALIAS "data_types")

install(TARGETS data_types EXPORT gyselalibxxTargets FILE_SET HEADERS)
