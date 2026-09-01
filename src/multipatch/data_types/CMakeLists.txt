
add_library("multipatch_data_types" INTERFACE)

target_compile_features("multipatch_data_types" INTERFACE cxx_std_20)

target_sources("multipatch_data_types"
    INTERFACE
        FILE_SET HEADERS
            FILES
                multipatch_field.hpp
                multipatch_field_mem.hpp
                multipatch_type.hpp
                types.hpp
)

target_link_libraries("multipatch_data_types" INTERFACE
    DDC::core
    DDC::splines

    gyselalibxx::data_types
    gyselalibxx::utils
)
add_library("gyselalibxx::multipatch_data_types" ALIAS "multipatch_data_types")

install(TARGETS multipatch_data_types EXPORT gyselalibxxTargets FILE_SET HEADERS)
