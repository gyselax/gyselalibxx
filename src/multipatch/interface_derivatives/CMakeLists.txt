
add_library("multipatch_interface_derivatives" INTERFACE)

target_compile_features("multipatch_interface_derivatives" INTERFACE cxx_std_20)

target_sources("multipatch_interface_derivatives"
    INTERFACE
        FILE_SET HEADERS
            FILES
                single_interface_derivatives_calculator.hpp
                single_interface_derivatives_calculator_collection.hpp
)

target_link_libraries("multipatch_interface_derivatives" INTERFACE
    DDC::core
    gyselalibxx::multipatch_data_types
    gyselalibxx::multipatch_connectivity
    gyselalibxx::utils

)
add_library("gyselalibxx::multipatch_interface_derivatives" ALIAS "multipatch_interface_derivatives")

install(TARGETS multipatch_interface_derivatives EXPORT gyselalibxxTargets FILE_SET HEADERS)
