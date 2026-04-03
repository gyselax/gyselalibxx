
add_library("multipatch_interface_derivatives" INTERFACE)

target_compile_features("multipatch_interface_derivatives" INTERFACE cxx_std_20)

target_sources("multipatch_interface_derivatives"
    INTERFACE
        FILE_SET HEADERS
            FILES
                single_interface_derivatives_calculator_collection.hpp
                single_interface_derivatives_calculator.hpp
)

target_link_libraries("multipatch_interface_derivatives" INTERFACE
    DDC::core
    gslx::multipatch_data_types
    gslx::multipatch_connectivity
    gslx::utils

)
add_library("gslx::multipatch_interface_derivatives" ALIAS "multipatch_interface_derivatives")

install(TARGETS multipatch_interface_derivatives EXPORT GyselalibxxTargets FILE_SET HEADERS)
