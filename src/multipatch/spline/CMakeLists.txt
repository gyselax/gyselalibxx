

add_library("multipatch_spline" INTERFACE)

target_compile_features("multipatch_spline" INTERFACE cxx_std_20)

target_sources("multipatch_spline"
    INTERFACE
        FILE_SET HEADERS
            FILES
                constant_extrapolation_rules_onion.hpp
                multipatch_spline_builder.hpp
                multipatch_spline_builder_2d.hpp
                multipatch_spline_evaluator_2d.hpp
                null_extrapolation_rules.hpp
)

target_link_libraries("multipatch_spline" INTERFACE
    DDC::core
    gyselalibxx::multipatch_data_types
    gyselalibxx::multipatch_connectivity
    gyselalibxx::utils

)
add_library("gyselalibxx::multipatch_spline" ALIAS "multipatch_spline")

install(TARGETS multipatch_spline EXPORT gyselalibxxTargets FILE_SET HEADERS)
