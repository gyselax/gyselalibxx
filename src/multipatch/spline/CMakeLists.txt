

add_library("multipatch_spline" INTERFACE)

target_sources("multipatch_spline"
    INTERFACE
        FILE_SET HEADERS
            FILES
                constant_extrapolation_rules_onion.hpp
                multipatch_spline_builder_2d.hpp
                multipatch_spline_builder.hpp
                multipatch_spline_evaluator_2d.hpp
                null_extrapolation_rules.hpp
)

target_link_libraries("multipatch_spline" INTERFACE
    DDC::core
    gslx::multipatch_data_types
    gslx::multipatch_connectivity
    gslx::utils

)
add_library("gslx::multipatch_spline" ALIAS "multipatch_spline")

install(TARGETS multipatch_spline EXPORT GyselalibxxTargets FILE_SET HEADERS)
