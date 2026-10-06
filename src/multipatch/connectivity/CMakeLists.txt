add_library("multipatch_connectivity" INTERFACE)

target_compile_features("multipatch_connectivity" INTERFACE cxx_std_20)

target_sources("multipatch_connectivity"
    INTERFACE
        FILE_SET HEADERS
            FILES
                connectivity.hpp
                connectivity_details.hpp
                edge.hpp
                edge_transformation.hpp
                interface.hpp
                matching_idx_slice.hpp
                onion_patch_locator.hpp
                patch.hpp
                utils_patch_locators.hpp
)

target_link_libraries("multipatch_connectivity" INTERFACE
    DDC::core
    gyselalibxx::coord_transformations
    gyselalibxx::data_types
    gyselalibxx::multipatch_data_types
    gyselalibxx::utils
)
add_library("gyselalibxx::multipatch_connectivity" ALIAS "multipatch_connectivity")

install(TARGETS multipatch_connectivity EXPORT gyselalibxxTargets FILE_SET HEADERS)
