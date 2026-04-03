add_library("multipatch_connectivity" INTERFACE)

target_compile_features("multipatch_connectivity" INTERFACE cxx_std_20)

target_sources("multipatch_connectivity"
    INTERFACE
        FILE_SET HEADERS
            FILES
                connectivity_details.hpp
                connectivity.hpp
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
    gslx::coord_transformations
    gslx::data_types
    gslx::multipatch_data_types
    gslx::utils
)
add_library("gslx::multipatch_connectivity" ALIAS "multipatch_connectivity")

install(TARGETS multipatch_connectivity EXPORT GyselalibxxTargets FILE_SET HEADERS)
