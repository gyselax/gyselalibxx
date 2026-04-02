# SPDX-License-Identifier: MIT

add_library("coord_transformations" INTERFACE)

target_sources("coord_transformations"
    INTERFACE
        FILE_SET HEADERS
            FILES
                barycentric_to_cartesian.hpp
                cartesian_to_barycentric.hpp
                cartesian_to_circular.hpp
                cartesian_to_cylindrical.hpp
                cartesian_to_czarny.hpp
                circular_to_cartesian.hpp
                combined_mapping.hpp
                coord_transformation_tools.hpp
                cylindrical_to_cartesian.hpp
                czarny_to_cartesian.hpp
                discrete_poloidal_cs_spline_mapping_builder.hpp
                discrete_poloidal_cs_spline_mapping.hpp
                geometry_pseudo_cartesian.hpp
                identity_coordinate_change.hpp
                inverse_jacobian_matrix.hpp
                inv_jacobian_o_point.hpp
                linear_coord_transform.hpp
                metric_tensor_evaluator.hpp
                orthogonal_coord_transforms.hpp
                toroidal_to_cylindrical.hpp
                vector_mapper.hpp
)

target_link_libraries("coord_transformations"
    INTERFACE
        DDC::core
        gslx::data_types
        gslx::math_tools
        gslx::utils
)

add_library("gslx::coord_transformations" ALIAS "coord_transformations")

install(TARGETS coord_transformations EXPORT GyselalibxxTargets FILE_SET HEADERS)
