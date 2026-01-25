# src/exact_splitting_solver/CMakeLists.txt

# Find the system-installed Eigen
find_package(Eigen3 3.3 REQUIRED NO_MODULE)

# Create a static library for the solver
add_library(exact_splitting_solver STATIC
    exact_splitting_solver.cpp
)

# Add the header file include path
target_include_directories(exact_splitting_solver
    PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}"
)

# Link Eigen and any other necessary libraries
target_link_libraries(exact_splitting_solver
    PUBLIC Eigen3::Eigen
)

# Optional: create an alias for easier linking in other targets
add_library(gslx::exact_splitting_solver ALIAS exact_splitting_solver)

