enum class BoundaryCondition {
    // User defined boundary condition
    sll_p_user_defined = -1,
    // Periodic boundary condition u(1)=u(n)
    sll_p_periodic = 0,
    // Dirichlet boundary condition
    sll_p_dirichlet = 1,
    // Neumann boundary condition
    sll_p_neumann = 2,
    // Hermite boundary condition
    sll_p_hermite = 3,
    // Neumann boundary condition
    sll_p_neumann_mode_0 = 4,
    // PLEASE ADD DOCUMENTATION
    sll_p_set_to_limit = 5,
    // Interior of domain
    sll_p_interior = 6,
    // Incoming wave boundar condition for Maxwell
    SLL_INCOMING_WAVE = 7,
    // Metallic boundary condition for Maxwell
    sll_p_conductor = 8,
    // Absorbing boundary condition fro Maxwell
    sll_p_silver_muller = 9,
    // Use a one-sided stencil at the boundary
    sll_p_one_sided = 10,
    // Values outside the domain are provided as halo cells (for domain
    // decomposition)
    sll_p_halo = 11,
    // Use Greville points instead of conditions on derivative for B-Spline
    // interpolation
    sll_p_greville = 12,
    // Duplicate boundary points to define knots for B-splines from grid
    sll_p_open = 13,
    // Mirror points around boundary to define knots for B-splines from grid
    sll_p_mirror = 14,
    // Treat boundary poas origin of polar coordinates (= center of a circle)
    sll_p_polar_origin = 15
};
