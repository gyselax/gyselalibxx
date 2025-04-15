#include <ddc/ddc.hpp>
// #include <ddc/kernels/splines.hpp>

#include <ginkgo/ginkgo.hpp>


// Fill in the matrix.
void generate_stencil_matrix(gko::matrix::Dense<>* matrix)
{
    const auto nb_points = matrix->get_size()[0];
    const double coefs[] = {-1, 2, -1};

    for (int i = 0; i < nb_points; ++i) {
        for (auto dofs : {-1, 0, 1}) {
            if (0 <= i + dofs && i + dofs < nb_points) {
                matrix->at(i, i + dofs) = coefs[dofs + 1];
            }
        }
    }
}


void generate_stencil_matrix(gko::matrix::Csr<>* matrix)
{
    Kokkos::View<
            int*,
            Kokkos::DefaultHostExecutionSpace> // , Kokkos::MemoryTraits<Kokkos::Unmanaged>>
            rows(matrix->get_row_ptrs(), matrix->get_size()[0]);
    Kokkos::View<
            int*,
            Kokkos::DefaultHostExecutionSpace> // , Kokkos::MemoryTraits<Kokkos::Unmanaged>>
            cols(matrix->get_col_idxs(), matrix->get_num_stored_elements());
    Kokkos::View<
            double*,
            Kokkos::DefaultHostExecutionSpace> // , Kokkos::MemoryTraits<Kokkos::Unmanaged>>
            vals(matrix->get_values(), matrix->get_num_stored_elements());


    const auto nb_points = matrix->get_size()[0];
    const double coefs[] = {-1, 2, -1};

    for (int i = 0; i < nb_points; ++i) {
        for (auto dofs : {-1, 0, 1}) {
            if (0 <= i + dofs && i + dofs < nb_points) {
                vals(3*i + dofs) = coefs[dofs + 1];
                cols(3*i + dofs) = i + dofs;
                rows(i) = i-1;
            }
        }
    }
}



int main(int argc, char* argv[])
{
    const unsigned int nb_points = 10;

    using mtx_dense = gko::matrix::Dense<double>;
    using mtx_csr = gko::matrix::Csr<double>;
    // using mtx = gko::matrix::Coo<double>;

    const auto exec = gko::ReferenceExecutor::create();

    // DENSE MATRIX -----------------------------------------------------------
    // Create matrix nb_points x nb_points
    std::shared_ptr<mtx_dense> matrix; 
    matrix = gko::share(mtx_dense::create(exec, gko::dim<2>(nb_points)));
    // std::unique_ptr<mtx_dense> matrix; 
    // matrix = std::make_unique(mtx_dense::create(exec, gko::dim<2>(nb_points)));
    // auto matrix = gko::share(mtx_dense::create(exec, gko::dim<2>(nb_points)));
    // auto vec = gko::share(mtx_dense::create(exec, gko::dim<2>(nb_points, 1)));

    auto vec = mtx_dense::create(exec, gko::dim<2>(nb_points, 1));
    auto sol = mtx_dense::create(exec, gko::dim<2>(nb_points, 1));

    // Use smart pointer.
    generate_stencil_matrix(gko::lend(matrix));
    // generate_stencil_matrix(matrix);

    for (int i = 0; i < nb_points; ++i) {
        for (int j = 0; j < nb_points; ++j) {
            std::cout << matrix->get_values()[i * (nb_points) + j] << "   ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;

    auto values = vec->get_values();
    for (int i = 0; i < nb_points; ++i) {
        values[i] = i * 0.1;
    }

    for (int i = 0; i < nb_points; ++i) {
        std::cout << vec->get_values()[i] << "   ";
        std::cout << std::endl;
    }
    std::cout << std::endl;


    for (int i = 0; i < sol->get_size()[0]; ++i) {
        sol->get_values()[i] = 0.0;
    }

    // SOLVERS ----------------------------------------------------------------
    using cg = gko::solver::Cg<>;
    using bj = gko::preconditioner::Jacobi<>;
    auto solver_factory
            = cg::build()
                      .with_criteria(
                              gko::stop::Iteration::build().with_max_iters(nb_points).on(exec),
                              gko::stop::ResidualNormReduction<>::build()
                                      .with_reduction_factor(1e-6)
                                      .on(exec))
                      .on(exec);
    auto solver_factory2
            = cg::build()
                      .with_criteria(
                              gko::stop::Iteration::build().with_max_iters(nb_points).on(exec),
                              gko::stop::ResidualNormReduction<>::build()
                                      .with_reduction_factor(1e-6)
                                      .on(exec))
                      .with_preconditioner(bj::build().with_max_block_size(8u).on(exec))
                      .on(exec);
    // Create solver
    auto solver = solver_factory->generate(matrix);

    // Solve system
    // solver->apply(gko::lend(vec), gko::lend(sol));
    solver->apply(vec, sol);
    for (int i = 0; i < nb_points; ++i) {
        std::cout << sol->get_values()[i] << "   ";
        std::cout << std::endl;
    }
    std::cout << std::endl;


    // CSR MATRIX -------------------------------------------------------------
    // Kokkos::View<int*, Kokkos::DefaultHostExecutionSpace> rows("rows", nb_points);
    int num_stored_elements = nb_points * 3 - 2; 
    auto matrix_csr = gko::share(mtx_csr::create(exec, gko::dim<2>(nb_points), num_stored_elements));

    // Use smart pointer.
    generate_stencil_matrix(gko::lend(matrix_csr));
    // generate_stencil_matrix(matrix_csr);


    for (int i = 0; i < nb_points; ++i) {
        for (int j = 0; j < nb_points; ++j) {
            if (j < i - 1 || j > i + 1) {
                std::cout << 0 << "   ";
            } else {
                std::cout << matrix_csr->get_values()[3 * i + j-i] << "   ";
            }
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}