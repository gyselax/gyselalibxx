#include <ddc/ddc.hpp>
#include <ddc/kernels/splines.hpp>

#include <gtest/gtest.h>

#include "discrete_poloidal_cs_spline_mapping_builder.hpp"
#include "initialize_splines_and_grids.hpp"
#include "interface.hpp"
#include "polar_spline_fem_poisson_like_solver.hpp"
#include "projection_matrix_utilities.hpp" // TODO: Change the name
#include "spline_interpolation.hpp"
using namespace onion_2d_2patches_uniform_unique_dims;

namespace {

class PoissonMatrixOnionTwoPatches : public ::testing::Test
{
protected:
    static constexpr Patch1::IdxStep1 r1_ncells = Patch1::IdxStep1(10);

    static constexpr std::size_t p_size = 10;
    static constexpr Patch1::IdxStep2 theta1_ncells = Patch1::IdxStep2(p_size);

    static constexpr Patch1::Coord1 r1_min {0.0};
    static constexpr Patch1::Coord1 r1_max {1.0};

    static constexpr Patch1::Coord2 theta1_min {0.0};
    static constexpr Patch1::Coord2 theta1_max {2 * M_PI};

    using MatrixType = gko::matrix::Csr<double>;

    using Mapping = CircularToCartesian<Patch1::Dim1, Patch1::Dim2, X, Y>;

    static Patch1::IdxRange12 idx_range_rtheta1;

    static std::shared_ptr<SplineBuilderOnPatch<Patch1>> builder1;
    MultipatchBuilder builders;

    Patch1::IdxRangeBS12 idx_range_bspline1_with_reps;
    IdxRange<PolarBSplinesPatch1> idx_range_spline_polar;

    Patch1::IdxRangeBS12 const idx_range_bspline1;

    DFieldMem<Patch1::IdxRangeBS12> spline_coeffs_with_reps_mem1;
    DField<Patch1::IdxRangeBS12> spline_coeffs_with_reps1;
    MultipatchSplineCoeff_host spline_coeffs_with_reps;

    DFieldMem<Patch1::IdxRangeBS12> spline_coeffs_mem1;
    DField<Patch1::IdxRangeBS12> spline_coeffs1;
    MultipatchSplineCoeff_host spline_coeffs;

    SplinePolarMemPatch1 polar_spline_mem;
    SplinePolarPatch1 polar_spline;
    DField<Patch1::IdxRangeBS12> spline_coeffs_polar_tensor;
    SplinePolarPatch1 spline_coeffs_polar_singular;
    DFieldMem<Patch1::IdxRange12> field_mem1;
    DField<Patch1::IdxRange12> field1;

    DMultipatchField fields;

    DFieldMem<Patch1::IdxRangeBS12> coeff_alpha_mem;
    DField<Patch1::IdxRangeBS12> coeff_alpha;
    DFieldMem<Patch1::IdxRangeBS12> coeff_beta_mem;
    DField<Patch1::IdxRangeBS12> coeff_beta;



    using SplineRThetaEvaluatorNullBound = ddc::SplineEvaluator2D<
            Kokkos::DefaultExecutionSpace,
            typename Kokkos::DefaultExecutionSpace::memory_space,
            Patch1::BSplines1,
            Patch1::BSplines2,
            Patch1::Grid1,
            Patch1::Grid2,
            ddc::NullExtrapolationRule, // boundary at r=0
            ddc::NullExtrapolationRule, // boundary at rmax
            ddc::PeriodicExtrapolationRule<Patch1::Dim2>,
            ddc::PeriodicExtrapolationRule<Patch1::Dim2>>;

    static std::shared_ptr<gko::ReferenceExecutor> exec;
    static Mapping const mapping;
    static SplineEvaluatorOnPatch<Patch1> const evaluator;

    using SplineInterpolatorPatch1 = SplineInterpolator<
            Kokkos::DefaultHostExecutionSpace,
            IdxRange<typename Patch1::BSplines1, typename Patch1::BSplines2>,
            Patch1::IdxRange12,
            ExtrapolationRule::Constant_Constant, // radial extrapolation
            ExtrapolationRule::Periodic, // poloidal extrapolation
            SplineBoundaryClosures<
                    ddc::SplineBuilderClosure::GREVILLE, // boundary at r=0
                    ddc::SplineBuilderClosure::GREVILLE>, // boundary at rmax
            SplineBoundaryClosures<
                    ddc::SplineBuilderClosure::PERIODIC,
                    ddc::SplineBuilderClosure::PERIODIC>>;

public:
    static void SetUpTestSuite()
    {
        initialize_splines_and_grids(
                r1_min,
                r1_max,
                r1_ncells,
                theta1_min,
                theta1_max,
                theta1_ncells,
                r2_min,
                r2_max,
                r2_ncells,
                theta2_min,
                theta2_max,
                theta2_ncells);

        idx_range_rtheta1 = Patch1::IdxRange12 {
                SplineInterpPointsR<1>::get_domain<Patch1::Grid1>(),
                SplineInterpPointsTheta<1>::get_domain<Patch1::Grid2>()};
        SplineInterpolatorPatch1 interpolator(idx_range_rtheta1);
        DiscretePoloidalCSSplineMappingBuilder<X, Y, SplineInterpolatorPatch1> const
                discrete_mapping_builder(
                        Kokkos::DefaultHostExecutionSpace(),
                        mapping,
                        interpolator);
        DiscretePoloidalCSSplineMapping const discrete_mapping = discrete_mapping_builder();

        ddc::init_discrete_space<PolarBSplinesPatch1>(discrete_mapping);
        builder1 = std::make_shared<SplineBuilderOnPatch<Patch1>>(idx_range_rtheta1);
    }

    PoissonMatrix()
        , idx_range_bspline1_with_reps(get_spline_idx_range(*builder1))
        , idx_range_spline_polar(ddc::discrete_space<PolarBSplinesPatch1>().full_domain())
        , idx_range_bspline1(
                  ddc::discrete_space<Patch1::BSplines1>().full_domain(),
                  ddc::discrete_space<Patch1::BSplines2>().full_domain().take_first(
                          IdxStep<Patch1::BSplines2> {
                                  ddc::discrete_space<Patch1::BSplines2>().nbasis()}))
        , spline_coeffs_with_reps_mem1(idx_range_bspline1_with_reps)
        , spline_coeffs_with_reps1(get_field(spline_coeffs_with_reps_mem1))
        , spline_coeffs_mem1(idx_range_bspline1)
        , spline_coeffs1(get_field(spline_coeffs_mem1))
        , spline_coeffs(spline_coeffs1, spline_coeffs2)
        , polar_spline_mem(idx_range_spline_polar)
        , polar_spline(get_field(polar_spline_mem))
        , spline_coeffs_polar_tensor(PolarBSplinesPatch1::get_tensor_product_subset(polar_spline))
        , spline_coeffs_polar_singular(PolarBSplinesPatch1::get_singular_subset(polar_spline))
        , field_mem1(idx_range_rtheta1)
        , field1(get_field(field_mem1))
        , fields(field1, field2)
        , coeff_alpha_mem {idx_range_bspline1_with_reps}
        , coeff_alpha {get_field(coeff_alpha_mem)}
        , coeff_beta_mem {idx_range_bspline1_with_reps}
        , coeff_beta {get_field(coeff_beta_mem)}
    {
        ddc::parallel_fill(coeff_alpha, 1.0);
        ddc::parallel_fill(coeff_beta, 0.0);
    }
};

std::shared_ptr<gko::ReferenceExecutor> PoissonMatrixOnionTwoPatches::exec
        = gko::ReferenceExecutor::create();
PoissonMatrixOnionTwoPatches::Mapping const PoissonMatrixOnionTwoPatches::mapping;

ddc::NullExtrapolationRule bv_r_min;
ddc::NullExtrapolationRule bv_r_max;
ddc::PeriodicExtrapolationRule<Theta<1>> bv_theta_min;
ddc::PeriodicExtrapolationRule<Theta<1>> bv_theta_max;
SplineEvaluator const evaluator(
        bv_r_min,
        bv_r_max,
        bv_theta_min,
        bv_theta_max);

Patch1::IdxRange12 PoissonMatrixOnionTwoPatches::idx_range_rtheta1;
std::shared_ptr<SplineBuilderOnPatch<Patch1>> PoissonMatrixOnionTwoPatches::builder1;


} // namespace


// TODO: Lots of duplication

TEST_F(PoissonMatrixOnionTwoPatches, H1NormConstantPolar)
{
    // Act
    Kokkos::View<double**, Kokkos::LayoutRight, Kokkos::HostSpace> m_values_csr_host;
    Kokkos::View<int*, Kokkos::LayoutRight, Kokkos::HostSpace> m_col_idx_csr_host;
    Kokkos::View<int*, Kokkos::LayoutRight, Kokkos::HostSpace> m_nnz_per_row_csr_host;
    // TODO: Make members of classes private!

    polar_matrix_assembler.
    operator()(get_const_field(coeff_alpha), get_const_field(coeff_beta), mapping, evaluator);
    // TODO: Should the polar matrix assembler really own the matrix data?
    auto& values = polar_matrix_assembler.m_values_csr_host;
    auto& row_ptrs_extended = polar_matrix_assembler.m_nnz_per_row_csr_host;
    auto& cols = polar_matrix_assembler.m_col_idx_csr_host;

    std::size_t n_rows = row_ptrs_extended.size() - 2;
    std::size_t nnz = values.size();

    std::unique_ptr<MatrixType> matrix
            = MatrixType::create(exec, gko::dim<2>(n_rows, n_rows), values.size());

    Kokkos::View<int*, Kokkos::DefaultHostExecutionSpace>
            row_ptrs_matrix {matrix->get_row_ptrs(), n_rows + 1};
    Kokkos::View<int*, Kokkos::DefaultHostExecutionSpace> cols_matrix {matrix->get_col_idxs(), nnz};
    Kokkos::View<double**, Kokkos::DefaultHostExecutionSpace>
            values_matrix {matrix->get_values(), 1, nnz};

    auto row_ptrs = Kokkos::subview(row_ptrs_extended, std::make_pair(0UL, n_rows + 1));

    Kokkos::deep_copy(row_ptrs_matrix, row_ptrs);
    Kokkos::deep_copy(cols_matrix, cols);
    Kokkos::deep_copy(values_matrix, values);

    // Assert
    std::size_t matrix_size
            = PolarBSplinesPatch1::singular_idx_range<PolarBSplinesPatch1>().size()
              + ddc::discrete_space<PolarBSplinesPatch1>().tensor_bspline_idx_range().size();

    ddc::parallel_fill(polar_spline, 5.0);

    std::unique_ptr<VectorType> spline_coeffs_vec
            = VectorType::create(exec, gko::dim<2>(matrix_size, 1));
    collect_spline_coeffs_in_vec<PolarBSplinesPatch1, Patch1>(
            spline_coeffs_vec.get(),
            spline_coeffs_polar_singular,
            spline_coeffs_polar_tensor);

    std::unique_ptr<VectorType> matrix_times_spline_coeffs
            = VectorType::create(exec, gko::dim<2>(matrix_size, 1));

    assert("Dimensions of spline_coeffs_vec and matrix do not match"
           && matrix->get_size()[1] == spline_coeffs_vec->get_size()[0]);
    assert("Dimensions of matrix and matrix_times_spline_coeffs do not match"
           && matrix->get_size()[0] == matrix_times_spline_coeffs->get_size()[0]);
    matrix->apply(spline_coeffs_vec, matrix_times_spline_coeffs);


    std::unique_ptr<VectorType> h1_norm = gko::initialize<VectorType>({1.0}, exec);
    spline_coeffs_vec->transpose()->apply(matrix_times_spline_coeffs, h1_norm);

    ASSERT_NEAR(h1_norm->at(0, 0), 0.0, 1e-12);
}


TEST_F(PoissonMatrixOnionTwoPatches, H1NormZero)
{
    // Act

    std::unique_ptr<MatrixType> poisson_matrix
            = get_matrix_assembler().
              operator()(get_const_field(coeff_alpha), get_const_field(coeff_beta), mapping);


    // Assert

    std::size_t matrix_size
            = PolarBSplinesPatch1::singular_idx_range<PolarBSplinesPatch1>().size()
              + ddc::discrete_space<PolarBSplinesPatch1>().tensor_bspline_idx_range().size()

    VectorType::create(exec, gko::dim<2>(matrix_size, 1));
    ddc::parallel_fill(spline_coeffs_polar_singular, 0.0);
    ddc::parallel_fill(spline_coeffs_polar_tensor, 0.0);

    // DUPLICATED LOGIC: Collect coefficients in vector and compute vector-matrix-vector product
    std::unique_ptr<VectorType> spline_coeffs_vec
            = VectorType::create(exec, gko::dim<2>(matrix_size, 1));
    collect_spline_coeffs_in_vec<PolarBSplines>(
            spline_coeffs_vec.get(),
            spline_coeffs_polar_singular,
            spline_coeffs_polar_tensor);

    std::unique_ptr<VectorType> broken_matrix_times_spline_coeffs
            = VectorType::create(exec, gko::dim<2>(matrix_size, 1));

    assert("Dimensions of spline_coeffs_vec and broken_matrix do not match"
           && broken_matrix->get_size()[1] == spline_coeffs_vec->get_size()[0]);
    assert("Dimensions of broken_matrix and broken_matrix_times_spline_coeffs do not match"
           && broken_matrix->get_size()[0] == broken_matrix_times_spline_coeffs->get_size()[0]);
    broken_matrix->apply(spline_coeffs_vec, broken_matrix_times_spline_coeffs);


    std::unique_ptr<VectorType> h1_norm = gko::initialize<VectorType>({1.0}, exec);
    spline_coeffs_vec->transpose()->apply(broken_matrix_times_spline_coeffs, h1_norm);

    ASSERT_NEAR(h1_norm->at(0, 0), 0.0, 1e-14);
}


TEST_F(PoissonMatrixOnionTwoPatches, H1NormConstant)
{
    // Act

    std::unique_ptr<MatrixType> broken_matrix
            = get_broken_matrix_assembler().
              operator()(get_const_field(coeff_alpha), get_const_field(coeff_beta), mapping);

    // Assert
    std::function<double(double)> field1_analytic
            = [](double r1) { return (r1 < 0.5 ? 0.0 : r1 - 0.5); };

    ddc::host_for_each(idx_range_rtheta1, [&](Patch1::Idx12 idx) {
        double const r1 = ddc::get<R<1>>(ddc::coordinate(idx));
        field1(idx) = field1_analytic(r1);
    });

    interpolate_field(spline_coeffs_with_reps, builder, field);

    std::size_t matrix_size = idx_range_spline_polar.size() + idx_range_bspline2.size();

    std::unique_ptr<VectorType> spline_coeffs_vec
            = VectorType::create(exec, gko::dim<2>(matrix_size, 1));
    collect_spline_coeffs_in_vec<PolarBSplinesPatch1, Patch1, Patch2>(
            spline_coeffs_vec.get(),
            spline_coeffs_polar_singular,
            spline_coeffs_polar_tensor,
            spline_coeffs2);

    std::unique_ptr<VectorType> broken_matrix_times_spline_coeffs
            = VectorType::create(exec, gko::dim<2>(matrix_size, 1));

    assert("Dimensions of spline_coeffs_vec and broken_matrix do not match"
           && broken_matrix->get_size()[1] == spline_coeffs_vec->get_size()[0]);
    assert("Dimensions of broken_matrix and broken_matrix_times_spline_coeffs do not match"
           && broken_matrix->get_size()[0] == broken_matrix_times_spline_coeffs->get_size()[0]);
    broken_matrix->apply(spline_coeffs_vec, broken_matrix_times_spline_coeffs);

    std::unique_ptr<VectorType> h1_norm = gko::initialize<VectorType>({1.0}, exec);
    spline_coeffs_vec->transpose()->apply(broken_matrix_times_spline_coeffs, h1_norm);

    ASSERT_NEAR(std::abs(h1_norm->at(0, 0)), 0.0, 1e-12);
}


TEST_F(PoissonMatrix, H1QuadraticInR)
{
    // Act

    std::unique_ptr<MatrixType> broken_matrix
            = get_broken_matrix_assembler().
              operator()(get_const_field(coeff_alpha), get_const_field(coeff_beta), mapping);

    // Assert
    std::function<double(double)> field1_analytic
            = [](double r1) { return (r1 < 0.5 ? 0.0 : r1 * r1 - 0.25); };

    ddc::host_for_each(idx_range_rtheta1, [&](Patch1::Idx12 idx) {
        double const r1 = ddc::get<R<1>>(ddc::coordinate(idx));
        field1(idx) = field1_analytic(r1);
    });

    interpolate_field(spline_coeffs_with_reps, builder, field);

    std::size_t matrix_size = idx_range_spline_polar.size() + idx_range_bspline2.size();

    std::unique_ptr<VectorType> spline_coeffs_vec
            = VectorType::create(exec, gko::dim<2>(matrix_size, 1));
    collect_spline_coeffs_in_vec<PolarBSplinesPatch1, Patch1, Patch2>(
            spline_coeffs_vec.get(),
            spline_coeffs_polar_singular,
            spline_coeffs_polar_tensor,
            spline_coeffs2);

    std::unique_ptr<VectorType> broken_matrix_times_spline_coeffs
            = VectorType::create(exec, gko::dim<2>(matrix_size, 1));

    assert("Dimensions of spline_coeffs_vec and broken_matrix do not match"
           && broken_matrix->get_size()[1] == spline_coeffs_vec->get_size()[0]);
    assert("Dimensions of broken_matrix and broken_matrix_times_spline_coeffs do not match"
           && broken_matrix->get_size()[0] == broken_matrix_times_spline_coeffs->get_size()[0]);
    broken_matrix->apply(spline_coeffs_vec, broken_matrix_times_spline_coeffs);

    std::unique_ptr<VectorType> h1_norm = gko::initialize<VectorType>({1.0}, exec);
    spline_coeffs_vec->transpose()->apply(broken_matrix_times_spline_coeffs, h1_norm);

    ASSERT_NEAR(std::abs(h1_norm->at(0, 0)), 0.0, 1e-12);
}

TEST_F(PoissonMatrix, H1NormHighDegreeR)
{
    // Act
    std::unique_ptr<MatrixType> broken_matrix
            = get_broken_matrix_assembler().
              operator()(get_const_field(coeff_alpha), get_const_field(coeff_beta), mapping);

    // Assert

    auto function_to_be_interpolated = [](double const r) {
        return r <= 0.5 ? -std::pow(1.5, 6) : std::pow(r - 0.5, 6) - std::pow(1.5, 6);
    };

    ddc::host_for_each(idx_range_rtheta1, [&](Patch1::IdxRange12::discrete_element_type idx) {
        auto coord = ddc::coordinate(idx);
        double const r1 = ddc::get<R<1>>(coord);
        field1(idx) = function_to_be_interpolated(r1);
    });
    ddc::host_for_each(idx_range_rtheta2, [&](Patch2::IdxRange12::discrete_element_type idx) {
        auto coord = ddc::coordinate(idx);
        double const r2 = ddc::get<R<2>>(coord);
        field2(idx) = function_to_be_interpolated(r2);
    });


    interpolate_fields<Patch1, Patch2>(spline_coeffs_with_reps, builders, fields);

    ddc::host_for_each(
            PolarBSplinesPatch1::singular_idx_range<PolarBSplinesPatch1>(),
            [&](Idx<PolarBSplinesPatch1> idx) { polar_spline(idx) = -std::pow(1.5, 6); });
    ddc::host_for_each(
            get_idx_range(spline_coeffs_polar_tensor),
            [&](Idx<Patch1::BSplines1, Patch1::BSplines2> idx) {
                spline_coeffs_polar_tensor(idx) = spline_coeffs_with_reps1(idx);
            });

    int const matrix_size = idx_range_spline_polar.size() + idx_range_bspline2.size();

    auto spline_coeffs_vec = VectorType::create(exec, gko::dim<2>(matrix_size, 1));
    ddc::parallel_deepcopy(spline_coeffs2, spline_coeffs_with_reps2[idx_range_bspline2]);
    collect_spline_coeffs_in_vec<PolarBSplinesPatch1, Patch1, Patch2>(
            spline_coeffs_vec.get(),
            spline_coeffs_polar_singular,
            spline_coeffs_polar_tensor,
            spline_coeffs2);

    auto stiffness_times_coeff = VectorType::create(exec, spline_coeffs_vec->get_size());
    broken_matrix->apply(spline_coeffs_vec, stiffness_times_coeff);
    auto norm_h1_broken = gko::initialize<VectorType>({-1.0}, exec);

    spline_coeffs_vec->transpose()->apply(stiffness_times_coeff, norm_h1_broken);
    double const norm_h1_ref_squared = 3334.9920704784536;
    double const norm_h1_broken_diff
            = std::sqrt(norm_h1_broken->at(0, 0)) - std::sqrt(norm_h1_ref_squared);
    ASSERT_NEAR(norm_h1_broken_diff, 1.25254e-4, 1e-7);
}


// TODO: Still lots of duplication in the tests
TEST_F(PoissonMatrixOnionTwoPatches, H1NormRThetaFunction)
{
    // Act
    std::unique_ptr<MatrixType> broken_matrix
            = get_broken_matrix_assembler().
              operator()(get_const_field(coeff_alpha), get_const_field(coeff_beta), mapping);

    // Assert
    auto field1_analytic = [](double r1, double theta1) {
        return (r1 < 0.5 ? 0.0 : (r1 - 0.5) * (r1 - 0.5) * std::sin(theta1));
    };
    auto field2_analytic
            = [](double r2, double theta2) { return (r2 - 0.5) * (r2 - 0.5) * std::sin(theta2); };

    ddc::parallel_for_each(
            idx_range_rtheta1,
            KOKKOS_LAMBDA(Patch1::Idx12 idx) {
                double const r1 = ddc::get<R<1>>(ddc::coordinate(idx));
                double const theta1 = ddc::get<Theta<1>>(ddc::coordinate(idx));
                field1(idx) = field1_analytic(r1, theta1);
            });
    ddc::parallel_for_each(
            idx_range_rtheta2,
            KOKKOS_LAMBDA(Patch2::Idx12 idx) {
                double const r2 = ddc::get<R<2>>(ddc::coordinate(idx));
                double const theta2 = ddc::get<Theta<2>>(ddc::coordinate(idx));
                field2(idx) = field2_analytic(r2, theta2);
            });

    interpolate_fields<Patch1, Patch2>(spline_coeffs_with_reps, builders, fields);
    std::size_t matrix_size = idx_range_spline_polar.size() + idx_range_bspline2.size();


    ddc::parallel_fill(polar_spline, 0.0);
    ddc::host_for_each(
            get_idx_range(spline_coeffs_polar_tensor),
            [&](Idx<Patch1::BSplines1, Patch1::BSplines2> idx) {
                spline_coeffs_polar_tensor(idx) = spline_coeffs_with_reps1(idx);
            });
    ddc::host_for_each(idx_range_bspline2, [&](Idx<Patch2::BSplines1, Patch2::BSplines2> idx) {
        spline_coeffs2(idx) = spline_coeffs_with_reps2(idx);
    });

    std::unique_ptr<VectorType> spline_coeffs_vec
            = VectorType::create(exec, gko::dim<2>(matrix_size, 1));
    collect_spline_coeffs_in_vec<PolarBSplinesPatch1, Patch1, Patch2>(
            spline_coeffs_vec.get(),
            spline_coeffs_polar_singular,
            spline_coeffs_polar_tensor,
            spline_coeffs2);

    std::unique_ptr<VectorType> broken_matrix_times_spline_coeffs
            = VectorType::create(exec, gko::dim<2>(matrix_size, 1));

    broken_matrix->apply(spline_coeffs_vec, broken_matrix_times_spline_coeffs);
    std::unique_ptr<VectorType> h1_norm_squared = gko::initialize<VectorType>({1.0}, exec);
    spline_coeffs_vec->transpose()->apply(broken_matrix_times_spline_coeffs, h1_norm_squared);

    double const h1_norm_squared_ref = 25.748551186492673;
    double const h1_norm_diff = std::abs(std::sqrt(h1_norm_squared->at(0,0)) - std::sqrt(h1_norm_squared_ref));
    ASSERT_NEAR(h1_norm_diff, 0.001215992, 1e-9);
}


TEST_F(PoissonMatrixOnionTwoPatches, MatrixComparisonFirstPatch)
{
    SplineInterpolatorPatch1 interpolator(idx_range_rtheta1);
    PolarSplineFEMPoissonLikeSolver<GridR<1>, GridTheta<1>, PolarBSplinesPatch1, SplineInterpolatorPatch1, Mapping>
            poisson_solver_polar_only(
                    mapping,
                    interpolator);

    DFieldMem<Patch1::IdxRange12> coeff_alpha_field_mem{idx_range_rtheta1};
    DField<Patch1::IdxRange12> coeff_alpha_field{coeff_alpha_field_mem};
    ddc::parallel_fill(coeff_alpha_field, 1.0);
    DFieldMem<Patch1::IdxRange12> coeff_beta_field_mem{idx_range_rtheta1};
    DField<Patch1::IdxRange12> coeff_beta_field{coeff_beta_field_mem};
    ddc::parallel_fill(coeff_beta_field, 0.0);
    poisson_solver_polar_only.update_coefficients(get_const_field(coeff_alpha_field), get_const_field(coeff_beta_field));

    auto [vals_polar, col_idxs_polar, row_ptrs_polar]
            = poisson_solver_polar_only.m_gko_matrix->get_batch_csr();


    // Act
    std::unique_ptr<MatrixType> broken_matrix
            = get_broken_matrix_assembler().
              operator()(get_const_field(coeff_alpha), get_const_field(coeff_beta), mapping);

    // Assert
    int max_error_row_ptr = 0;

    for (int i(0); i < 3; i++) {
        int row_ptr_broken = broken_matrix->get_row_ptrs()[i];
        int row_ptr_polar = row_ptrs_polar(i);
        int error = std::abs(row_ptr_broken - row_ptr_polar);
        max_error_row_ptr = (error > max_error_row_ptr ? error : max_error_row_ptr);
    }


    for (int r_idx(0); r_idx < 4; r_idx++) {
        for (int theta_idx(0); theta_idx < theta1_ncells; theta_idx++) {
            int row_idx = 3 + r_idx * theta1_ncells + theta_idx;
            int row_ptr_broken = broken_matrix->get_row_ptrs()[row_idx];
            int row_ptr_polar = row_ptrs_polar(row_idx);
            int error = std::abs(row_ptr_broken - row_ptr_polar);
            max_error_row_ptr = (error > max_error_row_ptr ? error : max_error_row_ptr);
        }
    }
    ASSERT_EQ(max_error_row_ptr, 0);

    int max_error_col_idx = 0;
    double max_error_value = 0.0;

    auto update_max_error_col_idx_value = [&](int row_idx) {
        for (int j(broken_matrix->get_row_ptrs()[row_idx]);
             j < broken_matrix->get_row_ptrs()[row_idx + 1];
             j++) {
            int col_idx_broken = broken_matrix->get_col_idxs()[j];
            int col_idx_poisson = col_idxs_polar(j);
            int error_col_idx = std::abs(col_idx_broken - col_idx_poisson);
            max_error_col_idx
                    = (error_col_idx > max_error_col_idx ? error_col_idx : max_error_col_idx);

            double value_broken = broken_matrix->get_values()[j];
            double value_polar = vals_polar(0, j);
            double error_val = std::abs(value_broken - value_polar);
            max_error_value = (error_val > max_error_value ? error_val : max_error_value);
        }
    };

    for (int i(0); i < 3; i++) {
        update_max_error_col_idx_value(i);
    }

    for (int r_idx(0); r_idx < 4; r_idx++) {
        for (int theta_idx(0); theta_idx < theta1_ncells; theta_idx++) {
            int row_idx = 3 + r_idx * theta1_ncells + theta_idx;
            update_max_error_col_idx_value(row_idx);
        }
    }

    ASSERT_EQ(max_error_col_idx, 0);
    ASSERT_NEAR(max_error_value, 0.0, 1.0e-14);
}