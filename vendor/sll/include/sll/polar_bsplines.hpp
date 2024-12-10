#pragma once
#include <array>
#include <vector>

#include <ddc/ddc.hpp>

#include <sll/bernstein.hpp>
#include <sll/mapping/cartesian_to_barycentric.hpp>
#include <sll/mapping/discrete_to_cartesian.hpp>
#include <sll/polar_spline.hpp>
#include <sll/view.hpp>

/**
 * A class containing all information describing polar bsplines.
 *
 * Polar bsplines are 2D bsplines with a special treatment for the central singular point
 * of a polar index range. At this singular point new bsplines are created which traverse the
 * singular point and ensure the desired continuity condition.
 *
 * @tparam BSplinesR  The basis of radial bsplines from which the polar bsplines are constructed.
 * @tparam BSplinesTheta  The poloidal bspline from which the polar bsplines are constructed.
 * @tparam C          The continuity condition. The resulting splines will be continuously
 *                    differentiable C times. If C == -1 then the resulting spline representation
 *                    will be discontinuous at the singular point.
 */
template <class BSplinesR, class BSplinesTheta, int C>
class PolarBSplines
{
    static_assert(C >= -1, "Parameter `C` cannot be less than -1");
    static_assert(C < 2, "Values larger than 1 are not implemented for parameter `C`");
    static_assert(!BSplinesR::is_periodic(), "Radial bsplines must not be periodic.");
    static_assert(!BSplinesR::is_uniform(), "Radial bsplines must have knots at the boundary.");
    static_assert(BSplinesTheta::is_periodic(), "Poloidal bsplines should be periodic.");

private:
    // Tags to determine what to evaluate
    struct eval_type
    {
    };

    struct eval_deriv_type
    {
    };

public:
    /// The radial bspline from which the polar bsplines are constructed.
    using BSplinesR_tag = BSplinesR;

    /// The poloidal bspline from which the polar bsplines are constructed.
    using BSplinesTheta_tag = BSplinesTheta;


    /// The tag for the radial direction of the bsplines.
    using DimR = typename BSplinesR::continuous_dimension_type;

    /// The tag for the poloidal direction of the bsplines.
    using DimTheta = typename BSplinesTheta::continuous_dimension_type;

public:
    /// The continuity enforced by the bsplines at the singular point.
    static int constexpr continuity = C;

public:
    /**
     * The tag denoting the discrete dimension described by this class.
     *
     * This is the tag which should be used to create a Chunk whose contents are each associated with a PolarBSpline.
     * In other words a spline defined on this basis would have the type:
     * ddc::Chunk<double, ddc::DiscreteDomain<PolarBSplines>;
     */
    using discrete_dimension_type = PolarBSplines;

    /**
     * The type of a 2D index for the subset of the polar bsplines which can be expressed as a tensor
     * product of 1D bsplines.
     */
    using tensor_product_index_type = ddc::DiscreteElement<BSplinesR, BSplinesTheta>;

    /**
     * The type of the 2D idx_range for the subset of the polar bsplines which can be expressed as a tensor
     * product of 1D bsplines.
     */
    using tensor_product_idx_range_type = ddc::DiscreteDomain<BSplinesR, BSplinesTheta>;

    /**
     * The type of a 2D vector for the subset of the polar bsplines which can be expressed as a tensor
     * product of 1D bsplines.
     */
    using tensor_product_idx_step_type = ddc::DiscreteVector<BSplinesR, BSplinesTheta>;

private:
    using IdxR = ddc::DiscreteElement<BSplinesR>;
    using IdxTheta = ddc::DiscreteElement<BSplinesTheta>;
    using IdxStepR = ddc::DiscreteVector<BSplinesR>;
    using IdxStepTheta = ddc::DiscreteVector<BSplinesTheta>;

public:
    /**
     * Get the number of singular bsplines i.e. bsplines which traverse the singular point.
     *
     * @returns The number of bsplines which traverse the singular point.
     */
    static constexpr std::size_t n_singular_basis()
    {
        return (C + 1) * (C + 2) / 2;
    }


    /**
     * @brief Get the ddc::DiscreteDomain containing the indices of the b-splines which traverse
     * the singular point.
     *
     * @returns The ddc::DiscreteDomain containing the indices of the b-splines which traverse
     * the singular point.
     */
    template <class DDim>
    static constexpr ddc::DiscreteDomain<DDim> singular_idx_range()
    {
        return ddc::DiscreteDomain<DDim>(
                ddc::DiscreteElement<DDim> {0},
                ddc::DiscreteVector<DDim> {n_singular_basis()});
    }

    /**
     * Get the index of the polar bspline which, when evaluated at the same point, returns the
     * same values as the 2D tensor product bspline indicated by the index passed as an argument.
     *
     * @param idx The index of a 2D BSpline which is expressed as a tensor product of 1D BSplines.
     *
     * @returns The index of the basis spline in the PolarBSpline index range.
     */
    template <class DDim>
    static KOKKOS_FUNCTION ddc::DiscreteElement<DDim> get_polar_index(
            tensor_product_index_type const& idx)
    {
        int const r_idx = ddc::select<BSplinesR>(idx).uid();
        int const theta_idx = ddc::select<BSplinesTheta>(idx).uid();
        assert(r_idx >= C + 1);
        int local_idx((r_idx - C - 1) * ddc::discrete_space<BSplinesTheta>().nbasis() + theta_idx);
        return ddc::DiscreteElement<DDim>(n_singular_basis() + local_idx);
    }

    /**
     * Get the 2D index of the tensor product bspline which, when evaluated at the same point,
     * returns the same values as the polar bspline indicated by the index passed as an argument.
     *
     * @param idx The index of the basis spline in the PolarBSpline index range.
     *
     * @returns The index of the equivalent 2D BSpline expressed as a 2D tensor product of 1D BSplines.
     */
    template <class DDim>
    static KOKKOS_FUNCTION tensor_product_index_type
    get_2d_index(ddc::DiscreteElement<DDim> const& idx)
    {
        assert(idx.uid() >= n_singular_basis());
        int const idx_2d = idx.uid() - n_singular_basis();
        int const r_idx = idx_2d / ddc::discrete_space<BSplinesTheta>().nbasis();
        int const theta_idx = idx_2d - r_idx * ddc::discrete_space<BSplinesTheta>().nbasis();
        ddc::DiscreteElement<BSplinesR> r_idx_elem(r_idx + C + 1);
        ddc::DiscreteElement<BSplinesTheta> theta_idx_elem(theta_idx);
        return ddc::DiscreteElement<BSplinesR, BSplinesTheta>(r_idx_elem, theta_idx_elem);
    }

public:
    /**
     * The Impl class holds the implementation of the PolarBSplines. The implementation is specific to the
     * memory space so that the Chunks can be defined with index ranges related to instances of this class.
     *
     * @tparam MemorySpace Indicates where the object is saved. This is either on the host or the device.
     */
    template <class DDim, class MemorySpace>
    class Impl
    {
        template <class ODDim, class OMemorySpace>
        friend class Impl;

        template <class ExecSpace, class PBSpl, class OMemorySpace>
        friend PolarSplineSpan<PBSpl, OMemorySpace> integrals(
                ExecSpace const& execution_space,
                PolarSplineSpan<PBSpl, OMemorySpace> int_vals);

    private:
        /**
         * The type of the index range for the linear combinations defining the bsplines which traverse the singular point.
         *
         * The bsplines which traverse the singular O-point are constructed from a linear combination of 2D
         * bsplines. These 2D bsplines can be expressed as a tensor product of 1D bsplines. This type
         * describes the index range on which the coefficients of these linear combinations are defined. There is
         * an index for the polar bspline being constructed, and 2 indices for the 2D bspline.
         */
        using singular_basis_linear_combination_idx_range_type
                = ddc::DiscreteDomain<DDim, BSplinesR, BSplinesTheta>;

        ddc::Chunk<
                double,
                singular_basis_linear_combination_idx_range_type,
                ddc::KokkosAllocator<double, MemorySpace>>
                m_singular_basis_elements_alloc;

        ddc::ChunkSpan<
                double,
                singular_basis_linear_combination_idx_range_type,
                std::experimental::layout_right,
                MemorySpace>
                m_singular_basis_elements;

    public:
        /// The tag for the first corner of the Barycentric coordinates
        struct Corner1Tag
        {
        };
        /// The tag for the second corner of the Barycentric coordinates
        struct Corner2Tag
        {
        };
        /// The tag for the third corner of the Barycentric coordinates
        struct Corner3Tag
        {
        };

        template <class DiscreteMapping>
        struct IntermediateBernsteinBasis
            : TriangularBernsteinPolynomialBasis<
                      typename DiscreteMapping::cartesian_tag_x,
                      typename DiscreteMapping::cartesian_tag_y,
                      Corner1Tag,
                      Corner2Tag,
                      Corner3Tag,
                      C>
        {
        };

        /// The tag which should be used to create a Chunk whose contents are each associated with a PolarBSpline.
        using discrete_dimension_type = PolarBSplines;

        /// The type of an index associated with a PolarBSpline.
        using discrete_element_type = ddc::DiscreteElement<DDim>;

        /// The type of a index range of PolarBSplines.
        using discrete_domain_type = ddc::DiscreteDomain<DDim>;

        /// The type of a vector associated with a PolarBSpline.
        using discrete_vector_type = ddc::DiscreteVector<DDim>;

        /**
         * A constructor for the PolarBSplines.
         *
         * @param curvilinear_to_cartesian  A mapping from curvilinear to cartesian coordinates. This is used to find the
         *                                  singular point and determine the Barycentric coordinates which are used to define
         *                                  the new basis splines which cross the singular point.
         */
        template <class DiscreteMapping>
        Impl(const DiscreteMapping& curvilinear_to_cartesian)
        {
            static_assert(std::is_same_v<MemorySpace, Kokkos::HostSpace>);
            using DimX = typename DiscreteMapping::cartesian_tag_x;
            using DimY = typename DiscreteMapping::cartesian_tag_y;
            using mapping_tensor_product_index_type = ddc::DiscreteElement<
                    typename DiscreteMapping::BSplineR,
                    typename DiscreteMapping::BSplineTheta>;
            if constexpr (C > -1) {
                const ddc::Coordinate<DimX, DimY> pole
                        = curvilinear_to_cartesian(ddc::Coordinate<DimR, DimTheta>(0.0, 0.0));
                const double x0 = ddc::get<DimX>(pole);
                const double y0 = ddc::get<DimY>(pole);
                double tau = 0.0;
                for (std::size_t i(0); i < ddc::discrete_space<BSplinesTheta>().size(); ++i) {
                    const ddc::Coordinate<DimX, DimY> point
                            = curvilinear_to_cartesian.control_point(
                                    mapping_tensor_product_index_type(1, i));

                    const double c_x = ddc::get<DimX>(point);
                    const double c_y = ddc::get<DimY>(point);

                    double tau1 = -2.0 * (c_x - x0);
                    double tau2 = c_x - x0 - sqrt(3.0) * (c_y - y0);
                    double tau3 = c_x - x0 + sqrt(3.0) * (c_y - y0);
                    tau = tau > tau1 ? tau : tau1;
                    tau = tau > tau2 ? tau : tau2;
                    tau = tau > tau3 ? tau : tau3;
                }
                // Determine the corners for the barycentric coordinates
                const ddc::Coordinate<DimX, DimY> corner1(x0 + tau, y0);
                const ddc::Coordinate<DimX, DimY>
                        corner2(x0 - 0.5 * tau, y0 + 0.5 * tau * sqrt(3.0));
                const ddc::Coordinate<DimX, DimY>
                        corner3(x0 - 0.5 * tau, y0 - 0.5 * tau * sqrt(3.0));

                const CartesianToBarycentric<DimX, DimY, Corner1Tag, Corner2Tag, Corner3Tag>
                        barycentric_coordinate_converter(corner1, corner2, corner3);

                using BernsteinBasis = IntermediateBernsteinBasis<DiscreteMapping>;

                ddc::init_discrete_space<BernsteinBasis>(barycentric_coordinate_converter);

                // The number of radial bases used to construct the bsplines traversing the singular point.
                constexpr IdxStepR nr_in_singular(C + 1);
                assert(nr_in_singular.value() < int(ddc::discrete_space<BSplinesR>().size()));

                // The number of poloidal bases used to construct the bsplines traversing the singular point.
                const IdxStepTheta np_in_singular(ddc::discrete_space<BSplinesTheta>().nbasis());

                // The number of elements of the poloidal basis which will have an associated coefficient
                // (This will be larger than np_in_singular as it includes the periodicity)
                const IdxStepTheta np_tot(ddc::discrete_space<BSplinesTheta>().size());

                // The index range of the 2D bsplines in the innermost circles from which the polar bsplines
                // traversing the singular point will be constructed.
                tensor_product_idx_range_type const dom_bsplines_inner(
                        tensor_product_index_type(0, 0),
                        tensor_product_idx_step_type(nr_in_singular, np_tot));

                // Initialise memory
                m_singular_basis_elements_alloc
                        = ddc::Chunk<double, singular_basis_linear_combination_idx_range_type>(
                                singular_basis_linear_combination_idx_range_type(
                                        singular_idx_range<DDim>(),
                                        dom_bsplines_inner));
                m_singular_basis_elements = m_singular_basis_elements_alloc.span_view();

                ddc::DiscreteDomain<BernsteinBasis> bernstein_idx_range(
                        ddc::DiscreteElement<BernsteinBasis> {0},
                        ddc::DiscreteVector<BernsteinBasis> {n_singular_basis()});

                ddc::DiscreteDomain<BSplinesTheta> poloidal_spline_idx_range
                        = ddc::discrete_space<BSplinesTheta>().full_domain();

                for (IdxR const ir : ddc::DiscreteDomain<BSplinesR>(IdxR(0), IdxStepR(C + 1))) {
                    for (IdxTheta const ip : poloidal_spline_idx_range.take_first(np_in_singular)) {
                        const ddc::Coordinate<DimX, DimY> point
                                = curvilinear_to_cartesian.control_point(
                                        mapping_tensor_product_index_type(ir, ip));
                        ddc::Chunk<double, ddc::DiscreteDomain<BernsteinBasis>> bernstein_vals(
                                bernstein_idx_range);
                        ddc::discrete_space<BernsteinBasis>().eval_basis(bernstein_vals, point);
                        // Fill spline coefficients
                        for (auto k : bernstein_idx_range) {
                            m_singular_basis_elements(discrete_element_type {k.uid()}, ir, ip)
                                    = bernstein_vals(k);
                        }
                    }
                    for (discrete_element_type k : singular_idx_range<DDim>()) {
                        for (IdxTheta const ip : poloidal_spline_idx_range.take_first(
                                     IdxStepTheta {BSplinesTheta::degree()})) {
                            m_singular_basis_elements(k, ir, ip + np_in_singular)
                                    = m_singular_basis_elements(k, ir, ip);
                        }
                    }
                }
            } else {
                // Initialise m_singular_basis_elements to avoid any problems in the copy constructor
                tensor_product_idx_range_type const empty_dom_bsplines(
                        tensor_product_index_type(0, 0),
                        tensor_product_idx_step_type(0, 0));
                m_singular_basis_elements_alloc
                        = ddc::Chunk<double, singular_basis_linear_combination_idx_range_type>(
                                singular_basis_linear_combination_idx_range_type(
                                        singular_idx_range<DDim>(),
                                        empty_dom_bsplines));
                m_singular_basis_elements = m_singular_basis_elements_alloc.span_view();
            }
        }

        /**
         * A copy constructor for the PolarBSplines.
         *
         * @param impl The PolarBSplines being copied.
         */
        template <class OriginMemorySpace>
        explicit Impl(Impl<DDim, OriginMemorySpace> const& impl)
            : m_singular_basis_elements_alloc(impl.m_singular_basis_elements.domain())
        {
            m_singular_basis_elements = m_singular_basis_elements_alloc.span_view();
            ddc::parallel_deepcopy(m_singular_basis_elements, impl.m_singular_basis_elements);
        }

        /**
         * A copy constructor for the PolarBSplines.
         *
         * @param x The PolarBSplines being copied.
         */
        Impl(Impl const& x) = default;

        /**
         * A copy constructor for the PolarBSplines taking a temporary r-value.
         *
         * @param x The PolarBSplines being copied.
         */
        Impl(Impl&& x) = default;

        /**
         * The destructor for the PolarBSplines.
         */
        ~Impl() = default;

        /**
         * A copy operator for the PolarBSplines.
         *
         * @param x The PolarBSplines being copied.
         *
         * @returns A reference to this PolarBSpline.
         */
        Impl& operator=(Impl const& x) = default;

        /**
         * A copy operator for the PolarBSplines taking a temporary r-value.
         *
         * @param x The PolarBSplines being copied.
         *
         * @returns A reference to this PolarBSpline.
         */
        Impl& operator=(Impl&& x) = default;

        /**
         * @brief Evaluate the polar basis splines at the coordinate p.
         *
         * Evaluate all the b-spline elements near the singular point which cannot be
         * expressed as a tensor product of 1D bsplines, as well as the non-zero b-spline
         * elements which can be expressed as a tensor product of 1D bsplines.
         *
         * @param[out] singular_values  The value of the b-spline elements near the singular point
         *                              which cannot be expressed as a tensor product of 1D bsplines,
         *                              evaluated at the coordinate p.
         * @param[out] values           The value of the non-zero b-spline elements which can be
         *                              expressed as a tensor product of 1D bsplines.
         * @param[in] p                 The coordinate where the basis functions are evaluated.
         *
         * @returns The 2D tensor product index of the first b-spline element in the values array.
         */
        tensor_product_index_type eval_basis(
                DSpan1D singular_values,
                DSpan2D values,
                ddc::Coordinate<DimR, DimTheta> p) const;

        /**
         * @brief Evaluate the radial derivative of the polar basis splines at the coordinate p.
         *
         * Evaluate the radial derivative of all the b-spline elements near the singular point which
         * cannot be expressed as a tensor product of 1D bsplines, as well as the non-zero b-spline
         * elements which can be expressed as a tensor product of 1D bsplines.
         *
         * @param[out] singular_derivs  The value of the radial derivative b-spline elements near the
         *                              singular point which cannot be expressed as a tensor product
         *                              of 1D bsplines, evaluated at the coordinate p.
         * @param[out] derivs           The value of the radial derivative of the non-zero b-spline
         *                              elements which can be expressed as a tensor product of 1D bsplines.
         * @param[in] p                 The coordinate where the basis functions are evaluated.
         *
         * @returns The 2D tensor product index of the first b-spline element in the values array.
         */
        tensor_product_index_type eval_deriv_r(
                DSpan1D singular_derivs,
                DSpan2D derivs,
                ddc::Coordinate<DimR, DimTheta> p) const;

        /**
         * @brief Evaluate the poloidal derivative of the polar basis splines at the coordinate p.
         *
         * Evaluate the poloidal derivative of all the b-spline elements near the singular point which
         * cannot be expressed as a tensor product of 1D bsplines, as well as the non-zero b-spline
         * elements which can be expressed as a tensor product of 1D bsplines.
         *
         * @param[out] singular_derivs  The value of the poloidal derivative b-spline elements near the
         *                              singular point which cannot be expressed as a tensor product
         *                              of 1D bsplines, evaluated at the coordinate p.
         * @param[out] derivs           The value of the poloidal derivative of the non-zero b-spline
         *                              elements which can be expressed as a tensor product of 1D bsplines.
         * @param[in] p                 The coordinate where the basis functions are evaluated.
         *
         * @returns The 2D tensor product index of the first b-spline element in the values array.
         */
        tensor_product_index_type eval_deriv_theta(
                DSpan1D singular_derivs,
                DSpan2D derivs,
                ddc::Coordinate<DimR, DimTheta> p) const;

        /**
         * @brief Evaluate the second order derivative of the polar basis splines in the radial and poloidal
         * directions, at the coordinate p.
         *
         * Evaluate the 2nd order derivative of all the b-spline elements near the singular point which
         * cannot be expressed as a tensor product of 1D bsplines, as well as the non-zero b-spline
         * elements which can be expressed as a tensor product of 1D bsplines.
         *
         * @param[out] singular_derivs  The value of the 2nd order derivative b-spline elements near the
         *                              singular point which cannot be expressed as a tensor product
         *                              of 1D bsplines, evaluated at the coordinate p.
         * @param[out] derivs           The value of the 2nd order derivative of the non-zero b-spline
         *                              elements which can be expressed as a tensor product of 1D bsplines.
         * @param[in] p                 The coordinate where the basis functions are evaluated.
         *
         * @returns The 2D tensor product index of the first b-spline element in the values array.
         */
        tensor_product_index_type eval_deriv_r_and_theta(
                DSpan1D singular_derivs,
                DSpan2D derivs,
                ddc::Coordinate<DimR, DimTheta> p) const;

        /**
         * Calculate the integrals of each of the basis splines.
         *
         * @param[out] int_vals The integrals of the basis splines.
         */
        template <class MemorySpace2>
        [[deprecated("Use `integrals` instead")]] void integrals(
                PolarSplineSpan<DDim, MemorySpace2> int_vals) const;

        /**
         * Get the total number of basis functions.
         *
         * @returns The number of basis functions.
         */
        std::size_t nbasis() const noexcept
        {
            std::size_t nr = ddc::discrete_space<BSplinesR>().nbasis() - C - 1;
            std::size_t ntheta = ddc::discrete_space<BSplinesTheta>().nbasis();
            return n_singular_basis() + nr * ntheta;
        }

        /**
         * Returns the index range containing the indices of all the polar b-splines.
         *
         * @returns The index range containing the indices of all the polar b-splines.
         */
        discrete_domain_type full_domain() const noexcept
        {
            return discrete_domain_type(discrete_element_type {0}, discrete_vector_type {nbasis()});
        }

        /**
         * @brief Returns the ddc::DiscreteDomain containing the indices of the b-splines which don't
         * traverse the singular point and can be expressed as a tensor-product of 1D b-splines.
         *
         * @returns The ddc::DiscreteDomain containing the indices of the b-splines which don't traverse
         * the singular point.
         */
        discrete_domain_type tensor_bspline_idx_range() const noexcept
        {
            return full_domain().remove_first(discrete_vector_type {n_singular_basis()});
        }

    private:
        template <class EvalTypeR, class EvalTypeTheta>
        ddc::DiscreteElement<BSplinesR, BSplinesTheta> eval(
                DSpan1D singular_values,
                DSpan2D values,
                ddc::Coordinate<DimR, DimTheta> coord_eval,
                EvalTypeR const,
                EvalTypeTheta const) const;
    };
};

template <class BSplinesR, class BSplinesTheta, int C>
template <class DDim, class MemorySpace>
ddc::DiscreteElement<BSplinesR, BSplinesTheta> PolarBSplines<BSplinesR, BSplinesTheta, C>::
        Impl<DDim, MemorySpace>::eval_basis(
                DSpan1D singular_values,
                DSpan2D values,
                ddc::Coordinate<DimR, DimTheta> p) const
{
    return eval(singular_values, values, p, eval_type(), eval_type());
}

template <class BSplinesR, class BSplinesTheta, int C>
template <class DDim, class MemorySpace>
ddc::DiscreteElement<BSplinesR, BSplinesTheta> PolarBSplines<BSplinesR, BSplinesTheta, C>::
        Impl<DDim, MemorySpace>::eval_deriv_r(
                DSpan1D singular_derivs,
                DSpan2D derivs,
                ddc::Coordinate<DimR, DimTheta> p) const
{
    return eval(singular_derivs, derivs, p, eval_deriv_type(), eval_type());
}

template <class BSplinesR, class BSplinesTheta, int C>
template <class DDim, class MemorySpace>
ddc::DiscreteElement<BSplinesR, BSplinesTheta> PolarBSplines<BSplinesR, BSplinesTheta, C>::
        Impl<DDim, MemorySpace>::eval_deriv_theta(
                DSpan1D singular_derivs,
                DSpan2D derivs,
                ddc::Coordinate<DimR, DimTheta> p) const
{
    return eval(singular_derivs, derivs, p, eval_type(), eval_deriv_type());
}

template <class BSplinesR, class BSplinesTheta, int C>
template <class DDim, class MemorySpace>
ddc::DiscreteElement<BSplinesR, BSplinesTheta> PolarBSplines<BSplinesR, BSplinesTheta, C>::
        Impl<DDim, MemorySpace>::eval_deriv_r_and_theta(
                DSpan1D singular_derivs,
                DSpan2D derivs,
                ddc::Coordinate<DimR, DimTheta> p) const
{
    return eval(singular_derivs, derivs, p, eval_deriv_type(), eval_deriv_type());
}

template <class BSplinesR, class BSplinesTheta, int C>
template <class DDim, class MemorySpace>
template <class EvalTypeR, class EvalTypeTheta>
ddc::DiscreteElement<BSplinesR, BSplinesTheta> PolarBSplines<BSplinesR, BSplinesTheta, C>::
        Impl<DDim, MemorySpace>::eval(
                DSpan1D singular_values,
                DSpan2D values,
                ddc::Coordinate<DimR, DimTheta> coord_eval,
                EvalTypeR const,
                EvalTypeTheta const) const
{
    assert(singular_values.extent(0) == n_singular_basis());
    assert(values.extent(0) == BSplinesR::degree() + 1);
    assert(values.extent(1) == BSplinesTheta::degree() + 1);
    static_assert(
            std::is_same_v<EvalTypeR, eval_type> || std::is_same_v<EvalTypeR, eval_deriv_type>);
    static_assert(
            std::is_same_v<
                    EvalTypeTheta,
                    eval_type> || std::is_same_v<EvalTypeTheta, eval_deriv_type>);

    ddc::DiscreteElement<BSplinesR> jmin_r;
    ddc::DiscreteElement<BSplinesTheta> jmin_theta;

    std::size_t constexpr nr = BSplinesR::degree() + 1;
    std::size_t constexpr ntheta = BSplinesTheta::degree() + 1;

    std::array<double, nr> vals_r_ptr;
    std::array<double, ntheta> vals_theta_ptr;
    DSpan1D const vals_r(vals_r_ptr.data(), nr);
    DSpan1D const vals_theta(vals_theta_ptr.data(), ntheta);

    if constexpr (std::is_same_v<EvalTypeR, eval_type>) {
        jmin_r = ddc::discrete_space<BSplinesR>().eval_basis(vals_r, ddc::select<DimR>(coord_eval));
    } else if constexpr (std::is_same_v<EvalTypeR, eval_deriv_type>) {
        jmin_r = ddc::discrete_space<BSplinesR>().eval_deriv(vals_r, ddc::select<DimR>(coord_eval));
    }
    if constexpr (std::is_same_v<EvalTypeTheta, eval_type>) {
        jmin_theta = ddc::discrete_space<BSplinesTheta>()
                             .eval_basis(vals_theta, ddc::select<DimTheta>(coord_eval));
    } else if constexpr (std::is_same_v<EvalTypeTheta, eval_deriv_type>) {
        jmin_theta = ddc::discrete_space<BSplinesTheta>()
                             .eval_deriv(vals_theta, ddc::select<DimTheta>(coord_eval));
    }

    std::size_t nr_done = 0;

    if (jmin_r.uid() < C + 1) {
        nr_done = C + 1 - jmin_r.uid();
        for (discrete_element_type k : singular_idx_range<DDim>()) {
            singular_values(k.uid()) = 0.0;
            for (std::size_t i(0); i < nr_done; ++i) {
                for (std::size_t j(0); j < ntheta; ++j) {
                    singular_values(k.uid())
                            += m_singular_basis_elements(k, jmin_r + i, jmin_theta + j) * vals_r[i]
                               * vals_theta[j];
                }
            }
        }
    } else {
        for (std::size_t k(0); k < n_singular_basis(); ++k) {
            singular_values(k) = 0.0;
        }
    }

    for (std::size_t i(0); i < nr - nr_done; ++i) {
        for (std::size_t j(0); j < ntheta; ++j) {
            values(i, j) = vals_r[i + nr_done] * vals_theta[j];
        }
    }
    for (std::size_t i(nr - nr_done); i < nr; ++i) {
        for (std::size_t j(0); j < ntheta; ++j) {
            values(i, j) = 0.0;
        }
    }
    return ddc::DiscreteElement<BSplinesR, BSplinesTheta>(jmin_r, jmin_theta);
}

template <class ExecSpace, class DDim, class MemorySpace>
PolarSplineSpan<DDim, MemorySpace> integrals(
        ExecSpace const& execution_space,
        PolarSplineSpan<DDim, MemorySpace> int_vals)
{
    static_assert(
            Kokkos::SpaceAccessibility<ExecSpace, MemorySpace>::accessible,
            "MemorySpace has to be accessible for ExecutionSpace.");
    using BSplinesR = typename DDim::BSplinesR_tag;
    using BSplinesTheta = typename DDim::BSplinesTheta_tag;
    using tensor_product_idx_range_type = ddc::DiscreteDomain<BSplinesR, BSplinesTheta>;
    using tensor_product_idx_type = ddc::DiscreteElement<BSplinesR, BSplinesTheta>;
    using IdxR = ddc::DiscreteElement<BSplinesR>;
    using IdxTheta = ddc::DiscreteElement<BSplinesTheta>;

    auto r_bspl_space = ddc::discrete_space<BSplinesR>();
    auto theta_bspl_space = ddc::discrete_space<BSplinesTheta>();

    assert(int_vals.singular_spline_coef.domain().extents() == DDim::n_singular_basis());
    assert(int_vals.spline_coef.domain().front().template uid<BSplinesR>() == DDim::continuity + 1);
    assert(int_vals.spline_coef.domain().back().template uid<BSplinesR>()
           == r_bspl_space.nbasis() - 1);
    assert(int_vals.spline_coef.domain().template extent<BSplinesTheta>()
                   == theta_bspl_space.nbasis()
           || int_vals.spline_coef.domain().template extent<BSplinesTheta>()
                      == theta_bspl_space.size());

    ddc::Chunk<double, ddc::DiscreteDomain<BSplinesR>, ddc::KokkosAllocator<double, MemorySpace>>
            r_integrals_alloc(r_bspl_space.full_domain().take_first(
                    ddc::DiscreteVector<BSplinesR> {r_bspl_space.nbasis()}));
    ddc::Chunk<
            double,
            ddc::DiscreteDomain<BSplinesTheta>,
            ddc::KokkosAllocator<double, MemorySpace>>
            theta_integrals_alloc(theta_bspl_space.full_domain().take_first(
                    ddc::DiscreteVector<BSplinesTheta> {theta_bspl_space.size()}));
    ddc::ChunkSpan r_integrals = r_integrals_alloc.span_view();
    ddc::ChunkSpan theta_integrals = theta_integrals_alloc.span_view();

    ddc::integrals(execution_space, r_integrals);
    ddc::integrals(execution_space, theta_integrals);

    ddc::DiscreteDomain<BSplinesR, BSplinesTheta> singular_2d_idx_range(
            ddc::discrete_space<DDim>().m_singular_basis_elements.domain());
    ddc::ChunkSpan singular_spline_integrals = int_vals.singular_spline_coef.span_view();

    ddc::DiscreteDomain<DDim> singular_idx_range = DDim::template singular_idx_range<DDim>();
    Kokkos::parallel_for(
            Kokkos::TeamPolicy<>(execution_space, singular_idx_range.size(), Kokkos::AUTO),
            KOKKOS_LAMBDA(const Kokkos::TeamPolicy<>::member_type& team) {
                const int idx = team.league_rank();
                ddc::DiscreteElement<DDim> k(idx);

                // Sum over quadrature dimensions
                double teamSum = 0;
                Kokkos::parallel_reduce(
                        Kokkos::TeamThreadMDRange(
                                team,
                                singular_2d_idx_range.template extent<BSplinesR>().value(),
                                singular_2d_idx_range.template extent<BSplinesTheta>().value()),
                        [&](int r_thread_index, int theta_thread_index, double& sum) {
                            IdxR i(r_thread_index);
                            IdxTheta j(theta_thread_index);
                            sum += ddc::discrete_space<DDim>().m_singular_basis_elements(k, i, j)
                                   * r_integrals(i) * theta_integrals(j);
                        },
                        teamSum);
                singular_spline_integrals(k) = teamSum;
            });


    ddc::DiscreteDomain<BSplinesR> r_tensor_product_dom(
            ddc::select<BSplinesR>(int_vals.spline_coef.domain()));
    tensor_product_idx_range_type
            tensor_bspline_idx_range(r_tensor_product_dom, theta_integrals.domain());
    ddc::ChunkSpan spline_integrals = int_vals.spline_coef.span_view();

    ddc::parallel_for_each(
            execution_space,
            tensor_bspline_idx_range,
            KOKKOS_LAMBDA(tensor_product_idx_type idx) {
                int_vals.spline_coef(idx) = r_integrals(ddc::select<BSplinesR>(idx))
                                            * theta_integrals(ddc::select<BSplinesTheta>(idx));
            });

    if (int_vals.spline_coef.domain().template extent<BSplinesTheta>() == theta_bspl_space.size()) {
        ddc::DiscreteDomain<BSplinesTheta> periodic_points(theta_integrals.domain().take_last(
                ddc::DiscreteVector<BSplinesTheta> {BSplinesTheta::degree()}));
        tensor_product_idx_range_type repeat_idx_range(r_tensor_product_dom, periodic_points);
        ddc::parallel_fill(execution_space, int_vals.spline_coef, 0.0);
    }
    return int_vals;
}

template <class BSplinesR, class BSplinesTheta, int C>
template <class DDim, class MemorySpace>
template <class MemorySpace2>
void PolarBSplines<BSplinesR, BSplinesTheta, C>::Impl<DDim, MemorySpace>::integrals(
        PolarSplineSpan<DDim, MemorySpace2> int_vals) const
{
    integrals(Kokkos::DefaultHostExecutionSpace(), int_vals);
}
