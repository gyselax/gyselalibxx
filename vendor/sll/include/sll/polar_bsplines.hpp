#pragma once
#include <array>
#include <vector>

#include <ddc/ddc.hpp>

#include <sll/bernstein.hpp>
#include <sll/bspline.hpp>
#include <sll/mapping/barycentric_coordinates.hpp>
#include <sll/mapping/discrete_mapping_to_cartesian.hpp>
#include <sll/null_boundary_value.hpp>
#include <sll/spline_builder.hpp>
#include <sll/spline_builder_2d.hpp>
#include <sll/spline_evaluator_2d.hpp>
#include <sll/view.hpp>

template <class BSplinesR, class BSplinesP, int C>
class PolarBSplines
{
    static_assert(C >= -1, "Parameter `C` cannot be less than -1");
    static_assert(C < 2, "Values larger than 1 are not implemented for parameter `C`");
    static_assert(!BSplinesR::is_periodic(), "");
    static_assert(!BSplinesR::is_uniform(), "Radial bsplines must have knots at the boundary");
    static_assert(BSplinesP::is_periodic(), "");

private:
    // Tags to determine what to evaluate
    struct eval_type
    {
    };

    struct eval_deriv_type
    {
    };

public:
    using BSplinesR_tag = BSplinesR;
    using BSplinesP_tag = BSplinesP;
    using DimR = typename BSplinesR::tag_type;
    using DimP = typename BSplinesP::tag_type;

public:
    static int constexpr continuity = C;

public:
    using discrete_dimension_type = PolarBSplines;

    using discrete_element_type = ddc::DiscreteElement<PolarBSplines>;

    using discrete_domain_type = ddc::DiscreteDomain<PolarBSplines>;

    using discrete_vector_type = ddc::DiscreteVector<PolarBSplines>;

public:
    // The number of bsplines describing the singular point
    static constexpr std::size_t n_singular_basis()
    {
        return (C + 1) * (C + 2) / 2;
    }

    static constexpr discrete_domain_type singular_domain()
    {
        return discrete_domain_type(
                discrete_element_type {0},
                discrete_vector_type {n_singular_basis()});
    }

    static discrete_element_type get_polar_index(
            ddc::DiscreteElement<BSplinesR, BSplinesP> const& idx)
    {
        int const r_idx = ddc::select<BSplinesR>(idx).uid();
        int const p_idx = ddc::select<BSplinesP>(idx).uid();
        assert(r_idx >= C + 1);
        int local_idx((r_idx - C - 1) * ddc::discrete_space<BSplinesP>().nbasis() + p_idx);
        return discrete_element_type(n_singular_basis() + local_idx);
    }

    static ddc::DiscreteElement<BSplinesR, BSplinesP> get_2d_index(discrete_element_type const& idx)
    {
        int const idx_2d = idx.uid() - n_singular_basis();
        int const r_idx = idx_2d / ddc::discrete_space<BSplinesP>().nbasis();
        int const p_idx = idx_2d - r_idx * ddc::discrete_space<BSplinesP>().nbasis();
        ddc::DiscreteElement<BSplinesR> r_idx_elem(r_idx + C + 1);
        ddc::DiscreteElement<BSplinesP> p_idx_elem(p_idx);
        return ddc::DiscreteElement<BSplinesR, BSplinesP>(r_idx_elem, p_idx_elem);
    }

private:
    using Spline2D = ddc::Chunk<double, ddc::DiscreteDomain<BSplinesR, BSplinesP>>;

public:
    template <class MemorySpace>
    class Impl
    {
        template <class OMemorySpace>
        friend class Impl;

    private:
        std::array<Spline2D, n_singular_basis()> m_singular_basis_elements;

        SplineEvaluator2D<BSplinesR, BSplinesP> m_spline_evaluator;

    public:
        struct Corner1Tag
        {
        };
        struct Corner2Tag
        {
        };
        struct Corner3Tag
        {
        };

        using discrete_dimension_type = PolarBSplines;

        template <class DimX, class DimY, class SplineBuilderR, class SplineBuilderP>
        Impl(const DiscreteToCartesian<DimX, DimY, SplineBuilder2D<SplineBuilderR, SplineBuilderP>>&
                     curvilinear_to_cartesian,
             SplineBuilderR const& spline_builder_r,
             SplineBuilderP const& spline_builder_p)
            : m_spline_evaluator(
                    g_null_boundary_2d<BSplinesR, BSplinesP>,
                    g_null_boundary_2d<BSplinesR, BSplinesP>,
                    g_null_boundary_2d<BSplinesR, BSplinesP>,
                    g_null_boundary_2d<BSplinesR, BSplinesP>)
        {
            if constexpr (C > -1) {
                const ddc::Coordinate<DimX, DimY> pole
                        = curvilinear_to_cartesian(ddc::Coordinate<DimR, DimP>(0.0, 0.0));
                const double x0 = ddc::get<DimX>(pole);
                const double y0 = ddc::get<DimY>(pole);
                double tau = 0.0;
                for (std::size_t i(0); i < ddc::discrete_space<BSplinesP>().size(); ++i) {
                    const ddc::Coordinate<DimX, DimY> point
                            = curvilinear_to_cartesian.control_point(
                                    ddc::DiscreteElement<BSplinesR, BSplinesP>(1, i));

                    const double c_x = ddc::get<DimX>(point);
                    const double c_y = ddc::get<DimY>(point);

                    double tau1 = -2.0 * (c_x - x0);
                    double tau2 = c_x - x0 - sqrt(3.0) * (c_y - y0);
                    double tau3 = c_x - x0 + sqrt(3.0) * (c_y - y0);
                    tau = tau > tau1 ? tau : tau1;
                    tau = tau > tau2 ? tau : tau2;
                    tau = tau > tau3 ? tau : tau3;
                }
                const ddc::Coordinate<DimX, DimY> corner1(x0 + tau, y0);
                const ddc::Coordinate<DimX, DimY>
                        corner2(x0 - 0.5 * tau, y0 + 0.5 * tau * sqrt(3.0));
                const ddc::Coordinate<DimX, DimY>
                        corner3(x0 - 0.5 * tau, y0 - 0.5 * tau * sqrt(3.0));

                const CartesianToBarycentricCoordinates<
                        DimX,
                        DimY,
                        Corner1Tag,
                        Corner2Tag,
                        Corner3Tag>
                        barycentric_coordinate_converter(corner1, corner2, corner3);
                using BernsteinBasis = BernsteinPolynomialBasis<
                        DimX,
                        DimY,
                        Corner1Tag,
                        Corner2Tag,
                        Corner3Tag,
                        C>;
                ddc::init_discrete_space<BernsteinBasis>(barycentric_coordinate_converter);

                using IndexR = ddc::DiscreteElement<BSplinesR>;
                using IndexP = ddc::DiscreteElement<BSplinesP>;
                using LengthR = ddc::DiscreteVector<BSplinesR>;
                using LengthP = ddc::DiscreteVector<BSplinesP>;
                using SplIndexRP = ddc::DiscreteElement<BSplinesR, BSplinesP>;

                constexpr LengthR nr(C + 1);
                const LengthP np(ddc::discrete_space<BSplinesP>().nbasis());
                const LengthP np_tot(ddc::discrete_space<BSplinesP>().size());
                assert(nr.value() < int(ddc::discrete_space<BSplinesR>().size()));

                ddc::DiscreteDomain<BSplinesR, BSplinesP> const dom_bsplines_inner(
                        ddc::DiscreteElement<BSplinesR, BSplinesP>(0, 0),
                        ddc::DiscreteVector<BSplinesR, BSplinesP>(nr, np_tot));

                for (std::size_t k(0); k < n_singular_basis(); ++k) {
                    // Initialise memory
                    m_singular_basis_elements[k] = Spline2D(dom_bsplines_inner);
                }

                ddc::DiscreteDomain<BernsteinBasis> bernstein_domain(
                        ddc::DiscreteElement<BernsteinBasis> {0},
                        ddc::DiscreteVector<BernsteinBasis> {n_singular_basis()});

                for (IndexR const ir : ddc::DiscreteDomain<BSplinesR>(IndexR(0), LengthR(C + 1))) {
                    for (IndexP const ip : spline_builder_p.spline_domain().take_first(np)) {
                        const ddc::Coordinate<DimX, DimY> point
                                = curvilinear_to_cartesian.control_point(SplIndexRP(ir, ip));
                        ddc::Chunk<double, ddc::DiscreteDomain<BernsteinBasis>> bernstein_vals(
                                bernstein_domain);
                        ddc::discrete_space<BernsteinBasis>().eval_basis(bernstein_vals, point);
                        // Fill spline coefficients
                        for (auto k : bernstein_domain) {
                            SplIndexRP const idx(ir.uid(), ip.uid());
                            m_singular_basis_elements[k.uid()](idx) = bernstein_vals(k);
                        }
                    }
                    for (std::size_t k(0); k < n_singular_basis(); ++k) {
                        for (std::size_t ip(0); ip < BSplinesP::degree(); ++ip) {
                            SplIndexRP const start_idx(ir.uid(), ip);
                            SplIndexRP const
                                    end_idx(ir.uid(),
                                            ddc::discrete_space<BSplinesP>().nbasis() + ip);
                            m_singular_basis_elements[k](end_idx)
                                    = m_singular_basis_elements[k](start_idx);
                        }
                    }
                }
            }
        }

        template <class OriginMemorySpace>
        explicit Impl(Impl<OriginMemorySpace> const& impl)
            : m_spline_evaluator(
                    g_null_boundary_2d<BSplinesR, BSplinesP>,
                    g_null_boundary_2d<BSplinesR, BSplinesP>,
                    g_null_boundary_2d<BSplinesR, BSplinesP>,
                    g_null_boundary_2d<BSplinesR, BSplinesP>)
        {
            for (std::size_t i = 0; i < m_singular_basis_elements.size(); ++i) {
                m_singular_basis_elements[i] = Spline2D(impl.m_singular_basis_elements[i].domain());
                ddc::parallel_deepcopy(
                        m_singular_basis_elements[i],
                        impl.m_singular_basis_elements[i]);
            }
        }

        Impl(Impl const& x) = default;

        Impl(Impl&& x) = default;

        ~Impl() = default;

        Impl& operator=(Impl const& x) = default;

        Impl& operator=(Impl&& x) = default;

        ddc::DiscreteElement<BSplinesR, BSplinesP> eval_basis(
                DSpan1D singular_values,
                DSpan2D values,
                ddc::Coordinate<DimR, DimP> p) const;
        ddc::DiscreteElement<BSplinesR, BSplinesP> eval_deriv_r(
                DSpan1D singular_derivs,
                DSpan2D derivs,
                ddc::Coordinate<DimR, DimP> p) const;
        ddc::DiscreteElement<BSplinesR, BSplinesP> eval_deriv_p(
                DSpan1D singular_derivs,
                DSpan2D derivs,
                ddc::Coordinate<DimR, DimP> p) const;
        ddc::DiscreteElement<BSplinesR, BSplinesP> eval_deriv_r_and_p(
                DSpan1D singular_derivs,
                DSpan2D derivs,
                ddc::Coordinate<DimR, DimP> p) const;
        void integrals(DSpan1D singular_int_vals, DSpan2D int_vals) const;

        std::size_t nbasis() const noexcept
        {
            std::size_t nr = ddc::discrete_space<BSplinesR>().nbasis() - C - 1;
            std::size_t np = ddc::discrete_space<BSplinesP>().nbasis();
            return n_singular_basis() + nr * np;
        }

        discrete_domain_type full_domain() const noexcept
        {
            return discrete_domain_type(discrete_element_type {0}, discrete_vector_type {nbasis()});
        }

        /**
         * @brief Returns the ddc::DiscreteDomain containing the indices of the non-singular b-splines.
         */
        discrete_domain_type non_singular_domain() const noexcept
        {
            return full_domain().remove_first(discrete_vector_type {n_singular_basis()});
        }

    private:
        template <class EvalTypeR, class EvalTypeP>
        ddc::DiscreteElement<BSplinesR, BSplinesP> eval(
                DSpan1D singular_values,
                DSpan2D values,
                ddc::Coordinate<DimR, DimP> coord_eval,
                EvalTypeR const,
                EvalTypeP const) const;
    };
};

template <class BSplinesR, class BSplinesP, int C>
template <class MemorySpace>
ddc::DiscreteElement<BSplinesR, BSplinesP> PolarBSplines<BSplinesR, BSplinesP, C>::
        Impl<MemorySpace>::eval_basis(
                DSpan1D singular_values,
                DSpan2D values,
                ddc::Coordinate<DimR, DimP> p) const
{
    return eval(singular_values, values, p, eval_type(), eval_type());
}

template <class BSplinesR, class BSplinesP, int C>
template <class MemorySpace>
ddc::DiscreteElement<BSplinesR, BSplinesP> PolarBSplines<BSplinesR, BSplinesP, C>::
        Impl<MemorySpace>::eval_deriv_r(
                DSpan1D singular_derivs,
                DSpan2D derivs,
                ddc::Coordinate<DimR, DimP> p) const
{
    return eval(singular_derivs, derivs, p, eval_deriv_type(), eval_type());
}

template <class BSplinesR, class BSplinesP, int C>
template <class MemorySpace>
ddc::DiscreteElement<BSplinesR, BSplinesP> PolarBSplines<BSplinesR, BSplinesP, C>::
        Impl<MemorySpace>::eval_deriv_p(
                DSpan1D singular_derivs,
                DSpan2D derivs,
                ddc::Coordinate<DimR, DimP> p) const
{
    return eval(singular_derivs, derivs, p, eval_type(), eval_deriv_type());
}

template <class BSplinesR, class BSplinesP, int C>
template <class MemorySpace>
ddc::DiscreteElement<BSplinesR, BSplinesP> PolarBSplines<BSplinesR, BSplinesP, C>::
        Impl<MemorySpace>::eval_deriv_r_and_p(
                DSpan1D singular_derivs,
                DSpan2D derivs,
                ddc::Coordinate<DimR, DimP> p) const
{
    return eval(singular_derivs, derivs, p, eval_deriv_type(), eval_deriv_type());
}

template <class BSplinesR, class BSplinesP, int C>
template <class MemorySpace>
template <class EvalTypeR, class EvalTypeP>
ddc::DiscreteElement<BSplinesR, BSplinesP> PolarBSplines<BSplinesR, BSplinesP, C>::
        Impl<MemorySpace>::eval(
                DSpan1D singular_values,
                DSpan2D values,
                ddc::Coordinate<DimR, DimP> coord_eval,
                EvalTypeR const,
                EvalTypeP const) const
{
    assert(singular_values.extent(0) == n_singular_basis());
    assert(values.extent(0) == BSplinesR::degree() + 1);
    assert(values.extent(1) == BSplinesP::degree() + 1);
    ddc::DiscreteElement<BSplinesR> jmin_r;
    ddc::DiscreteElement<BSplinesP> jmin_p;
    double vals_r_data[BSplinesR::degree() + 1];
    double vals_p_data[BSplinesP::degree() + 1];
    DSpan1D vals_r(vals_r_data, BSplinesR::degree() + 1);
    DSpan1D vals_p(vals_p_data, BSplinesP::degree() + 1);
    static_assert(
            std::is_same_v<EvalTypeR, eval_type> || std::is_same_v<EvalTypeR, eval_deriv_type>);
    static_assert(
            std::is_same_v<EvalTypeP, eval_type> || std::is_same_v<EvalTypeP, eval_deriv_type>);

    if constexpr (std::is_same_v<EvalTypeR, eval_type>) {
        jmin_r = ddc::discrete_space<BSplinesR>().eval_basis(vals_r, ddc::select<DimR>(coord_eval));
    } else if constexpr (std::is_same_v<EvalTypeR, eval_deriv_type>) {
        jmin_r = ddc::discrete_space<BSplinesR>().eval_deriv(vals_r, ddc::select<DimR>(coord_eval));
    }
    if constexpr (std::is_same_v<EvalTypeP, eval_type>) {
        jmin_p = ddc::discrete_space<BSplinesP>().eval_basis(vals_p, ddc::select<DimP>(coord_eval));
    } else if constexpr (std::is_same_v<EvalTypeP, eval_deriv_type>) {
        jmin_p = ddc::discrete_space<BSplinesP>().eval_deriv(vals_p, ddc::select<DimP>(coord_eval));
    }
    std::size_t nr = vals_r.size();
    std::size_t np = vals_p.size();
    std::size_t nr_done = 0;

    if (jmin_r.uid() < C + 1) {
        nr_done = C + 1 - jmin_r.uid();
        std::size_t np_eval = BSplinesP::degree() + 1;
        for (std::size_t k(0); k < n_singular_basis(); ++k) {
            singular_values(k) = 0.0;
            for (std::size_t i(0); i < nr_done; ++i) {
                for (std::size_t j(0); j < np_eval; ++j) {
                    ddc::DiscreteElement<BSplinesR, BSplinesP> icoeff(jmin_r + i, jmin_p + j);
                    singular_values(k)
                            += m_singular_basis_elements[k](icoeff) * vals_r(i) * vals_p(j);
                }
            }
        }
    } else {
        for (std::size_t k(0); k < n_singular_basis(); ++k) {
            singular_values(k) = 0.0;
        }
    }

    for (std::size_t i(0); i < nr - nr_done; ++i) {
        for (std::size_t j(0); j < np; ++j) {
            values(i, j) = vals_r(i + nr_done) * vals_p(j);
        }
    }
    for (std::size_t i(nr - nr_done); i < nr; ++i) {
        for (std::size_t j(0); j < np; ++j) {
            values(i, j) = 0.0;
        }
    }
    return ddc::DiscreteElement<BSplinesR, BSplinesP>(jmin_r, jmin_p);
}

template <class BSplinesR, class BSplinesP, int C>
template <class MemorySpace>
void PolarBSplines<BSplinesR, BSplinesP, C>::Impl<MemorySpace>::integrals(
        DSpan1D singular_int_vals,
        DSpan2D int_vals) const
{
    const int nr = ddc::discrete_space<BSplinesR>().ncells() + BSplinesR::degree() - C - 1;
    const int np = ddc::discrete_space<BSplinesP>().ncells() + BSplinesP::degree();
    assert(singular_int_vals.extent(0) == n_singular_basis());
    assert(int_vals.extent(0) == nr);
    assert(int_vals.extent(1) == np
           || int_vals.extent(1) == ddc::discrete_space<BSplinesP>().ncells());

    std::vector<double> r_integrals_data(nr);
    std::vector<double> p_integrals_data(ddc::discrete_space<BSplinesP>().ncells());
    DSpan1D r_integrals(r_integrals_data.data(), nr);
    DSpan1D p_integrals(p_integrals_data.data(), ddc::discrete_space<BSplinesP>().ncells());
    ddc::discrete_space<BSplinesR>().integrals(r_integrals);
    ddc::discrete_space<BSplinesP>().integrals(p_integrals);

    for (int k(0); k < n_singular_basis(); ++k) {
        singular_int_vals(k) = 0.0;
        ddc::for_each(
                m_singular_basis_elements[k].domain(),
                [=](ddc::DiscreteElement<BSplinesR, BSplinesP> const i) {
                    singular_int_vals(k) += m_singular_basis_elements[k](i)
                                            * r_integrals(ddc::select<BSplinesR>(i))
                                            * p_integrals(ddc::select<BSplinesP>(i));
                });
    }
    for (int i(n_singular_basis()); i < nr; ++i) {
        for (int j(0); j < ddc::discrete_space<BSplinesP>().ncells(); ++j) {
            int_vals(i, j) = r_integrals(i) * p_integrals(j);
        }
    }
    if (int_vals.extent(1) == np) {
        for (int i(n_singular_basis()); i < nr; ++i) {
            for (int j(0); j < BSplinesP::degree(); ++j) {
                int_vals(i, j) = int_vals(i, j + ddc::discrete_space<BSplinesP>().ncells());
            }
        }
    }
}
