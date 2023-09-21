#pragma once

#include <sll/spline_boundary_value.hpp>

#include "sll/view.hpp"

/**
 * @brief A class for describing a spline boundary value by a constant extrapolation for 1D evaluator.
 *
 * To define the value of a function on B-splines out of the domain, we here  use a constant
 * extrapolation on the edge.
 *
 * @see SplineBoundaryValue
 */
template <class BSplines>
class ConstantExtrapolationBoundaryValue : public SplineBoundaryValue<BSplines>
{
public:
    /**
	 * @brief Indicate the dimension we are working on.
	 */
    using tag_type = typename BSplines::tag_type;
    /**
	 * @brief Indicate the coordinate type in the dimension of the boundary condition.
	 */
    using coord_type = ddc::Coordinate<tag_type>;

private:
    coord_type m_eval_pos;

public:
    /**
     * @brief Instantiate a ConstantExtrapolationBoundaryValue.
     *
     * The boundary value will be the same as at the coordinate eval_pos given.
     *
     * @param[in] eval_pos
     * 			Coordinate inside the domain where we will evaluate each points outside the domain.
     */
    explicit ConstantExtrapolationBoundaryValue(coord_type eval_pos) : m_eval_pos(eval_pos) {}

    ~ConstantExtrapolationBoundaryValue() override = default;

    /**
     * @brief Get the value of the function on B-splines at a coordinate outside the domain.
     *
     * @param[in] pos
     * 			The coordinate where we want to evaluate the function on B-splines.
     * @param[in] spline_coef
     *			The coefficients of the function on B-splines.
     *
     *@return A double with the value of the function on B-splines evaluated at the coordinate.
     */
    double operator()(
            coord_type const pos,
            ddc::ChunkSpan<double const, ddc::DiscreteDomain<BSplines>> const spline_coef)
            const final
    {
        std::array<double, BSplines::degree() + 1> values;
        DSpan1D const vals = as_span(values);

        ddc::DiscreteElement<BSplines> idx
                = ddc::discrete_space<BSplines>().eval_basis(vals, m_eval_pos);

        double y = 0.0;
        for (std::size_t i = 0; i < BSplines::degree() + 1; ++i) {
            y += spline_coef(idx + i) * vals(i);
        }
        return y;
    }
};



/**
 * @brief A class for describing a spline boundary value by a constant extrapolation for 2D evaluator.
 *
 * To define the value of a function on B-splines out of the domain, we here use a constant
 * extrapolation on the edge.
 *
 * @see SplineBoundaryValue
 */
template <class BSplines1, class BSplines2, class BCDim>
class ConstantExtrapolationBoundaryValue2D : public SplineBoundaryValue2D<BSplines1, BSplines2>
{
public:
    /**
	 * @brief Indicate the first dimension we are working on.
	 */
    using Dim1 = typename BSplines1::tag_type;
    /**
	 * @brief Indicate the second dimension we are working on.
	 */
    using Dim2 = typename BSplines2::tag_type;
    /**
	 * @brief Indicate coordinate type in the first dimension we are working on.
	 */
    using coord_type1 = ddc::Coordinate<Dim1>;
    /**
	 * @brief Indicate coordinate type in the second dimension we are working on.
	 */
    using coord_type2 = ddc::Coordinate<Dim2>;
    /**
	 * @brief Indicate coordinate type in the dimension of the boundary condition.
	 */
    using coord_type_bc = ddc::Coordinate<BCDim>;
    /**
     * @brief Indicate the complementary dimension of the boundary condition dimension.
     */
    using NoBCDim = typename std::conditional_t<std::is_same_v<Dim1, BCDim>, Dim2, Dim1>;
    /**
     * @brief Indicate the coordinate type of the complementary dimension of the boundary condition dimension.
     */
    using coord_type_no_bc = ddc::Coordinate<NoBCDim>;
    /**
     * @brief Boolean set at True if the Bsplines on the first dimension is periodic.
     */
    static constexpr bool BSplines1_is_periodic = BSplines1::is_periodic();
    /**
     * @brief Boolean set at True if the Bsplines on the second dimension is periodic.
     */
    static constexpr bool BSplines2_is_periodic = BSplines2::is_periodic();
    /**
     * @brief Boolean set at True if the Bsplines on the complementary dimension of the boundary
     * condition dimension is periodic.
     */
    static constexpr bool BSplinesNoBC_is_periodic
            = std::is_same_v<Dim1, BCDim> ? BSplines2_is_periodic : BSplines1_is_periodic;



    static_assert(std::is_same_v<BCDim, Dim1> || std::is_same_v<BCDim, Dim2>);
    static_assert(std::is_same_v<NoBCDim, Dim1> || std::is_same_v<NoBCDim, Dim2>);

private:
    coord_type_bc const m_eval_pos_bc;
    coord_type_no_bc const m_eval_pos_no_bc_min;
    coord_type_no_bc const m_eval_pos_no_bc_max;

public:
    /**
     * @brief Instantiate a ConstantExtrapolationBoundaryValue2D.
     *
     * The boundary value will be the same as at the coordinate given in a dimension given.
     * The dimension of the input defines the dimension of the boundary condition.
     * The second and the third parameters are needed in case of non-periodic splines on the
     * no-boundary condition dimension (the complementary dimension of the boundary condition),
     * because the evaluator can receive coordinates outside the domain in both dimension.
     *
     * @param[in] eval_pos_bc
     * 			Coordinate in the dimension given inside the domain where we will evaluate
     * 			each points outside the domain.
     * @param[in] eval_pos_no_bc_min
     * 			The minimum coordinate inside the domain on the complementary dimension of the boundary condition.
     * @param[in] eval_pos_no_bc_max
     * 			The maximum coordinate inside the domain on the complementary dimension of the boundary condition.
     */
    template <typename T = std::enable_if<!BSplinesNoBC_is_periodic>>
    explicit ConstantExtrapolationBoundaryValue2D(
            coord_type_bc const eval_pos_bc,
            coord_type_no_bc const eval_pos_no_bc_min,
            coord_type_no_bc const eval_pos_no_bc_max)
        : m_eval_pos_bc(eval_pos_bc)
        , m_eval_pos_no_bc_min(eval_pos_no_bc_min)
        , m_eval_pos_no_bc_max(eval_pos_no_bc_max) {};

    /**
     * @brief Instantiate a ConstantExtrapolationBoundaryValue2D with periodic splines
     * on the no-boundary condition dimension.
     *
     * The boundary value will be the same as at the coordinate given in a dimension given.
     * The dimension of the input defines the dimension of the boundary condition.
     * This constructor can only be used with periodic splines on the no-boundary condition
     * dimension. Otherwise, we have to use the previous constructor.
     *
     * @param[in] eval_pos_bc
     *          Coordinate in the dimension given inside the domain where we will evaluate
     *          each point outside the domain.
     */
    template <typename T = std::enable_if<BSplinesNoBC_is_periodic>>
    explicit ConstantExtrapolationBoundaryValue2D(coord_type_bc const eval_pos_bc)
        : m_eval_pos_bc(eval_pos_bc)
        , m_eval_pos_no_bc_min(coord_type_no_bc(0.))
        , m_eval_pos_no_bc_max(coord_type_no_bc(0.)) {};


    ~ConstantExtrapolationBoundaryValue2D() override = default;

    /**
     * @brief Get the value of the function on B-splines at a coordinate outside the domain.
     *
     * In the dimension defined in the constructor Dim1 (or Dim2), it sets the coordinate pos_1 (or pos_2)
     * given at the m_eval_pos_bc coordinate if it is outside the domain.
     * If the coordinate on the complementary dimension of the boundary condition dimension pos_2 (or pos_1) is
     * outside the domain, then it also sets the coordinate at eval_pos_no_bc_min
     * (if pos_2 (or pos_1) @f$ < @f$ eval_pos_no_bc_min) or
     * at eval_pos_no_bc_max (if pos_2 (or pos_1) @f$ > @f$ eval_pos_no_bc_max).
     *
     * @param[in] pos_1
     * 			The coordinate in the first dimension where we want to evaluate the function on B-splines
     * @param[in] pos_2
     * 			The coordinate in the second dimension where we want to evaluate the function on B-splines.
     * @param[in] spline_coef
     *			The coefficients of the function on B-splines.
     *
     *@return A double with the value of the function on B-splines evaluated at the coordinate.
     */
    double operator()(
            coord_type1 const pos_1,
            coord_type2 const pos_2,
            ddc::ChunkSpan<double const, ddc::DiscreteDomain<BSplines1, BSplines2>> const
                    spline_coef) const final
    {
        coord_type1 eval_pos_1;
        coord_type2 eval_pos_2;

        if constexpr (BSplinesNoBC_is_periodic) {
            if constexpr (std::is_same_v<BCDim, Dim1>) {
                eval_pos_1 = m_eval_pos_bc;
                eval_pos_2 = pos_2;
            } else {
                eval_pos_1 = pos_1;
                eval_pos_2 = m_eval_pos_bc;
            }

        } else {
            double bc_min = m_eval_pos_no_bc_min;
            double bc_max = m_eval_pos_no_bc_max;

            if constexpr (std::is_same_v<BCDim, Dim1>) {
                eval_pos_1 = m_eval_pos_bc;
                // Take the maximum between m_eval_pos_no_bc_min and pos_2
                // and the minimum between m_eval_pos_no_bc_max and pos_2
                // to be inside the domain.
                eval_pos_2 = coord_type2(std::max(bc_min, std::min(bc_max, double(pos_2))));
            } else {
                // Take the maximum between m_eval_pos_no_bc_min and pos_1
                // and the minimum between m_eval_pos_no_bc_max and pos_1
                // to be inside the domain.
                eval_pos_1 = coord_type1(std::max(bc_min, std::min(bc_max, double(pos_1))));
                eval_pos_2 = m_eval_pos_bc;
            }
        }

        std::array<double, BSplines1::degree() + 1> values1;
        DSpan1D const vals1 = as_span(values1);
        std::array<double, BSplines2::degree() + 1> values2;
        DSpan1D const vals2 = as_span(values2);

        ddc::DiscreteElement<BSplines1> idx1
                = ddc::discrete_space<BSplines1>().eval_basis(vals1, eval_pos_1);
        ddc::DiscreteElement<BSplines2> idx2
                = ddc::discrete_space<BSplines2>().eval_basis(vals2, eval_pos_2);

        double y = 0.0;
        for (std::size_t i = 0; i < BSplines1::degree() + 1; ++i) {
            for (std::size_t j = 0; j < BSplines2::degree() + 1; ++j) {
                y += spline_coef(idx1 + i, idx2 + j) * vals1(i) * vals2(j);
            }
        }
        return y;
    }
};
