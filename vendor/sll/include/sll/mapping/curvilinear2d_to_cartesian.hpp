#pragma once

#include <array>

#include <ddc/ddc.hpp>

template <class DimX, class DimY, class DimR, class DimP>
class Curvilinear2DToCartesian
{
public:
    using cartesian_tag_x = DimX;
    using cartesian_tag_y = DimY;
    using curvilinear_tag_r = DimR;
    using curvilinear_tag_p = DimP;
    using Matrix_2x2 = std::array<std::array<double, 2>, 2>;

public:
    Curvilinear2DToCartesian() = default;

    Curvilinear2DToCartesian(Curvilinear2DToCartesian const& other) = default;

    Curvilinear2DToCartesian(Curvilinear2DToCartesian&& x) = default;

    ~Curvilinear2DToCartesian() = default;

    Curvilinear2DToCartesian& operator=(Curvilinear2DToCartesian const& x) = default;

    Curvilinear2DToCartesian& operator=(Curvilinear2DToCartesian&& x) = default;

    virtual ddc::Coordinate<DimX, DimY> operator()(
            ddc::Coordinate<DimR, DimP> const& coord) const = 0;

    virtual double jacobian(ddc::Coordinate<DimR, DimP> const& coord) const
    {
        const double j_rr = jacobian_11(coord);
        const double j_rp = jacobian_12(coord);
        const double j_pr = jacobian_21(coord);
        const double j_pp = jacobian_22(coord);
        return j_rr * j_pp - j_rp * j_pr;
    }

    virtual double jacobian_11(ddc::Coordinate<DimR, DimP> const& coord) const = 0;
    virtual double jacobian_12(ddc::Coordinate<DimR, DimP> const& coord) const = 0;
    virtual double jacobian_21(ddc::Coordinate<DimR, DimP> const& coord) const = 0;
    virtual double jacobian_22(ddc::Coordinate<DimR, DimP> const& coord) const = 0;

    virtual void metric_tensor(ddc::Coordinate<DimR, DimP> const& coord, Matrix_2x2& matrix) const
    {
        const double J_rr = jacobian_11(coord);
        const double J_rp = jacobian_12(coord);
        const double J_pr = jacobian_21(coord);
        const double J_pp = jacobian_22(coord);
        matrix[0][0] = (J_rr * J_rr + J_pr * J_pr);
        matrix[0][1] = (J_rr * J_rp + J_pr * J_pp);
        matrix[1][0] = (J_rr * J_rp + J_pr * J_pp);
        matrix[1][1] = (J_rp * J_rp + J_pp * J_pp);
    }

    virtual void inverse_metric_tensor(ddc::Coordinate<DimR, DimP> const& coord, Matrix_2x2& matrix)
            const
    {
        const double J_rr = jacobian_11(coord);
        const double J_rp = jacobian_12(coord);
        const double J_pr = jacobian_21(coord);
        const double J_pp = jacobian_22(coord);
        const double jacob_2 = jacobian(coord) * jacobian(coord);
        matrix[0][0] = (J_rp * J_rp + J_pp * J_pp) / jacob_2;
        matrix[0][1] = (-J_rr * J_rp - J_pr * J_pp) / jacob_2;
        matrix[1][0] = (-J_rr * J_rp - J_pr * J_pp) / jacob_2;
        matrix[1][1] = (J_rr * J_rr + J_pr * J_pr) / jacob_2;
    }

    std::array<double, 2> to_covariant(
            std::array<double, 2> const& contravariant_vector,
            ddc::Coordinate<DimR, DimP> const& coord) const
    {
        Matrix_2x2 inv_metric_tensor;
        inverse_metric_tensor(coord, inv_metric_tensor);
        std::array<double, 2> covariant_vector;
        covariant_vector[0] = inv_metric_tensor[0][0] * contravariant_vector[0]
                              + inv_metric_tensor[0][1] * contravariant_vector[1];
        covariant_vector[1] = inv_metric_tensor[1][0] * contravariant_vector[0]
                              + inv_metric_tensor[1][1] * contravariant_vector[1];
        return covariant_vector;
    }
};
