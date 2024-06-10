// SPDX-License-Identifier: MIT

#pragma once

#include <ddc/ddc.hpp>

#include <ddc_helper.hpp>

#include "interface.hpp"

/**
 * @brief Transform a coordinate from one edge to a coordinate on the other edge.
 * 
 * According to the orientation stored in the interface, we compute 
 * * if True, 
 * 
 *  @f$ t \mapsto min_2 + \frac{t - min_1}{max_1 - min_1}(max_2 - min_2) @f$
 * 
 * * if False, 
 * 
 *  @f$ t \mapsto max_2 - \frac{t - min_1}{max_1 - min_1}(max_2 - min_2) @f$
 * 
 * where @f$ min_i @f$ and @f$ max_i @f$ are the minimum and maximum 
 * coordinates of the edge @f$ i @f$. 
 * 
 * @tparam IDim1
 *          The discrete dimension of the first edge. 
 * @tparam IDim2
 *          The discrete dimension of the second edge. 
 * 
*/
template <class IDim1, class IDim2>
class EdgeCoordinatesTransformation
{
    using RDim1 = typename IDim1::continuous_dimension_type;
    using RDim2 = typename IDim2::continuous_dimension_type;

    Interface<IDim1, IDim2> const m_interface;

public:
    /**
     * @brief Instantiate an EdgeCoordinatesTransformation.
     * 
     * @param interface 
     *      An Interface object between the two patches. 
    */
    EdgeCoordinatesTransformation(Interface<IDim1, IDim2> const interface)
        : m_interface(interface) {};
    ~EdgeCoordinatesTransformation() = default;


    /**
     * @brief Transform a coordinate on the edge in the dimension 
     * of the current patch to the analogous coordinate on the target patch. 
     * 
     * @param current_coord
     *      A coordinate on the edge of the current patch.
     * 
     * @tparam CurrentRDim 
     *      The current continuous dimension of the given coordinate coord. 
     * 
     * @return The analogous coordinate on the target patch. 
    */
    template <class CurrentRDim>
    ddc::Coordinate<std::conditional_t<std::is_same_v<CurrentRDim, RDim1>, RDim2, RDim1>>
    operator()(ddc::Coordinate<CurrentRDim> current_coord) const
    {
        // The discrete dimension of CurrentRDim
        using CurrentIDim = std::conditional_t<std::is_same_v<CurrentRDim, RDim1>, IDim1, IDim2>;
        // The other continuous dimension
        using ORDim = std::conditional_t<std::is_same_v<CurrentRDim, RDim1>, RDim2, RDim1>;
        // The other discrete dimension
        using OIDim = std::conditional_t<std::is_same_v<CurrentIDim, IDim1>, IDim2, IDim1>;


        bool const orientations_agree = m_interface.orientations_agree;

        ddc::Coordinate<CurrentRDim> current_min;
        double current_length;
        get_min_and_length<CurrentIDim>(current_min, current_length);

        ddc::Coordinate<ORDim> target_min;
        double target_length;
        get_min_and_length<OIDim>(target_min, target_length);

        double rescale_x = (current_coord - current_min) / current_length * target_length;

        if (!orientations_agree) {
            rescale_x = target_length - rescale_x;
        }
        return target_min + rescale_x; // This has type ddc::Coordinate<ODim>.
    };


private:
    /**
     * @brief Get the minimum coordinate of a patch on the edge and its length.
     * 
     * @tparam IDim 
     *      The discrete dimension of the edge of the current patch. 
    */
    template <class IDim>
    void get_min_and_length(
            ddc::Coordinate<typename IDim::continuous_dimension_type>& min,
            double& length) const
    {
        Edge<IDim> const edge = m_interface.template get_edge<IDim>();

        ddc::DiscreteDomain<IDim> const dom = edge.domain;

        min = ddc::coordinate(dom.front());
        length = ddcHelper::total_interval_length(dom);
    };
};
