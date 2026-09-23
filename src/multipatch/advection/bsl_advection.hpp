// SPDX-License-Identifier: MIT

#pragma once

#include "ddc_alias_inline_functions.hpp"
#include "ddc_aliases.hpp"
#include "directional_tag.hpp"
#include "multipatch_field.hpp"
#include "multipatch_type.hpp"
#include "types.hpp"
#include "vector_field.hpp"


/// @brief Advection operator with a backward semi-Lagrangian method for multi-patch domains.
/// @tparam MultipatchInterpolator Multi-patch interpolator operator.
/// @tparam FootFinder Operator computing the feet of the characteristics of the advection.
///                    It has to be adapated for multi-patch domain.
/// @tparam X First continuous dimension of the physical domain.
/// @tparam Y Second continuous dimension of the physical domain.
/// @tparam ...Patches  Tags for each patch of the multi-patch domain.
template <class MultipatchInterpolator, class FootFinder, class X, class Y, class... Patches>
class MultipatchBslAdvection
{
public:
    template <class Patch>
    using DConstVectorField_XY_OnPatch_host
            = host_t<VectorConstField<double, typename Patch::IdxRange12, NDTag<X, Y>>>;

private:
    MultipatchInterpolator& m_interpolator;
    FootFinder& m_find_feet;


public:
    /// @brief Instantiate the MultipatchBslAdvection.
    /// @param interpolator Reference to a interpolator on multi-patch domain.
    /// @param foot_finder Reference to an operator computing the feet of
    /// the characteristics of the advection on a multi-patch domain.
    MultipatchBslAdvection(MultipatchInterpolator& interpolator, FootFinder& foot_finder)
        : m_interpolator(interpolator)
        , m_find_feet(foot_finder)
    {
    }

    ~MultipatchBslAdvection() = default;

    /// @brief Advect the givem function on dt along the given advection field. 
    /// @param[in,out] allfdistribus Given function we want to advect.
    /// @param[in] advection_fields_xy Given advection field corresponding to a vector along the physical axes X and Y. 
    /// @param[in] dt Time step. 
    /// @return Advected function given in allfdistribus. 
    MultipatchField<DFieldOnPatch_host, Patches...> operator()(
            MultipatchField<DFieldOnPatch_host, Patches...> allfdistribus,
            MultipatchField<DConstVectorField_XY_OnPatch_host, Patches...> advection_fields_xy,
            double const dt) const
    {
        // Initialise the feet
        MultipatchFieldMem<CoordFieldMemOnPatch_host, Patches...> feet_alloc(
                advection_fields_xy.idx_range());
        MultipatchField<CoordFieldOnPatch_host, Patches...> feet(feet_alloc);
        intialise_the_feet_on_grids(feet);

        // Compute the characteristic feet at tn ----------------------------------------------------
        m_find_feet(feet, advection_fields_xy, dt);

        // Interpolate the function on the characteristic feet. -------------------------------------
        MultipatchField<CoordConstFieldOnPatch, Patches...> const_feet
                = get_const_field(MultipatchField);

        m_interpolator(allfdistribus, const_feet);

        return allfdistribus;
    }

private:
    template <class... PatchesP>
    void intialise_the_feet_on_grids(MultipatchField<CoordFieldOnPatch_host, PatchesP...> feet)
    {
        (intialise_the_feet_on_patch_grid<PatchesP>(feet.template get<PatchesP>()), ...);
    };

    template <class PatchP>
    void intialise_the_feet_on_patch_grid(CoordFieldOnPatch_host<PatchP> feet)
    {
        ddc::for_each(get_idx_range(feet), [&](typename PatchP::Idx12 const& idx) {
            feet(idx) = ddc::coordinate(idx);
        });
    };
};
