// SPDX-License-Identifier: MIT
#pragma once
#include <functional>

#include <ddc/kernels/splines.hpp>

#include "ddc_alias_inline_functions.hpp"
#include "ddc_aliases.hpp"
#include "multipatch_field.hpp"
#include "multipatch_field_mem.hpp"
#include "multipatch_type.hpp"
#include "types.hpp"
#include "utils_tools.hpp"


template <
        class MultipatchTimeStepper,
        class InterfaceDerivativesCalculator,
        class MultipatchSplineBuilder,
        class MultipatchSplineEvaluator,
        class Mappings,
        class... Patches>
class MultipatchSplineFootFinder
{
private:
    using X = typename std::tuple_element<0, Mappings>::type::cartesian_tag_x;
    using Y = typename std::tuple_element<0, Mappings>::type::cartesian_tag_y;

    using PatchSeq = ddc::detail::TypeSeq<Patches...>;
    // TODO: adapt to real geometries.
    using PatchOpt = ddc::type_seq_element_t<0, PatchSeq>;

public:
    template <class Patch>
    using DVectorField_XY_OnPatch_host
            = host_t<VectorField<double, typename Patch::IdxRange12, NDTag<X, Y>>>;

    template <class Patch>
    using DConstVectorField_XY_OnPatch_host
            = host_t<VectorConstField<double, typename Patch::IdxRange12, NDTag<X, Y>>>;


    template <class Patch>
    using VectorSplineCoeffsMem2D_XY_OnPatch_host
            = host_t<VectorFieldMem<double, typename Patch::IdxRangeBS12, NDTag<X, Y>>>;

    template <class Patch>
    using VectorSplineCoeffs2D_XY_OnPatch_host
            = host_t<VectorField<double, typename Patch::IdxRangeBS12, NDTag<X, Y>>>;


private:
    MultipatchTimeStepper const& m_time_steppers;

    InterfaceDerivativesCalculator& m_matrices_derivatives_interfaces;

    MultipatchSplineBuilder& m_builders;
    MultipatchSplineEvaluator const& m_evaluators_advection_field;

    Mappings const& m_mappings;


public:
    MultipatchSplineFootFinder(
            MultipatchTimeStepper const& time_steppers,
            InterfaceDerivativesCalculator& matrices_derivatives_interfaces,
            MultipatchSplineBuilder& builders,
            MultipatchSplineEvaluator const& evaluators_advection_field,
            Mappings const& mappings)
        : m_time_steppers(time_steppers)
        , m_matrices_derivatives_interfaces(matrices_derivatives_interfaces)
        , m_builders(builders)
        , m_evaluators_advection_field(evaluators_advection_field)
        , m_mappings(mappings)
    {
    }

    ~MultipatchSplineFootFinder() {};


    void operator()(
            MultipatchField<CoordFieldOnPatch_host, Patches...> feet,
            MultipatchField<DConstVectorField_XY_OnPatch_host, Patches...> const_advection_fields,
            double const dt) const
    {
        // Define splines for the advection field
        // TODO: implement MultipatchSplineBuilder::spline_idx_ranges().
        MultipatchField<VectorSplineCoeffsMem2D_XY_OnPatch_host, Patches...>
                advection_field_splines_alloc(m_builders.spline_idx_ranges());

        MultipatchField<VectorSplineCoeffs2D_XY_OnPatch_host, Patches...> advection_field_splines(
                advection_field_splines_alloc);


        // Compute derivatives of the advection field
        // and build spline representations of the advection field.
        // TODO: implement MultipatchSplineBuilder::operator()().
        m_builders(advection_field_splines, const_advection_fields);

        // Apply time stepper to each coordinate field on each patch.
        (apply_time_stepper_to_patch<Patches>(feet, advection_field_splines, dt), ...);
    }



private:
    template <class Patch>
    void apply_time_stepper_to_patch(
            MultipatchField<CoordFieldOnPatch_host, Patches...> feet,
            MultipatchField<VectorSplineCoeffs2D_XY_OnPatch_host, Patches...>
                    advection_field_splines,
            double const dt) const
    {
        // The function describing how the derivative of the evolve function is calculated.
        std::function<void(DVectorField_XY_OnPatch_host<Patch>, CoordConstFieldOnPatch_host<Patch>)>
                dy = [&](DVectorField_XY_OnPatch_host<Patch> updated_advection_field,
                         CoordConstFieldOnPatch_host<Patch> feet_on_patch) {
                    ddc::parallel_deepcopy(feet.template get<Patch>(), feet_on_patch);
                    MultipatchField<CoordConstFieldOnPatch_host, Patches...> const_feet
                            = get_const_field(feet);

                    // Split into two: for X and for Y projections.
                    // TODO: maybe implement a MultipatchVectorField to be able to easily split?
                    MultipatchFieldMem<DFieldMemOnPatch_host, Patches...> advection_fields_x_alloc(
                            feet.idx_range());
                    MultipatchFieldMem<DFieldMemOnPatch_host, Patches...> advection_fields_y_alloc(
                            feet.idx_range());
                    MultipatchField<DFieldOnPatch_host, Patches...> advection_fields_x(
                            advection_fields_x_alloc);
                    MultipatchField<DFieldOnPatch_host, Patches...> advection_fields_y(
                            advection_fields_y_alloc);

                    MultipatchFieldMem<SplineCoeffMemOnPatch_2D_host, Patches...>
                            advection_field_splines_x_alloc(m_builders.spline_idx_ranges());
                    MultipatchFieldMem<SplineCoeffMemOnPatch_2D_host, Patches...>
                            advection_field_splines_y_alloc(m_builders.spline_idx_ranges());
                    MultipatchField<ConstSplineCoeffOnPatch_2D_host, Patches...>
                            advection_field_splines_x(advection_field_splines_x_alloc);
                    MultipatchField<ConstSplineCoeffOnPatch_2D_host, Patches...>
                            advection_field_splines_y(advection_field_splines_y_alloc);

                    copy(advection_field_splines_x,
                         advection_field_splines_y,
                         advection_field_splines);

                    m_evaluators_advection_field(
                            advection_fields_x,
                            const_feet,
                            advection_field_splines_x);
                    m_evaluators_advection_field(
                            advection_fields_y,
                            const_feet,
                            advection_field_splines_y);

                    ddc::for_each(get_idx_range(updated_advection_field), [&](auto const& idx) {
                        ddcHelper::get<X>(updated_advection_field)(idx)
                                = advection_fields_x.template get<Patch>()(idx);
                        ddcHelper::get<Y>(updated_advection_field)(idx)
                                = advection_fields_y.template get<Patch>()(idx);
                    });
                };

        // The function describing how the value(s) are updated using the derivative.
        std::function<
                void(CoordFieldOnPatch_host<Patch>,
                     DConstVectorField_XY_OnPatch_host<Patch>,
                     double)>
                update_function = [&](CoordFieldOnPatch_host<Patch> feet_on_patch,
                                      DConstVectorField_XY_OnPatch_host<Patch> advection_field,
                                      double dt) {
                    advect_feet_on_patch<Patch>(feet_on_patch, advection_field, dt);
                };

        // Solve the characteristic equation
        CoordFieldMemOnPatch_host<Patch> feet_copy_on_patch_alloc(
                get_idx_range(feet.template get<Patch>()));
        CoordFieldOnPatch_host<Patch> feet_copy_on_patch(feet_copy_on_patch_alloc);
        ddc::parallel_deepcopy(feet_copy_on_patch, feet.template get<Patch>());
        m_time_steppers.template get<Patch>()
                .update(Kokkos::DefaultHostExecutionSpace(),
                        feet_copy_on_patch,
                        dt,
                        dy,
                        update_function);
        ddc::parallel_deepcopy(feet.template get<Patch>(), feet_copy_on_patch);

        is_unified<Patch>(feet.template get<Patch>());
    }



    /**
     * Advect in the physical domain, but the feet are computed back in the logical 
     * domain. 
     * TODO: separate the two operations? 
     */
    template <class Patch>
    void advect_feet_on_patch(
            CoordFieldOnPatch_host<Patch> const& feet_log,
            DConstVectorField_XY_OnPatch_host<Patch> const& advection_field,
            double const dt) const
    {
        IdxRangeOnPatch<Patch> const idx_range
                = get_idx_range<typename Patch::Grid1, typename Patch::Grid2>(feet_log);

        // TODO: inverse mapping.
        constexpr int PatchIdx = ddc::type_seq_rank_v<Patch, PatchSeq>;
        auto mapping = std::get<PatchIdx>(m_mappings);
        auto inv_mapping;

        constexpr bool has_an_Opoint = std::is_same_v<Patch, PatchOpt>;

        ddc::for_each(idx_range_rp, [&](typename Patch::Idx12 const idx) {
            Coord<X, Y> const coord_xy = mapping(feet_log(idx));
            Coord<X, Y> const feet_xy = coord_xy - dt * advection_field(idx);

            if constexpr (has_an_Opoint) {
                typename Patch::Coord12 coord_Opoint_log(0, 0);
                if (norm_inf(feet_xy - mapping(coord_Opoint_log)) < 1e-15) {
                    feet_log(idx) = coord_Opoint_log;
                }
            }

            feet_log(idx) = inv_mapping(feet_xy);
        });


        // Correction for the periodicity.
        if constexpr (typename Patch::Dim1::PERIODIC) {
            using PDim1 = typename Patch::Dim1;
            ddc::for_each(idx_range, [&](typename Patch::Idx12 const idx) {
                ddc::select<PDim1>(feet_log(idx)) = ddcHelper::restrict_to_idx_range(
                        ddc::select<PDim1>(feet_log(idx)),
                        typename Patch::IdxRange1(idx_range));
            });
        }
        if constexpr (typename Patch::Dim2::PERIODIC) {
            using PDim2 = typename Patch::Dim2;
            ddc::for_each(idx_range, [&](typename Patch::Idx12 const idx) {
                ddc::select<PDim2>(feet_log(idx)) = ddcHelper::restrict_to_idx_range(
                        ddc::select<PDim2>(feet_log(idx)),
                        typename Patch::IdxRange2(idx_range));
            });
        }


        if constexpr (has_an_Opoint) {
            // Treatment to conserve the C0 property of the advected function:
            unify_value_at_center_pt<Patch>(feet_log);
            // Test if the values are the same at the center point
            is_unified<Patch>(feet_log);
        }
    }


    // TODO: maybe remove this function.
    template <class Patch, class T>
    void is_unified(Field<T, typename Patch::IdxRange12, Kokkos::HostSpace> const& values) const
    {
        typename Patch::IdxRange1 const idx_range_1(get_idx_range(values));
        typename Patch::IdxRange2 const idx_range_2(get_idx_range(values));
        typename Patch::Idx1 const idx1_min = idx_range_1.front();
        typename Patch::Idx2 const idx2_min = idx_range_2.front();
        // We check that there is an O-point?
        assert(std::fabs(ddc::coordinate(idx1_min)) < 1e-15);
        ddc::for_each(idx_range_2, [&](const typename Patch::Idx2 idx2) {
            if (norm_inf(values(idx1_min, idx2) - values(idx1_min, idx2_min)) > 1e-15) {
                std::cout << "WARNING ! -> Discontinous at the center point." << std::endl;
            }
            assert(values(idx1_min, idx2) == values(idx1_min, idx2_min));
        });
    }



    template <class Patch, class T>
    void unify_value_at_center_pt(Field<T, typename Patch::IdxRange12> values) const
    {
        typename Patch::IdxRange1 const idx_range_1(get_idx_range(values));
        typename Patch::IdxRange2 const idx_range_2(get_idx_range(values));
        typename Patch::Idx1 const idx1_min = idx_range_1.front();
        typename Patch::Idx2 const idx2_min = idx_range_2.front();
        // We check that there is an O-point?
        assert(std::fabs(ddc::coordinate(idx1_min)) < 1e-15);
        ddc::for_each(idx_range_2, [&](const typename Patch::Idx2 idx2) {
            values(idx1_min, idx2) = values(idx1_min, idx2_min);
        });
    }

    template <class... PatchesP>
    void copy(
            MultipatchField<ConstSplineCoeffOnPatch_2D_host, PatchesP...> field_x,
            MultipatchField<ConstSplineCoeffOnPatch_2D_host, PatchesP...> field_y,
            MultipatchField<VectorSplineCoeffs2D_XY_OnPatch_host, PatchesP...> field_xy)
    {
        (copy_on_patch<PatchesP>(
                 field_x.template get<PatchesP>(),
                 field_y.template get<PatchesP>(),
                 field_xy.template get<PatchesP>()),
         ...);
    }

    template <class PatchP>
    void copy_on_patch(
            ConstSplineCoeffOnPatch_2D_host<PatchP> field_x,
            ConstSplineCoeffOnPatch_2D_host<PatchP> field_y,
            VectorSplineCoeffs2D_XY_OnPatch_host<PatchP> field_xy)
    {
        deep_copy(field_x, ddcHelper::get<X>(field_xy));
        deep_copy(field_y, ddcHelper::get<Y>(field_xy));
    }
};
