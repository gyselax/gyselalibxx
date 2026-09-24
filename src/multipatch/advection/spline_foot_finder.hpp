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
        // class SplineRThetaBuilder1,
        // class SplineRThetaBuilder2,
        // class SplineRThetaBuilder3,
        class MultipatchSplineEvaluator,
        // class MultipatchSplineEvaluatorType,
        class Mappings,
        // class Mapping,
        class... Patches>
// class Patch1,
// class Patch2,
// class Patch3>
class MultipatchSplineFootFinder
{
private:
    using X = typename std::tuple_element<0, Mappings>::type::cartesian_tag_x;
    using Y = typename std::tuple_element<0, Mappings>::type::cartesian_tag_y;

    // using X = typename Mapping::cartesian_tag_x;
    // using Y = typename Mapping::cartesian_tag_y;

    // using R = typename Patch1::Dim1;
    // using Theta = typename Patch1::Dim2;

public:
    template <class Patch>
    using DVectorFieldMem_XY_OnPatch_host
            = host_t<VectorFieldMem<double, typename Patch::IdxRange12, NDTag<X, Y>>>;

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

    // SplineRThetaBuilder1 const& m_builder_1;
    // SplineRThetaBuilder2 const& m_builder_2;
    // SplineRThetaBuilder3 const& m_builder_3;

    MultipatchSplineBuilder& m_builders;
    MultipatchSplineEvaluator const& m_evaluators_advection_field;

    // Mapping const& m_mapping;
    Mappings const& m_mappings;


public:
    MultipatchSplineFootFinder(
            MultipatchTimeStepper const& time_steppers,
            InterfaceDerivativesCalculator& matrices_derivatives_interfaces,
            MultipatchSplineBuilder& builders,
            // SplineRThetaBuilder1 const& builder_1,
            // SplineRThetaBuilder2 const& builder_2,
            // SplineRThetaBuilder3 const& builder_3,
            MultipatchSplineEvaluator const& evaluators_advection_field,
            // Mapping const& mapping)
            Mappings const& mappings)
        : m_time_steppers(time_steppers)
        , m_matrices_derivatives_interfaces(matrices_derivatives_interfaces)
        , m_builders(builders)
        // , m_builder_1(builder_1)
        // , m_builder_2(builder_2)
        // , m_builder_3(builder_3)
        , m_evaluators_advection_field(evaluators_advection_field)
        // , m_mapping(mapping)
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
        // build_spline_from_values(advection_field_splines, const_advection_fields);

        // Apply time stepper to each coordinate field on each patch.
        apply_time_stepper_to_patch<Patch1>(feet, advection_field_splines, dt);
        apply_time_stepper_to_patch<Patch2>(feet, advection_field_splines, dt);
        apply_time_stepper_to_patch<Patch3>(feet, advection_field_splines, dt);
    }



private:
    // void build_spline_from_values(
    //         MultipatchField<VectorSplineCoeffs2D_XY_OnPatch_host, Patch1, Patch2, Patch3> const&
    //                 advection_field_splines,
    //         MultipatchField<DConstVectorField_XY_OnPatch_host, Patch1, Patch2, Patch3> const&
    //                 advection_fields) const
    // {
    //     MultipatchType<IdxRangeOnPatch, Patch1, Patch2, Patch3> const all_idx_ranges
    //             = advection_fields.idx_range();
    //     typename Patch1::IdxRange2 idx_range_theta1(all_idx_ranges.template get<Patch1>());
    //     typename Patch2::IdxRange2 idx_range_theta2(all_idx_ranges.template get<Patch2>());
    //     typename Patch3::IdxRange2 idx_range_theta3(all_idx_ranges.template get<Patch3>());


    //     // STEP 1: DEFINE AND COMPUTE THE DERIVATIVES ================================================
    //     // Instantiate -------------------------------
    //     // --- patch 1
    //     Idx<ddc::Deriv<R>> first_deriv_r1(1);
    //     IdxStep<ddc::Deriv<R>> n_deriv_r1(1);
    //     IdxRange<ddc::Deriv<R>> deriv_r_idx_range1(first_deriv_r1, n_deriv_r1);

    //     IdxRange<ddc::Deriv<R>, typename Patch1::Grid2>
    //             derivs_r_idx_range_1(deriv_r_idx_range1, idx_range_theta1);

    //     // ------ for x component of the advection field
    //     DFieldMem<
    //             IdxRange<ddc::Deriv<typename Patch1::Dim1>, typename Patch1::Grid2>,
    //             Kokkos::HostSpace>
    //             derivs_min_1_alloc_x(derivs_r_idx_range_1);
    //     DField<IdxRange<ddc::Deriv<typename Patch1::Dim1>, typename Patch1::Grid2>,

    //            Kokkos::HostSpace>
    //             derivs_min_1_x = get_field(derivs_min_1_alloc_x);

    //     DFieldMem<
    //             IdxRange<ddc::Deriv<typename Patch1::Dim1>, typename Patch1::Grid2>,
    //             Kokkos::HostSpace>
    //             derivs_max_1_alloc_x(derivs_r_idx_range_1);
    //     DField<IdxRange<ddc::Deriv<typename Patch1::Dim1>, typename Patch1::Grid2>,

    //            Kokkos::HostSpace>
    //             derivs_max_1_x = get_field(derivs_max_1_alloc_x);

    //     // ------ for y component of the advection field
    //     DFieldMem<
    //             IdxRange<ddc::Deriv<typename Patch1::Dim1>, typename Patch1::Grid2>,
    //             Kokkos::HostSpace>
    //             derivs_min_1_alloc_y(derivs_r_idx_range_1);
    //     DField<IdxRange<ddc::Deriv<typename Patch1::Dim1>, typename Patch1::Grid2>,

    //            Kokkos::HostSpace>
    //             derivs_min_1_y = get_field(derivs_min_1_alloc_y);

    //     DFieldMem<
    //             IdxRange<ddc::Deriv<typename Patch1::Dim1>, typename Patch1::Grid2>,
    //             Kokkos::HostSpace>
    //             derivs_max_1_alloc_y(derivs_r_idx_range_1);
    //     DField<IdxRange<ddc::Deriv<typename Patch1::Dim1>, typename Patch1::Grid2>,

    //            Kokkos::HostSpace>
    //             derivs_max_1_y = get_field(derivs_max_1_alloc_y);

    //     // --- patch 2
    //     Idx<ddc::Deriv<R>> first_deriv_r2(1);
    //     IdxStep<ddc::Deriv<R>> n_deriv_r2(1);
    //     IdxRange<ddc::Deriv<R>> deriv_r_idx_range2(first_deriv_r2, n_deriv_r2);

    //     IdxRange<ddc::Deriv<R>, typename Patch2::Grid2>
    //             derivs_r_idx_range_2(deriv_r_idx_range2, idx_range_theta2);

    //     // ------ for x component of the advection field
    //     DFieldMem<
    //             IdxRange<ddc::Deriv<typename Patch2::Dim1>, typename Patch2::Grid2>,
    //             Kokkos::HostSpace>
    //             derivs_min_2_alloc_x(derivs_r_idx_range_2);
    //     DField<IdxRange<ddc::Deriv<typename Patch2::Dim1>, typename Patch2::Grid2>,

    //            Kokkos::HostSpace>
    //             derivs_min_2_x = get_field(derivs_min_2_alloc_x);

    //     DFieldMem<
    //             IdxRange<ddc::Deriv<typename Patch2::Dim1>, typename Patch2::Grid2>,
    //             Kokkos::HostSpace>
    //             derivs_max_2_alloc_x(derivs_r_idx_range_2);
    //     DField<IdxRange<ddc::Deriv<typename Patch2::Dim1>, typename Patch2::Grid2>,

    //            Kokkos::HostSpace>
    //             derivs_max_2_x = get_field(derivs_max_2_alloc_x);

    //     // ------ for y component of the advection field
    //     DFieldMem<
    //             IdxRange<ddc::Deriv<typename Patch2::Dim1>, typename Patch2::Grid2>,
    //             Kokkos::HostSpace>
    //             derivs_min_2_alloc_y(derivs_r_idx_range_2);
    //     DField<IdxRange<ddc::Deriv<typename Patch2::Dim1>, typename Patch2::Grid2>,

    //            Kokkos::HostSpace>
    //             derivs_min_2_y = get_field(derivs_min_2_alloc_y);

    //     DFieldMem<
    //             IdxRange<ddc::Deriv<typename Patch2::Dim1>, typename Patch2::Grid2>,
    //             Kokkos::HostSpace>
    //             derivs_max_2_alloc_y(derivs_r_idx_range_2);
    //     DField<IdxRange<ddc::Deriv<typename Patch2::Dim1>, typename Patch2::Grid2>,

    //            Kokkos::HostSpace>
    //             derivs_max_2_y = get_field(derivs_max_2_alloc_y);

    //     // --- patch 3
    //     Idx<ddc::Deriv<R>> first_deriv_r3(1);
    //     IdxStep<ddc::Deriv<R>> n_deriv_r3(1);
    //     IdxRange<ddc::Deriv<R>> deriv_r_idx_range3(first_deriv_r3, n_deriv_r3);

    //     IdxRange<ddc::Deriv<R>, typename Patch3::Grid2>
    //             derivs_r_idx_range_3(deriv_r_idx_range3, idx_range_theta3);

    //     // ------ for x component of the advection field
    //     DFieldMem<
    //             IdxRange<ddc::Deriv<typename Patch3::Dim1>, typename Patch3::Grid2>,
    //             Kokkos::HostSpace>
    //             derivs_min_3_alloc_x(derivs_r_idx_range_3);
    //     DField<IdxRange<ddc::Deriv<typename Patch3::Dim1>, typename Patch3::Grid2>,

    //            Kokkos::HostSpace>
    //             derivs_min_3_x = get_field(derivs_min_3_alloc_x);

    //     DFieldMem<
    //             IdxRange<ddc::Deriv<typename Patch3::Dim1>, typename Patch3::Grid2>,
    //             Kokkos::HostSpace>
    //             derivs_max_3_alloc_x(derivs_r_idx_range_3);
    //     DField<IdxRange<ddc::Deriv<typename Patch3::Dim1>, typename Patch3::Grid2>,

    //            Kokkos::HostSpace>
    //             derivs_max_3_x = get_field(derivs_max_3_alloc_x);

    //     // ------ for y component of the advection field
    //     DFieldMem<
    //             IdxRange<ddc::Deriv<typename Patch3::Dim1>, typename Patch3::Grid2>,
    //             Kokkos::HostSpace>
    //             derivs_min_3_alloc_y(derivs_r_idx_range_3);
    //     DField<IdxRange<ddc::Deriv<typename Patch3::Dim1>, typename Patch3::Grid2>,

    //            Kokkos::HostSpace>
    //             derivs_min_3_y = get_field(derivs_min_3_alloc_y);

    //     DFieldMem<
    //             IdxRange<ddc::Deriv<typename Patch3::Dim1>, typename Patch3::Grid2>,
    //             Kokkos::HostSpace>
    //             derivs_max_3_alloc_y(derivs_r_idx_range_3);
    //     DField<IdxRange<ddc::Deriv<typename Patch3::Dim1>, typename Patch3::Grid2>,

    //            Kokkos::HostSpace>
    //             derivs_max_3_y = get_field(derivs_max_3_alloc_y);


    //     // Initialise only the boundary values: default null.
    //     ddc::for_each(derivs_r_idx_range_1, [&](auto const& idx) {
    //         derivs_min_1_x(idx) = 0;
    //         derivs_min_1_y(idx) = 0;
    //     });
    //     ddc::for_each(derivs_r_idx_range_3, [&](auto const& idx) {
    //         derivs_max_3_x(idx) = 0;
    //         derivs_max_3_y(idx) = 0;
    //     });


    //     MultipatchField<Derivs1OnPatch_host, Patch1, Patch2, Patch3>
    //             derivs_min_x(derivs_min_1_x, derivs_min_2_x, derivs_min_3_x);
    //     MultipatchField<Derivs1OnPatch_host, Patch1, Patch2, Patch3>
    //             derivs_max_x(derivs_max_1_x, derivs_max_2_x, derivs_max_3_x);

    //     MultipatchField<Derivs1OnPatch_host, Patch1, Patch2, Patch3>
    //             derivs_min_y(derivs_min_1_y, derivs_min_2_y, derivs_min_3_y);
    //     MultipatchField<Derivs1OnPatch_host, Patch1, Patch2, Patch3>
    //             derivs_max_y(derivs_max_1_y, derivs_max_2_y, derivs_max_3_y);

    //     MultipatchField<DConstFieldOnPatch_host, Patch1, Patch2, Patch3> const_advection_fields_x(
    //             get_const_field(ddcHelper::get<X>(advection_fields.template get<Patch1>())),
    //             get_const_field(ddcHelper::get<X>(advection_fields.template get<Patch2>())),
    //             get_const_field(ddcHelper::get<X>(advection_fields.template get<Patch3>())));

    //     MultipatchField<DConstFieldOnPatch_host, Patch1, Patch2, Patch3> const_advection_fields_y(
    //             get_const_field(ddcHelper::get<Y>(advection_fields.template get<Patch1>())),
    //             get_const_field(ddcHelper::get<Y>(advection_fields.template get<Patch2>())),
    //             get_const_field(ddcHelper::get<Y>(advection_fields.template get<Patch3>())));

    //     m_matrices_derivatives_interfaces
    //             .solve(const_advection_fields_x, derivs_min_x, derivs_max_x);
    //     m_matrices_derivatives_interfaces
    //             .solve(const_advection_fields_y, derivs_min_y, derivs_max_y);


    //     // STEP 2: BUILD SPLINE REPRESENTATION =======================================================
    //     // --- X component
    //     VectorSplineCoeffs2D_XY_OnPatch_host<Patch1> advection_field_splines_1
    //             = advection_field_splines.template get<Patch1>();
    //     VectorSplineCoeffs2D_XY_OnPatch_host<Patch2> advection_field_splines_2
    //             = advection_field_splines.template get<Patch2>();
    //     VectorSplineCoeffs2D_XY_OnPatch_host<Patch3> advection_field_splines_3
    //             = advection_field_splines.template get<Patch3>();

    //     m_builder_1(
    //             ddcHelper::get<X>(advection_field_splines_1),
    //             const_advection_fields_x.template get<Patch1>(),
    //             std::optional<DConstField<
    //                     IdxRange<ddc::Deriv<typename Patch1::Dim1>, typename Patch1::Grid2>,

    //                     Kokkos::HostSpace>> {std::nullopt},
    //             std::optional(get_const_field(derivs_max_1_x)));

    //     m_builder_2(
    //             ddcHelper::get<X>(advection_field_splines_2),
    //             const_advection_fields_x.template get<Patch2>(),
    //             std::optional(get_const_field(derivs_min_2_x)),
    //             std::optional(get_const_field(derivs_max_2_x)));

    //     m_builder_3(
    //             ddcHelper::get<X>(advection_field_splines_3),
    //             const_advection_fields_x.template get<Patch3>(),
    //             std::optional(get_const_field(derivs_min_3_x)),
    //             std::optional<DConstField<
    //                     IdxRange<ddc::Deriv<typename Patch3::Dim1>, typename Patch3::Grid2>,

    //                     Kokkos::HostSpace>> {std::nullopt});

    //     // --- Y component
    //     m_builder_1(
    //             ddcHelper::get<Y>(advection_field_splines_1),
    //             const_advection_fields_y.template get<Patch1>(),
    //             std::optional<DConstField<
    //                     IdxRange<ddc::Deriv<typename Patch1::Dim1>, typename Patch1::Grid2>,

    //                     Kokkos::HostSpace>> {std::nullopt},
    //             std::optional(get_const_field(derivs_max_1_y)));

    //     m_builder_2(
    //             ddcHelper::get<Y>(advection_field_splines_2),
    //             const_advection_fields_y.template get<Patch2>(),
    //             std::optional(get_const_field(derivs_min_2_y)),
    //             std::optional(get_const_field(derivs_max_2_y)));

    //     m_builder_3(
    //             ddcHelper::get<Y>(advection_field_splines_3),
    //             const_advection_fields_y.template get<Patch3>(),
    //             std::optional(get_const_field(derivs_min_3_y)),
    //             std::optional<DConstField<
    //                     IdxRange<ddc::Deriv<typename Patch3::Dim1>, typename Patch3::Grid2>,

    //                     Kokkos::HostSpace>> {std::nullopt});
    // }



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
        CoordFieldMemOnPatch_host<Patch> feet_copy_on_patch(
                get_idx_range(feet.template get<Patch>()));
        ddc::parallel_deepcopy(get_field(feet_copy_on_patch), feet.template get<Patch>());
        m_time_steppers.template get<Patch>()
                .update(Kokkos::DefaultHostExecutionSpace(),
                        get_field(feet_copy_on_patch),
                        dt,
                        dy,
                        update_function);
        ddc::parallel_deepcopy(feet.template get<Patch>(), get_field(feet_copy_on_patch));

        is_unified<Patch>(feet.template get<Patch>());
    }



    template <class Patch>
    void advect_feet_on_patch(
            CoordFieldOnPatch_host<Patch> const& feet_log,
            DConstVectorField_XY_OnPatch_host<Patch> const& advection_field,
            double const dt) const
    {
        using CoordLog = typename Patch::Coord12;
        using Dim2 = typename Patch::Dim2;

        IdxRangeOnPatch<Patch> const idx_range
                = get_idx_range<typename Patch::Grid1, typename Patch::Grid2>(feet_log);

        ddc::for_each(idx_range_rp, [&](typename Patch::Idx12 const idx) {
            CoordLog const coord_log(feet_log(idx));
            Coord<X, Y> const coord_xy = m_mapping(coord_log);

            Coord<X, Y> const feet_xy = coord_xy - dt * advection_field(idx);
        });

        // Correction for the periodicity.
        ddc::for_each(idx_range, [&](typename Patch::Idx12 const idx) {
            feet_log(idx) = m_mapping(feet_xy);
            ddc::select<Dim2>(feet_log(idx)) = ddcHelper::restrict_to_idx_range(
                    ddc::select<Dim2>(feet_log(idx)),
                    typename Patch::IdxRange2(idx_range));
        });
    }


    template <class Patch>
    void advect_feet_on_patch_with_Opoint(
            CoordFieldOnPatch_host<Patch> const& feet_rp,
            DConstVectorField_XY_OnPatch_host<Patch> const& advection_field,
            double const dt) const
    {
        using Coord12 = typename Patch::Coord12;
        using Dim2 = typename Patch::Dim2;

        IdxRangeOnPatch<Patch> const idx_range_rp
                = get_idx_range<typename Patch::Grid1, typename Patch::Grid2>(feet_rp);
        // auto mapping = m_mappings.template get<Patch>();
        Coord<X, Y> coord_Opoint(m_mapping(Coord12(0, 0)));

        ddc::for_each(idx_range_rp, [&](typename Patch::Idx12 const irp) {
            Coord12 const coord_rp(feet_rp(irp));
            Coord<X, Y> const coord_xy = m_mapping(coord_rp);

            Coord<X, Y> const feet_xy = coord_xy - dt * advection_field(irp);

            // Correction on theta.
            if (norm_inf(feet_xy - coord_Opoint) < 1e-15) {
                feet_rp(irp) = Coord12(0, 0);
            } else {
                feet_rp(irp) = m_mapping(feet_xy);
                ddc::select<Dim2>(feet_rp(irp)) = ddcHelper::restrict_to_idx_range(
                        ddc::select<Dim2>(feet_rp(irp)),
                        typename Patch::IdxRange2(idx_range_rp));
            }
        });



        template <class Patch, class T>
        void is_unified(
                Field<T,
                      typename Patch::IdxRange12,

                      Kokkos::HostSpace> const& values) const
        {
            typename Patch::IdxRange1 const r_idx_range
                    = get_idx_range<typename Patch::Grid1>(values);
            typename Patch::IdxRange2 const theta_idx_range
                    = get_idx_range<typename Patch::Grid2>(values);
            if (std::fabs(ddc::coordinate(r_idx_range.front())) < 1e-15) {
                ddc::for_each(theta_idx_range, [&](const typename Patch::Idx2 ip) {
                    if (norm_inf(
                                values(r_idx_range.front(), ip)
                                - values(r_idx_range.front(), theta_idx_range.front()))
                        > 1e-15) {
                        std::cout << "WARNING ! -> Discontinous at the center point." << std::endl;
                    }
                    assert(values(r_idx_range.front(), ip)
                           == values(r_idx_range.front(), theta_idx_range.front()));
                });
            }
        }



        template <class Patch, class T>
        void unify_value_at_center_pt(Field<T, typename Patch::IdxRange12> values) const
        {
            typename Patch::IdxRange1 const r_idx_range
                    = get_idx_range<typename Patch::Grid1>(values);
            typename Patch::IdxRange2 const theta_idx_range
                    = get_idx_range<typename Patch::Grid2>(values);
            if (std::fabs(ddc::coordinate(r_idx_range.front())) < 1e-15) {
                ddc::for_each(theta_idx_range, [&](const typename Patch::Idx2 ip) {
                    values(r_idx_range.front(), ip)
                            = values(r_idx_range.front(), theta_idx_range.front());
                });
            }
        }

        template <class... PatchesP>
        void
        copy(MultipatchField<ConstSplineCoeffOnPatch_2D_host, PatchesP...> field_x,
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
