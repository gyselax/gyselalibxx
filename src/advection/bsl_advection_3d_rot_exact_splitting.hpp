// SPDX-License-Identifier: MIT
#pragma once
#include <ddc/ddc.hpp>

#include "ddc_alias_inline_functions.hpp"
#include "ddc_aliases.hpp"
#include "ddc_helper.hpp"
#include "iadvectionrot3d.hpp"
#include "iinterpolator.hpp"
#include "species_info.hpp"
#include <cmath>

#include "exact_splitting_solver.h"
#include <Eigen/Dense>
#include <iostream>

using namespace std;


/**
 * @brief A class which computes the velocity advections of a 3D shifted rotation. Working for every cartesian geometry.
 */

 template <class Geometry, class GridX, class GridVx, class GridVy, class GridVz>
 class BslAdvectionVelocityRot3DExact : public IAdvectionVelocityRot3D<Geometry, GridVx, GridVy, GridVz>
 {
     using IdxRangeFdistribu = typename Geometry::IdxRangeFdistribu;
     using IdxRangeSpatial = typename Geometry::IdxRangeSpatial;
     using IdxRangeSpSpatial = typename Geometry::IdxRangeSpSpatial;
     using IdxSpatial = typename IdxRangeSpatial::discrete_element_type;
     using IdxVx = Idx<GridVx>;
     using IdxVy = Idx<GridVy>;
     using IdxVz = Idx<GridVz>;
     using IdxSpX = Idx<Species, GridX>;
     using DimVx = typename GridVx::continuous_dimension_type;
     using DimVy = typename GridVy::continuous_dimension_type;
     using DimVz = typename GridVz::continuous_dimension_type;
     using IdxRangeSpaceVelocity
             = ddc::remove_dims_of_t<typename Geometry::IdxRangeFdistribu, Species>;
     using IdxRangeVelocity
             = ddc::remove_dims_of_t<typename Geometry::IdxRangeFdistribu, Species, GridX>;
 
     //using GridVy = typename Geometry::template the_other_velocity_dim<GridV>;
 
 private:
    // vx
     using PreallocatableInterpolatorType_Vx = interpolator_on_idx_range_t<
             IPreallocatableInterpolator,
             GridVx,
             IdxRangeVelocity>;
     using InterpolatorType_Vx
             = interpolator_on_idx_range_t<IInterpolator, GridVx, IdxRangeVelocity>;
     PreallocatableInterpolatorType_Vx const& m_interpolator_vx;

     // vy
    using PreallocatableInterpolatorType_Vy = interpolator_on_idx_range_t<
             IPreallocatableInterpolator,
             GridVy,
             IdxRangeVelocity>;
     using InterpolatorType_Vy
             = interpolator_on_idx_range_t<IInterpolator, GridVy, IdxRangeVelocity>;
     PreallocatableInterpolatorType_Vy const& m_interpolator_vy;

     // vz
     using PreallocatableInterpolatorType_Vz = interpolator_on_idx_range_t<
             IPreallocatableInterpolator,
             GridVz,
             IdxRangeVelocity>;
     using InterpolatorType_Vz
             = interpolator_on_idx_range_t<IInterpolator, GridVz, IdxRangeVelocity>;
     PreallocatableInterpolatorType_Vz const& m_interpolator_vz;


 public:
     /**
      * @brief Constructor 
      * @param[in] interpolator_v interpolator along the GridV direction which refers to the velocity space.  
      */
     explicit BslAdvectionVelocityRot3DExact(PreallocatableInterpolatorType_Vx const& interpolator_vx, 
                                             PreallocatableInterpolatorType_Vy const& interpolator_vy,
                                             PreallocatableInterpolatorType_Vz const& interpolator_vz)
         : m_interpolator_vx(interpolator_vx)
         , m_interpolator_vy(interpolator_vy)
         , m_interpolator_vz(interpolator_vz)
     {
     }
 
     ~BslAdvectionVelocityRot3DExact() override = default;
 
     /**
      * @brief Advects fdistribu along GridV for a duration dt.
      * @param[in, out] allfdistribu Reference to the whole distribution function for one species, allocated on the device (ie it lets the choice of the location depend on the build configuration).
      * @param[in] magnetic_field_z Reference to the magnetic field.
      * @param[in] mean_velocity_x Reference to the velocity shift in vx.
      * @param[in] mean_velocity_y Reference to the velocity shift in vy.
      * @param[in] dt Time step
      * @return A reference to the allfdistribu array containing the value of the function at the coordinates.
      */
     Field<double, IdxRangeFdistribu> operator()(
             Field<double, IdxRangeFdistribu> const allfdistribu,
             Field<const double, IdxRangeSpatial> const magnetic_field_x,
             Field<const double, IdxRangeSpatial> const magnetic_field_y,
             Field<const double, IdxRangeSpatial> const magnetic_field_z,
             Field<const double, IdxRangeSpatial> const mean_velocity_x,
             Field<const double, IdxRangeSpatial> const mean_velocity_y,
             Field<const double, IdxRangeSpatial> const mean_velocity_z,
             double const dt) const override
     {
         using IdxRangeBatch_Vx = ddc::remove_dims_of_t<IdxRangeFdistribu, Species, GridX, GridVx>;
         using IdxRangeBatch_Vy = ddc::remove_dims_of_t<IdxRangeFdistribu, Species, GridX, GridVy>;
         using IdxRangeBatch_Vz = ddc::remove_dims_of_t<IdxRangeFdistribu, Species, GridX, GridVz>;
         
         using IdxBatch_Vx = typename IdxRangeBatch_Vx::discrete_element_type;
         using IdxBatch_Vy = typename IdxRangeBatch_Vy::discrete_element_type;
         using IdxBatch_Vz = typename IdxRangeBatch_Vz::discrete_element_type;
 
         Kokkos::Profiling::pushRegion("BslAdvectionVelocityRot3DExact");
         IdxRangeFdistribu const idx_range = get_idx_range(allfdistribu);
         IdxRange<GridVx> const idx_range_vx = ddc::select<GridVx>(idx_range);
         IdxRange<GridVy> const idx_range_vy = ddc::select<GridVy>(idx_range);
         IdxRange<GridVz> const idx_range_vz = ddc::select<GridVz>(idx_range);
         IdxRange<Species> const idx_range_sp = ddc::select<Species>(idx_range);
         IdxRange<GridX> const idx_range_x = ddc::select<GridX>(idx_range);

         IdxRange<Species,GridX> const idx_range_spx = ddc::select<Species,GridX>(idx_range);
        // vx
         FieldMem<double, typename InterpolatorType_Vx::batched_derivs_idx_range_type> derivs_min_vx(
                 m_interpolator_vx.batched_derivs_idx_range_xmin(
                         ddc::remove_dims_of<Species, GridX>(idx_range)));
         FieldMem<double, typename InterpolatorType_Vx::batched_derivs_idx_range_type> derivs_max_vx(
                 m_interpolator_vx.batched_derivs_idx_range_xmax(
                         ddc::remove_dims_of<Species, GridX>(idx_range)));
         ddc::parallel_fill(derivs_min_vx, 0.);
         ddc::parallel_fill(derivs_max_vx, 0.);
        // vy 
         FieldMem<double, typename InterpolatorType_Vy::batched_derivs_idx_range_type> derivs_min_vy(
                 m_interpolator_vy.batched_derivs_idx_range_xmin(
                         ddc::remove_dims_of<Species, GridX>(idx_range)));
         FieldMem<double, typename InterpolatorType_Vy::batched_derivs_idx_range_type> derivs_max_vy(
                 m_interpolator_vy.batched_derivs_idx_range_xmax(
                         ddc::remove_dims_of<Species, GridX>(idx_range)));
         ddc::parallel_fill(derivs_min_vy, 0.);
         ddc::parallel_fill(derivs_max_vy, 0.);
        // vz
         FieldMem<double, typename InterpolatorType_Vz::batched_derivs_idx_range_type> derivs_min_vz(
                 m_interpolator_vz.batched_derivs_idx_range_xmin(
                         ddc::remove_dims_of<Species, GridX>(idx_range)));
         FieldMem<double, typename InterpolatorType_Vz::batched_derivs_idx_range_type> derivs_max_vz(
                 m_interpolator_vz.batched_derivs_idx_range_xmax(
                         ddc::remove_dims_of<Species, GridX>(idx_range)));
         ddc::parallel_fill(derivs_min_vz, 0.);
         ddc::parallel_fill(derivs_max_vz, 0.);
  
         // pre-allocate some memory to prevent allocation later in loop
         
         IdxRangeVelocity batched_feet_idx_range(idx_range);
         // vx 
         FieldMem<Coord<DimVx>, IdxRangeVelocity> feet_coords_alloc_vx(batched_feet_idx_range);
         Field<Coord<DimVx>, IdxRangeVelocity> feet_coords_vx(get_field(feet_coords_alloc_vx));
         std::unique_ptr<InterpolatorType_Vx> const interpolator_v_ptr_vx = m_interpolator_vx.preallocate();
         InterpolatorType_Vx const& interpolator_vx = *interpolator_v_ptr_vx;
         // vy 
         FieldMem<Coord<DimVy>, IdxRangeVelocity> feet_coords_alloc_vy(batched_feet_idx_range);
         Field<Coord<DimVy>, IdxRangeVelocity> feet_coords_vy(get_field(feet_coords_alloc_vy));
         std::unique_ptr<InterpolatorType_Vy> const interpolator_v_ptr_vy = m_interpolator_vy.preallocate();
         InterpolatorType_Vy const& interpolator_vy = *interpolator_v_ptr_vy;
         // vx 
         FieldMem<Coord<DimVz>, IdxRangeVelocity> feet_coords_alloc_vz(batched_feet_idx_range);
         Field<Coord<DimVz>, IdxRangeVelocity> feet_coords_vz(get_field(feet_coords_alloc_vz));
         std::unique_ptr<InterpolatorType_Vz> const interpolator_v_ptr_vz = m_interpolator_vz.preallocate();
         InterpolatorType_Vz const& interpolator_vz = *interpolator_v_ptr_vz;
 
         IdxRangeSpatial const idx_range_spatial(get_idx_range(allfdistribu));
 
         IdxRangeBatch_Vx batch_idx_range_vx(idx_range);
         IdxRangeBatch_Vy batch_idx_range_vy(idx_range);
         IdxRangeBatch_Vz batch_idx_range_vz(idx_range);


         // define magnetic field in host
         ddc::Chunk host_magnetic_x = ddc::create_mirror(magnetic_field_x.span_cview());
         ddc::Chunk host_magnetic_y = ddc::create_mirror(magnetic_field_y.span_cview());
         ddc::Chunk host_magnetic_z = ddc::create_mirror(magnetic_field_z.span_cview());

         ddc::parallel_deepcopy(
            get_field(host_magnetic_x),
            magnetic_field_x);
         ddc::parallel_deepcopy(
            get_field(host_magnetic_y),
            magnetic_field_y);
         ddc::parallel_deepcopy(
            get_field(host_magnetic_z),
            magnetic_field_z);

 
        ddc::for_each(idx_range_spx, [&](IdxSpX const ispx) {
            IdxSp const isp(ispx);
            IdxSpatial const ix(ispx);
            double q_over_m = -charge(isp)/mass(isp);
        //ddc::for_each(idx_range_x, [&](IdxSpatial const ix) {

            double Bfield[3] = {q_over_m * host_magnetic_x(ix), q_over_m * host_magnetic_y(ix), q_over_m * host_magnetic_z(ix)};
            //double B_abs = std::sqrt(Bfield[0] * Bfield[0] + Bfield[1] * Bfield[1] + Bfield[2] * Bfield[2]);
            double tol = 1e-13;
            int count = 0;
            for (int i = 0; i < 3; ++i)
            {
                if (std::abs(Bfield[i]) > tol)
                ++count;
            }
            //std::cout << "non-zero number is: " << count << std::endl;
            // 原始索引
            std::vector<int> idx_order = {0, 1, 2};

            // 按绝对值从小到大排序索引
            std::sort(idx_order.begin(), idx_order.end(),
                  [&Bfield](int i1, int i2) {
                      return std::abs(Bfield[i1]) < std::abs(Bfield[i2]);
                  });

            //std::cout << "idx_order x is: " << idx_order[0] << std::endl;
            //std::cout << "idx_order y is: " << idx_order[1] << std::endl;
            //std::cout << "idx_order z is: " << idx_order[2] << std::endl;

            
            // module vx 
            auto module_vx = [&](double coef1, double coef2){
            ddc::parallel_for_each(Kokkos::DefaultExecutionSpace(), batch_idx_range_vx,
                KOKKOS_LAMBDA(IdxBatch_Vx const ib){
                    IdxVy const ivy(ib);
                    IdxVz const ivz(ib);
                    double dvx = (ddc::coordinate(ivy)-mean_velocity_y(ix))*coef1
                           + (ddc::coordinate(ivz)-mean_velocity_z(ix))*coef2;
                    for (IdxVx const iv : idx_range_vx)
                        feet_coords_vx(iv, ib) = Coord<DimVx>(ddc::coordinate(iv)+dvx);
                });
            interpolator_vx(allfdistribu[isp][ix],
                        get_const_field(feet_coords_vx),
                        get_const_field(derivs_min_vx),
                        get_const_field(derivs_max_vx));
            };

            // module vy 
            auto module_vy = [&](double coef1, double coef2){
            ddc::parallel_for_each(Kokkos::DefaultExecutionSpace(), batch_idx_range_vy,
                KOKKOS_LAMBDA(IdxBatch_Vy const ib){
                    IdxVx const ivx(ib);
                    IdxVz const ivz(ib);
                    double dvy = (ddc::coordinate(ivx)-mean_velocity_x(ix))*coef1
                           + (ddc::coordinate(ivz)-mean_velocity_z(ix))*coef2;
                    for (IdxVy const iv : idx_range_vy)
                        feet_coords_vy(iv, ib) = Coord<DimVy>(ddc::coordinate(iv)+dvy);
                });
            interpolator_vy(allfdistribu[isp][ix],
                        get_const_field(feet_coords_vy),
                        get_const_field(derivs_min_vy),
                        get_const_field(derivs_max_vy));
            };

            // module vz 
            auto module_vz = [&](double coef1, double coef2){
            ddc::parallel_for_each(Kokkos::DefaultExecutionSpace(), batch_idx_range_vz,
                KOKKOS_LAMBDA(IdxBatch_Vz const ib){
                    IdxVx const ivx(ib);
                    IdxVy const ivy(ib);
                    double dvz = (ddc::coordinate(ivx)-mean_velocity_x(ix))*coef1
                           + (ddc::coordinate(ivy)-mean_velocity_y(ix))*coef2;
                    for (IdxVz const iv : idx_range_vz)
                        feet_coords_vz(iv, ib) = Coord<DimVz>(ddc::coordinate(iv)+dvz);
                });
            interpolator_vz(allfdistribu[isp][ix],
                        get_const_field(feet_coords_vz),
                        get_const_field(derivs_min_vz),
                        get_const_field(derivs_max_vz));
            };
            

            // ======= compute the rotation matrix =======
            /*
            double coe1 = std::sin(dt*B_abs) / B_abs; 
            double coe2 = 2.0 * std::sin(0.5 * dt * B_abs) * std::sin(0.5 * dt * B_abs) / (B_abs * B_abs);
            double M11 = 1.0 + coe2 * (-Bfield[1]*Bfield[1] - Bfield[2]*Bfield[2]);
            double M12 = coe1 *    Bfield[2] + coe2 * Bfield[0] * Bfield[1];
            double M13 = coe1 * (-Bfield[1]) + coe2 * Bfield[0] * Bfield[2];
            double M21 = coe1 * (-Bfield[2]) + coe2 * Bfield[0] * Bfield[1];
            double M22 = 1.0 + coe2 * (-Bfield[0]*Bfield[0] - Bfield[2]*Bfield[2]);
            double M23 = coe1 *    Bfield[0] + coe2 * Bfield[1] * Bfield[2];
            double M31 = coe1 *    Bfield[1] + coe2 * Bfield[0] * Bfield[2];
            double M32 = coe1 * (-Bfield[0]) + coe2 * Bfield[1] * Bfield[2];
            double M33 = 1.0 + coe2 * (-Bfield[0]*Bfield[0] - Bfield[1]*Bfield[1]); 
            */

            
            
            std::array<double, 3> magnetic_field = {Bfield[0], Bfield[1], Bfield[2]};
            std::array<double, 3> yl, y2, y3, yr;
            
            // start the advections based on the rotation 
            if (count >= 2 && idx_order[0] == 0 && idx_order[1] == 1 && idx_order[2] == 2)
            {
                std::array<int, 3> order = {0, 1, 2};

                exact_splitting_solver(magnetic_field, dt, yl, y2, y3, yr, order);

                double a1 = dt*yl[1];
                double a2 = dt*yl[2];
                double a3 = dt*y2[0];
                double a4 = dt*y2[2];
                double a5 = dt*y3[0];
                double a6 = dt*y3[1];
                double a7 = dt*yr[1];
                double a8 = dt*yr[2];

                module_vx(a7,a8);
                module_vz(a5,a6);
                module_vy(a3,a4);
                module_vx(a1,a2);
                /*
                // module vx 
                ddc::parallel_for_each(Kokkos::DefaultExecutionSpace(), batch_idx_range_vx,
                    KOKKOS_LAMBDA(IdxBatch_Vx const ib){
                        IdxVy const ivy(ib);
                        IdxVz const ivz(ib);
                        double dvx = (ddc::coordinate(ivy)-mean_velocity_y(ix))*a7
                            + (ddc::coordinate(ivz)-mean_velocity_z(ix))*a8;
                        for (IdxVx const iv : idx_range_vx)
                            feet_coords_vx(iv, ib) = Coord<DimVx>(ddc::coordinate(iv)+dvx);
                    });
                interpolator_vx(allfdistribu[isp][ix],
                            get_const_field(feet_coords_vx),
                            get_const_field(derivs_min_vx),
                            get_const_field(derivs_max_vx));
                
                // module vz
                ddc::parallel_for_each(Kokkos::DefaultExecutionSpace(), batch_idx_range_vz,
                    KOKKOS_LAMBDA(IdxBatch_Vz const ib){
                        IdxVx const ivx(ib);
                        IdxVy const ivy(ib);
                        double dvz = (ddc::coordinate(ivx)-mean_velocity_x(ix))*a5
                            + (ddc::coordinate(ivy)-mean_velocity_y(ix))*a6;
                        for (IdxVz const iv : idx_range_vz)
                            feet_coords_vz(iv, ib) = Coord<DimVz>(ddc::coordinate(iv)+dvz);
                    });
                interpolator_vz(allfdistribu[isp][ix],
                            get_const_field(feet_coords_vz),
                            get_const_field(derivs_min_vz),
                            get_const_field(derivs_max_vz));

                // module vy 
                ddc::parallel_for_each(Kokkos::DefaultExecutionSpace(), batch_idx_range_vy,
                    KOKKOS_LAMBDA(IdxBatch_Vy const ib){
                        IdxVx const ivx(ib);
                        IdxVz const ivz(ib);
                        double dvy = (ddc::coordinate(ivx)-mean_velocity_x(ix))*a3
                            + (ddc::coordinate(ivz)-mean_velocity_z(ix))*a4;
                        for (IdxVy const iv : idx_range_vy)
                            feet_coords_vy(iv, ib) = Coord<DimVy>(ddc::coordinate(iv)+dvy);
                    });
                interpolator_vy(allfdistribu[isp][ix],
                            get_const_field(feet_coords_vy),
                            get_const_field(derivs_min_vy),
                            get_const_field(derivs_max_vy));

                // module vx 
                ddc::parallel_for_each(Kokkos::DefaultExecutionSpace(), batch_idx_range_vx,
                    KOKKOS_LAMBDA(IdxBatch_Vx const ib){
                        IdxVy const ivy(ib);
                        IdxVz const ivz(ib);
                        double dvx = (ddc::coordinate(ivy)-mean_velocity_y(ix))*a1
                            + (ddc::coordinate(ivz)-mean_velocity_z(ix))*a2;
                        for (IdxVx const iv : idx_range_vx)
                            feet_coords_vx(iv, ib) = Coord<DimVx>(ddc::coordinate(iv)+dvx);
                    });
                interpolator_vx(allfdistribu[isp][ix],
                            get_const_field(feet_coords_vx),
                            get_const_field(derivs_min_vx),
                            get_const_field(derivs_max_vx));
            
               */
            
            }

            if (count >= 2 && idx_order[0] == 0 && idx_order[1] == 2 && idx_order[2] == 1)
            {
               std::array<int, 3> order = {0, 2, 1};

                exact_splitting_solver(magnetic_field, dt, yl, y2, y3, yr, order);

                double a1 = dt*yl[1];
                double a2 = dt*yl[2];
                double a3 = dt*y2[0];
                double a4 = dt*y2[1];
                double a5 = dt*y3[0];
                double a6 = dt*y3[2];
                double a7 = dt*yr[1];
                double a8 = dt*yr[2];

                module_vx(a7,a8);
                module_vy(a5,a6);
                module_vz(a3,a4);
                module_vx(a1,a2);
            

            }

            if (count >= 2 && idx_order[0] == 1 && idx_order[1] == 0 && idx_order[2] == 2)
            {
                std::array<int, 3> order = {1, 0, 2};

                exact_splitting_solver(magnetic_field, dt, yl, y2, y3, yr, order);

                double a1 = dt*yl[0];
                double a2 = dt*yl[2];
                double a3 = dt*y2[1];
                double a4 = dt*y2[2];
                double a5 = dt*y3[0];
                double a6 = dt*y3[1];
                double a7 = dt*yr[0];
                double a8 = dt*yr[2];


                module_vy(a7,a8);
                module_vz(a5,a6);
                module_vx(a3,a4);
                module_vy(a1,a2);

               
               

            }

            if (count >= 2 && idx_order[0] == 1 && idx_order[1] == 2 && idx_order[2] == 0)
            {
                
                std::array<int, 3> order = {1, 2, 0};

                exact_splitting_solver(magnetic_field, dt, yl, y2, y3, yr, order);

                double a1 = dt*yl[0];
                double a2 = dt*yl[2];
                double a3 = dt*y2[0];
                double a4 = dt*y2[1];
                double a5 = dt*y3[1];
                double a6 = dt*y3[2];
                double a7 = dt*yr[0];
                double a8 = dt*yr[2];
                
                module_vy(a7,a8);
                module_vx(a5,a6);
                module_vz(a3,a4);
                module_vy(a1,a2);

                
            }

            if (count >= 2 && idx_order[0] == 2 && idx_order[1] == 0 && idx_order[2] == 1)
            {
                std::array<int, 3> order = {2, 0, 1};

                exact_splitting_solver(magnetic_field, dt, yl, y2, y3, yr, order);

                double a1 = dt*yl[0];
                double a2 = dt*yl[1];
                double a3 = dt*y2[1];
                double a4 = dt*y2[2];
                double a5 = dt*y3[0];
                double a6 = dt*y3[2];
                double a7 = dt*yr[0];
                double a8 = dt*yr[1];


                module_vz(a7,a8);
                module_vy(a5,a6);
                module_vx(a3,a4);
                module_vz(a1,a2);
               

            }

            if (count >= 2 && idx_order[0] == 2 && idx_order[1] == 1 && idx_order[2] == 0)
            {

                
                std::array<int, 3> order = {2, 1, 0};

                exact_splitting_solver(magnetic_field, dt, yl, y2, y3, yr, order);
                //std::cout << "The exact splitting yl is x: " << yl[0] << "y" << yl[1] << "z" << yl[2] << std::endl;
                //std::cout << "The exact splitting y2 is x: " << y2[0] << "y" << y2[1] << "z" << y2[2] << std::endl;
                //std::cout << "The exact splitting y3 is x: " << y3[0] << "y" << y3[1] << "z" << y3[2] << std::endl;
                //std::cout << "The exact splitting yr is x: " << yr[0] << "y" << yr[1] << "z" << yr[2] << std::endl;

                double a1 = dt*yl[0];
                double a2 = dt*yl[1];
                double a3 = dt*y2[0];
                double a4 = dt*y2[2];
                double a5 = dt*y3[1];
                double a6 = dt*y3[2];
                double a7 = dt*yr[0];
                double a8 = dt*yr[1];


                module_vz(a7,a8);
                module_vx(a5,a6);
                module_vy(a3,a4);
                module_vz(a1,a2);


            }

            if (count == 1 && idx_order[2] == 0)
            {
                // module vy 
                ddc::parallel_for_each(Kokkos::DefaultExecutionSpace(), batch_idx_range_vy,
                    KOKKOS_LAMBDA(IdxBatch_Vy const ib){
                        IdxVz const ivz(ib);
                        double const dvy
                                 = (ddc::coordinate(ivz) - mean_velocity_z(ix)) * std::tan(-0.5*dt*Bfield[0]);
                        for (IdxVy const iv : idx_range_vy)
                            feet_coords_vy(iv, ib) = Coord<DimVy>(ddc::coordinate(iv)+dvy);
                    });
                interpolator_vy(allfdistribu[isp][ix],
                            get_const_field(feet_coords_vy),
                            get_const_field(derivs_min_vy),
                            get_const_field(derivs_max_vy));

                // module vz 
                ddc::parallel_for_each(Kokkos::DefaultExecutionSpace(), batch_idx_range_vz,
                    KOKKOS_LAMBDA(IdxBatch_Vz const ib){
                        IdxVy const ivy(ib);
                        double const dvz
                                = (ddc::coordinate(ivy) - mean_velocity_y(ix)) * std::sin(dt*Bfield[0]);
                        for (IdxVz const iv : idx_range_vz)
                            feet_coords_vz(iv, ib) = Coord<DimVz>(ddc::coordinate(iv)+dvz);
                    });
                interpolator_vz(allfdistribu[isp][ix],
                            get_const_field(feet_coords_vz),
                            get_const_field(derivs_min_vz),
                            get_const_field(derivs_max_vz));

                // module vy 
                ddc::parallel_for_each(Kokkos::DefaultExecutionSpace(), batch_idx_range_vy,
                    KOKKOS_LAMBDA(IdxBatch_Vy const ib){
                        IdxVz const ivz(ib);
                        double const dvy
                                 = (ddc::coordinate(ivz) - mean_velocity_z(ix)) * std::tan(-0.5*dt*Bfield[0]);
                        for (IdxVy const iv : idx_range_vy)
                            feet_coords_vy(iv, ib) = Coord<DimVy>(ddc::coordinate(iv)+dvy);
                    });
                interpolator_vy(allfdistribu[isp][ix],
                            get_const_field(feet_coords_vy),
                            get_const_field(derivs_min_vy),
                            get_const_field(derivs_max_vy));

            }
            
            if (count == 1 && idx_order[2] == 1)
            {
                // module vz
                ddc::parallel_for_each(Kokkos::DefaultExecutionSpace(), batch_idx_range_vz,
                    KOKKOS_LAMBDA(IdxBatch_Vz const ib){
                        IdxVx const ivx(ib);
                        double const dvz
                                 = (ddc::coordinate(ivx) - mean_velocity_x(ix)) * std::tan(-0.5*dt*Bfield[1]);
                        for (IdxVz const iv : idx_range_vz)
                            feet_coords_vz(iv, ib) = Coord<DimVz>(ddc::coordinate(iv)+dvz);
                    });
                interpolator_vz(allfdistribu[isp][ix],
                            get_const_field(feet_coords_vz),
                            get_const_field(derivs_min_vz),
                            get_const_field(derivs_max_vz));

                // module vx
                ddc::parallel_for_each(Kokkos::DefaultExecutionSpace(), batch_idx_range_vx,
                    KOKKOS_LAMBDA(IdxBatch_Vx const ib){
                        IdxVz const ivz(ib);
                        double const dvx
                                = (ddc::coordinate(ivz) - mean_velocity_z(ix)) * std::sin(dt*Bfield[1]);
                        for (IdxVx const iv : idx_range_vx)
                            feet_coords_vx(iv, ib) = Coord<DimVx>(ddc::coordinate(iv)+dvx);
                    });
                interpolator_vx(allfdistribu[isp][ix],
                            get_const_field(feet_coords_vx),
                            get_const_field(derivs_min_vx),
                            get_const_field(derivs_max_vx));

                // module vz
                ddc::parallel_for_each(Kokkos::DefaultExecutionSpace(), batch_idx_range_vz,
                    KOKKOS_LAMBDA(IdxBatch_Vz const ib){
                        IdxVx const ivx(ib);
                        double const dvz
                                 = (ddc::coordinate(ivx) - mean_velocity_x(ix)) * std::tan(-0.5*dt*Bfield[1]);
                        for (IdxVz const iv : idx_range_vz)
                            feet_coords_vz(iv, ib) = Coord<DimVz>(ddc::coordinate(iv)+dvz);
                    });
                interpolator_vz(allfdistribu[isp][ix],
                            get_const_field(feet_coords_vz),
                            get_const_field(derivs_min_vz),
                            get_const_field(derivs_max_vz));

            }
            
            if (count == 1 && idx_order[2] == 2)
            {
                // module vx 
                ddc::parallel_for_each(Kokkos::DefaultExecutionSpace(), batch_idx_range_vx,
                    KOKKOS_LAMBDA(IdxBatch_Vx const ib){
                        IdxVy const ivy(ib);
                        double const dvx
                                 = (ddc::coordinate(ivy) - mean_velocity_y(ix)) * std::tan(-0.5*dt*Bfield[2]);
                        for (IdxVx const iv : idx_range_vx)
                            feet_coords_vx(iv, ib) = Coord<DimVx>(ddc::coordinate(iv)+dvx);
                    });
                interpolator_vx(allfdistribu[isp][ix],
                            get_const_field(feet_coords_vx),
                            get_const_field(derivs_min_vx),
                            get_const_field(derivs_max_vx));

                // module vy 
                ddc::parallel_for_each(Kokkos::DefaultExecutionSpace(), batch_idx_range_vy,
                    KOKKOS_LAMBDA(IdxBatch_Vy const ib){
                        IdxVx const ivx(ib);
                        double const dvy
                                = (ddc::coordinate(ivx) - mean_velocity_x(ix)) * std::sin(dt*Bfield[2]);
                        for (IdxVy const iv : idx_range_vy)
                            feet_coords_vy(iv, ib) = Coord<DimVy>(ddc::coordinate(iv)+dvy);
                    });
                interpolator_vy(allfdistribu[isp][ix],
                            get_const_field(feet_coords_vy),
                            get_const_field(derivs_min_vy),
                            get_const_field(derivs_max_vy));

                // module vx 
                ddc::parallel_for_each(Kokkos::DefaultExecutionSpace(), batch_idx_range_vx,
                    KOKKOS_LAMBDA(IdxBatch_Vx const ib){
                        IdxVy const ivy(ib);
                        double const dvx
                                 = (ddc::coordinate(ivy) - mean_velocity_y(ix)) * std::tan(-0.5*dt*Bfield[2]);
                        for (IdxVx const iv : idx_range_vx)
                            feet_coords_vx(iv, ib) = Coord<DimVx>(ddc::coordinate(iv)+dvx);
                    });
                interpolator_vx(allfdistribu[isp][ix],
                            get_const_field(feet_coords_vx),
                            get_const_field(derivs_min_vx),
                            get_const_field(derivs_max_vx));

            }


         //});
         });

        

         Kokkos::Profiling::popRegion();
         return allfdistribu;
     }
 };
 
