!> @ingroup splines
!> @brief   Interpolator for 2D tensor-product splines of arbitrary degree,
!>          on uniform and non-uniform grids (directions are independent)
!> @author  Yaman Güçlü  - IPP Garching
!> @author  Edoardo Zoni - IPP Garching

module sll_m_spline_interpolator_2d
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "logging.inc"
  use iso_C_binding, only: c_ptr, c_double, c_int, c_bool, c_loc, c_null_ptr

  use prec_const, only: F64

  use sll_m_bsplines, only: sll_bsplines_t
  use sll_m_spline_2d, only: sll_spline_2d_t

  implicit none

  public :: &
    sll_spline_interpolator_2d_t, &
    sll_spline_2d_boundary_data, &
    sll_spline_2d_compute_num_cells, &
    sll_spline_interpolator_2d__new, &
    sll_spline_interpolator_2d__free, &
    sll_spline_interpolator_2d__get_interp_points, &
    sll_spline_interpolator_2d__compute_interpolant

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !-----------------------------------------------------------------------------
  !> Container for 2D boundary condition data:
  !>  . x1-derivatives at x1_min and x1_max, for all values of x2;
  !>  . x2-derivatives at x2_min and x2_max, for all values of x1;
  !>  . mixed derivatives at the four corners a,b,c,d.
  !>
  !>  x2_max  ____________
  !>         |            |
  !>         | c        d |
  !>         |            |
  !>         |            |
  !>         |            |
  !>         | a        b |
  !>  x2_min |____________|
  !>       x1_min       x1_max
  !>
  !-----------------------------------------------------------------------------
  type :: sll_spline_2d_boundary_data
    real(F64), allocatable :: derivs_x1_min (:,:)
    real(F64), allocatable :: derivs_x1_max (:,:)
    real(F64), allocatable :: derivs_x2_min (:,:)
    real(F64), allocatable :: derivs_x2_max (:,:)
    real(F64), allocatable :: mixed_derivs_a(:,:)
    real(F64), allocatable :: mixed_derivs_b(:,:)
    real(F64), allocatable :: mixed_derivs_c(:,:)
    real(F64), allocatable :: mixed_derivs_d(:,:)
  end type sll_spline_2d_boundary_data

  !-----------------------------------------------------------------------------
  !> 2D tensor-product spline interpolator
  !-----------------------------------------------------------------------------
  type :: sll_spline_interpolator_2d_t
      type(c_ptr) :: spl_interp_obj
  end type sll_spline_interpolator_2d_t

  interface
      subroutine compute_num_cells_2d( &
              degree_1   , &
              bc_xmin_1  , &
              bc_xmax_1  , &
              nipts_1    , &
              degree_2   , &
              bc_xmin_2  , &
              bc_xmax_2  , &
              nipts_2    , &
              ncell_1    , &
              ncell_2    ) bind(C)
          use iso_C_binding, only: c_int
          integer(c_int), value , intent(in   ) ::degree_1
          integer(c_int), value , intent(in   ) ::bc_xmin_1
          integer(c_int), value , intent(in   ) ::bc_xmax_1
          integer(c_int), value , intent(in   ) ::nipts_1
          integer(c_int), value , intent(in   ) ::degree_2
          integer(c_int), value , intent(in   ) ::bc_xmin_2
          integer(c_int), value , intent(in   ) ::bc_xmax_2
          integer(c_int), value , intent(in   ) ::nipts_2
          integer(c_int),         intent(  out) ::ncell_1
          integer(c_int),         intent(  out) ::ncell_2
      end subroutine compute_num_cells_2d
  end interface

  interface
      function new_spline_interpolator_2d( &
              bspl_1       , &
              bc_xmin_1    , &
              bc_xmax_1    , &
              bspl_2       , &
              bc_xmin_2    , &
              bc_xmax_2    ) bind(C)
          use iso_C_binding, only: c_int, c_ptr
          type(c_ptr) :: new_spline_interpolator_2d
          type(c_ptr)   , value , intent(in   ) :: bspl_1
          integer(c_int), value , intent(in   ) ::bc_xmin_1
          integer(c_int), value , intent(in   ) ::bc_xmax_1
          type(c_ptr)   , value , intent(in   ) :: bspl_2
          integer(c_int), value , intent(in   ) ::bc_xmin_2
          integer(c_int), value , intent(in   ) ::bc_xmax_2
      end function new_spline_interpolator_2d
  end interface

  interface
      subroutine free_spline_interpolator_2d( &
              spl_interp) bind (C)

          use iso_C_binding, only: c_ptr
          type(c_ptr), value, intent(in) :: spl_interp
      end subroutine
  end interface

  interface 
      subroutine get_interp_points_2d( &
              spl_interp       , &
              interp_points_1  , &
              npts_1           , &
              interp_points_2  , &
              npts_2           ) bind(C)
          use iso_C_binding, only: c_int, c_ptr, c_double
          type(c_ptr)   , value, intent(in   ) :: spl_interp
          integer(c_int), value, intent(in   ) :: npts_1
          real(c_double),        intent(inout) :: interp_points_1(0:npts_1-1)
          integer(c_int), value, intent(in   ) :: npts_2
          real(c_double),        intent(inout) :: interp_points_2(0:npts_2-1)
      end subroutine get_interp_points_2d
  end interface

  interface 
      subroutine get_n_interp_points_2d( &
              spl_interp       , &
              npts_1           , &
              npts_2           ) bind(C)
          use iso_C_binding, only: c_int, c_ptr, c_double
          type(c_ptr)   , value, intent(in   ) :: spl_interp
          integer(c_int),        intent(  out) :: npts_1
          integer(c_int),        intent(  out) :: npts_2
      end subroutine get_n_interp_points_2d
  end interface

  interface
      subroutine compute_interpolant_2d( &
              spl_interp    , &
              spline        , &
              vals          , nvals_1    , nvals_2    , &
              derivs_x1_min , nd_x1_min_1, nd_x1_min_2, &
              derivs_x1_max , nd_x1_max_1, nd_x1_max_2, &
              derivs_x2_min , nd_x2_min_1, nd_x2_min_2, &
              derivs_x2_max , nd_x2_max_1, nd_x2_max_2, &
              mixed_derivs_a, nmd_a_1    , nmd_a_2    , &
              mixed_derivs_b, nmd_b_1    , nmd_b_2    , &
              mixed_derivs_c, nmd_c_1    , nmd_c_2    , &
              mixed_derivs_d, nmd_d_1    , nmd_d_2    ) bind(C)
          use iso_C_binding, only: c_int, c_ptr, c_double
          type(c_ptr)    , value, intent(in   ) :: spl_interp
          type(c_ptr)    , value, intent(in   ) :: spline
          integer(c_int) , value, intent(in   ) :: nvals_1, nvals_2
          real(c_double)        , intent(in   ) :: vals(0:nvals_1-1,0:nvals_2-1)
          type(c_ptr)    , value, intent(in   ) :: derivs_x1_min
          integer(c_int) , value, intent(in   ) :: nd_x1_min_1, nd_x1_min_2
          type(c_ptr)    , value, intent(in   ) :: derivs_x1_max
          integer(c_int) , value, intent(in   ) :: nd_x1_max_1, nd_x1_max_2
          type(c_ptr)    , value, intent(in   ) :: derivs_x2_min
          integer(c_int) , value, intent(in   ) :: nd_x2_min_1, nd_x2_min_2
          type(c_ptr)    , value, intent(in   ) :: derivs_x2_max
          integer(c_int) , value, intent(in   ) :: nd_x2_max_1, nd_x2_max_2
          type(c_ptr)    , value, intent(in   ) :: mixed_derivs_a
          integer(c_int) , value, intent(in   ) :: nmd_a_1, nmd_a_2
          type(c_ptr)    , value, intent(in   ) :: mixed_derivs_b
          integer(c_int) , value, intent(in   ) :: nmd_b_1, nmd_b_2
          type(c_ptr)    , value, intent(in   ) :: mixed_derivs_c
          integer(c_int) , value, intent(in   ) :: nmd_c_1, nmd_c_2
          type(c_ptr)    , value, intent(in   ) :: mixed_derivs_d
          integer(c_int) , value, intent(in   ) :: nmd_d_1, nmd_d_2
      end subroutine compute_interpolant_2d
  end interface

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !-----------------------------------------------------------------------------
  !> @brief      Calculate number of cells from number of interpolation points
  !> @details    Important for parallelization: for given spline degree and BCs,
  !>             calculate the numbers of grid cells along x1 and x2 that yield
  !>             the desired number of interpolation points along x1 and x2
  !>
  !> @param[in]  degree   spline degrees along x1 and x2
  !> @param[in]  bc_xmin  boundary conditions at x1_min and x2_min
  !> @param[in]  bc_xmax  boundary conditions at x1_max and x2_max
  !> @param[in]  nipts    desired number of interpolation points along x1 and x2
  !> @param[out] ncells   calculated number of grid cells along x1 and x2
  !-----------------------------------------------------------------------------
  subroutine sll_spline_2d_compute_num_cells( &
      degree , &
      bc_xmin, &
      bc_xmax, &
      nipts  , &
      ncells )

    integer, intent(in   ) :: degree (2)
    integer, intent(in   ) :: bc_xmin(2)
    integer, intent(in   ) :: bc_xmax(2)
    integer, intent(in   ) :: nipts  (2)
    integer, intent(  out) :: ncells (2)

    call compute_num_cells_2d( degree(1), bc_xmin(1), bc_xmax(1), nipts(1), &
        degree(1), bc_xmin(1), bc_xmax(1), nipts(1), ncells(1), ncells(2))

  end subroutine sll_spline_2d_compute_num_cells

  !-----------------------------------------------------------------------------
  !> @brief      Initialize a 2D tensor-product spline interpolator object
  !> @param[out] self     2D tensor-product spline interpolator
  !> @param[in]  bspl1    B-splines (basis) along x1 direction
  !> @param[in]  bspl2    B-splines (basis) along x2 direction
  !> @param[in]  bc_xmin  boundary conditions at x1_min and x2_min
  !> @param[in]  bc_xmax  boundary conditions at x1_max and x2_max
  !-----------------------------------------------------------------------------
  function sll_spline_interpolator_2d__new( &
    bspl1  , &
    bspl2  , &
    bc_xmin, &
    bc_xmax ) result( self )

    type(sll_spline_interpolator_2d_t)                 :: self
    type(sll_bsplines_t)               , intent(in   ) :: bspl1
    type(sll_bsplines_t)               , intent(in   ) :: bspl2
    integer                            , intent(in   ) :: bc_xmin(2)
    integer                            , intent(in   ) :: bc_xmax(2)

    self % spl_interp_obj = new_spline_interpolator_2d( bspl1 % bspline_obj, bc_xmin(1), bc_xmax(1), &
        bspl2 % bspline_obj, bc_xmin(2), bc_xmax(2))

  end function sll_spline_interpolator_2d__new

  !-----------------------------------------------------------------------------
  !> @brief        Destroy local objects and free allocated memory
  !> @param[inout] self  2D tensor-product spline interpolator
  !-----------------------------------------------------------------------------
  subroutine sll_spline_interpolator_2d__free( self )

    type(sll_spline_interpolator_2d_t), intent(inout) :: self

    call free_spline_interpolator_2d( self % spl_interp_obj )

  end subroutine sll_spline_interpolator_2d__free

  !-----------------------------------------------------------------------------
  !> @brief      Get coordinates of interpolation points (2D tensor grid)
  !> @param[in]  self  2D tensor-product spline interpolator
  !> @param[out] tau1  x1 coordinates of interpolation points
  !> @param[out] tau2  x2 coordinates of interpolation points
  !-----------------------------------------------------------------------------
  subroutine sll_spline_interpolator_2d__get_interp_points( self, tau1, tau2 )

    type(sll_spline_interpolator_2d_t), intent(in   ) :: self
    real(F64),               allocatable, intent(  out) :: tau1(:)
    real(F64),               allocatable, intent(  out) :: tau2(:)

    integer(c_int) :: n1,n2

    call get_n_interp_points_2d( self % spl_interp_obj, n1, n2 )

    allocate( tau1(n1) )
    allocate( tau2(n2) )

    call get_interp_points_2d( self % spl_interp_obj, tau1, n1, tau2, n2 )

  end subroutine sll_spline_interpolator_2d__get_interp_points

  !-----------------------------------------------------------------------------
  !> @brief        Compute interpolating 2D spline
  !> @details      Compute coefficients of 2D tensor-product spline that
  !>               interpolates function values on grid. If Hermite BCs are used,
  !>               function derivatives at appropriate boundaries are also needed.
  !>
  !> @param[inout] self           2D tensor-product spline interpolator
  !> @param[inout] spline         2D tensor-product spline
  !> @param[in]    gtau           function values of interpolation points
  !> @param[in]    boundary_data  (optional) structure with boundary conditions
  !-----------------------------------------------------------------------------
  subroutine sll_spline_interpolator_2d__compute_interpolant( self, &
      spline, gtau, boundary_data )

    type(sll_spline_interpolator_2d_t), intent(inout)           :: self
    type(sll_spline_2d_t)             , intent(inout)           :: spline
    real(F64)                         , intent(in   )           :: gtau(:,:)
    type(sll_spline_2d_boundary_data), target, intent(in   ), optional :: boundary_data

    type(c_ptr)    :: derivs_x1_min 
    integer(c_int) :: n_derivs_x1_min_1
    integer(c_int) :: n_derivs_x1_min_2
    type(c_ptr)    :: derivs_x1_max 
    integer(c_int) :: n_derivs_x1_max_1
    integer(c_int) :: n_derivs_x1_max_2
    type(c_ptr)    :: derivs_x2_min 
    integer(c_int) :: n_derivs_x2_min_1
    integer(c_int) :: n_derivs_x2_min_2
    type(c_ptr)    :: derivs_x2_max 
    integer(c_int) :: n_derivs_x2_max_1
    integer(c_int) :: n_derivs_x2_max_2
    type(c_ptr)    :: mixed_derivs_a
    integer(c_int) :: n_mixed_derivs_a_1
    integer(c_int) :: n_mixed_derivs_a_2
    type(c_ptr)    :: mixed_derivs_b
    integer(c_int) :: n_mixed_derivs_b_1
    integer(c_int) :: n_mixed_derivs_b_2
    type(c_ptr)    :: mixed_derivs_c
    integer(c_int) :: n_mixed_derivs_c_1
    integer(c_int) :: n_mixed_derivs_c_2
    type(c_ptr)    :: mixed_derivs_d
    integer(c_int) :: n_mixed_derivs_d_1
    integer(c_int) :: n_mixed_derivs_d_2

    character(len=*), parameter :: &
      this_sub_name = "sll_spline_interpolator_2d_t % compute_interpolant"

   if ( .not. present(boundary_data) ) then
       call compute_interpolant_2d( self % spl_interp_obj, spline % spline_obj, &
           gtau, size(gtau,1), size(gtau,2), &
           c_null_ptr, 0, 0, c_null_ptr, 0, 0, &
           c_null_ptr, 0, 0, c_null_ptr, 0, 0, &
           c_null_ptr, 0, 0, c_null_ptr, 0, 0, &
           c_null_ptr, 0, 0, c_null_ptr, 0, 0)
   else
       if ( allocated( boundary_data % derivs_x1_min ) ) then
           derivs_x1_min = c_loc(boundary_data % derivs_x1_min)
           n_derivs_x1_min_1 = size(boundary_data % derivs_x1_min, 1)
           n_derivs_x1_min_2 = size(boundary_data % derivs_x1_min, 2)
       else
           derivs_x1_min = c_null_ptr
           n_derivs_x1_min_1 = 0
           n_derivs_x1_min_2 = 0
       endif
       if ( allocated( boundary_data % derivs_x1_max ) ) then
           derivs_x1_max = c_loc(boundary_data % derivs_x1_max)
           n_derivs_x1_max_1 = size(boundary_data % derivs_x1_max, 1)
           n_derivs_x1_max_2 = size(boundary_data % derivs_x1_max, 2)
       else
           derivs_x1_max = c_null_ptr
           n_derivs_x1_max_1 = 0
           n_derivs_x1_max_2 = 0
       endif
       if ( allocated( boundary_data % derivs_x2_min ) ) then
           derivs_x2_min = c_loc(boundary_data % derivs_x2_min)
           n_derivs_x2_min_1 = size(boundary_data % derivs_x2_min, 1)
           n_derivs_x2_min_2 = size(boundary_data % derivs_x2_min, 2)
       else
           derivs_x2_min = c_null_ptr
           n_derivs_x2_min_1 = 0
           n_derivs_x2_min_2 = 0
       endif
       if ( allocated( boundary_data % derivs_x2_max ) ) then
           derivs_x2_max = c_loc(boundary_data % derivs_x2_max)
           n_derivs_x2_max_1 = size(boundary_data % derivs_x2_max, 1)
           n_derivs_x2_max_2 = size(boundary_data % derivs_x2_max, 2)
       else
           derivs_x2_max = c_null_ptr
           n_derivs_x2_max_1 = 0
           n_derivs_x2_max_2 = 0
       endif
       if ( allocated( boundary_data % mixed_derivs_a ) ) then
           mixed_derivs_a = c_loc(boundary_data % mixed_derivs_a)
           n_mixed_derivs_a_1 = size(boundary_data % mixed_derivs_a, 1)
           n_mixed_derivs_a_2 = size(boundary_data % mixed_derivs_a, 2)
       else
           mixed_derivs_a = c_null_ptr
           n_mixed_derivs_a_1 = 0
           n_mixed_derivs_a_2 = 0
       endif
       if ( allocated( boundary_data % mixed_derivs_b ) ) then
           mixed_derivs_b = c_loc(boundary_data % mixed_derivs_b)
           n_mixed_derivs_b_1 = size(boundary_data % mixed_derivs_b, 1)
           n_mixed_derivs_b_2 = size(boundary_data % mixed_derivs_b, 2)
       else
           mixed_derivs_b = c_null_ptr
           n_mixed_derivs_b_1 = 0
           n_mixed_derivs_b_2 = 0
       endif
       if ( allocated( boundary_data % mixed_derivs_c ) ) then
           mixed_derivs_c = c_loc(boundary_data % mixed_derivs_c)
           n_mixed_derivs_c_1 = size(boundary_data % mixed_derivs_c, 1)
           n_mixed_derivs_c_2 = size(boundary_data % mixed_derivs_c, 2)
       else
           mixed_derivs_c = c_null_ptr
           n_mixed_derivs_c_1 = 0
           n_mixed_derivs_c_2 = 0
       endif
       if ( allocated( boundary_data % mixed_derivs_d ) ) then
           mixed_derivs_d = c_loc(boundary_data % mixed_derivs_d)
           n_mixed_derivs_d_1 = size(boundary_data % mixed_derivs_d, 1)
           n_mixed_derivs_d_2 = size(boundary_data % mixed_derivs_d, 2)
       else
           mixed_derivs_d = c_null_ptr
           n_mixed_derivs_d_1 = 0
           n_mixed_derivs_d_2 = 0
       endif

       call compute_interpolant_2d( self % spl_interp_obj, spline % spline_obj, gtau, size(gtau,1), size(gtau,2), &
           derivs_x1_min , n_derivs_x1_min_1 , n_derivs_x1_min_2, &
           derivs_x1_max , n_derivs_x1_max_1 , n_derivs_x1_max_2, &
           derivs_x2_min , n_derivs_x2_min_1 , n_derivs_x2_min_2, &
           derivs_x2_max , n_derivs_x2_max_1 , n_derivs_x2_max_2, &
           mixed_derivs_a, n_mixed_derivs_a_1, n_mixed_derivs_a_2, &
           mixed_derivs_b, n_mixed_derivs_b_1, n_mixed_derivs_b_2, &
           mixed_derivs_c, n_mixed_derivs_c_1, n_mixed_derivs_c_2, &
           mixed_derivs_d, n_mixed_derivs_d_1, n_mixed_derivs_d_2  )
    endif
  end subroutine sll_spline_interpolator_2d__compute_interpolant

end module sll_m_spline_interpolator_2d
