!> @ingroup splines
!> @brief   Interpolator for 1D splines of arbitrary degree,
!>          on uniform and non-uniform grids
!> @author  Yaman Güçlü  - IPP Garching
!> @author  Edoardo Zoni - IPP Garching

module sll_m_spline_interpolator_1d
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "logging.inc"
  use iso_C_binding, only: c_ptr, c_double, c_int, c_bool, c_loc, c_null_ptr

  use sll_m_bsplines, only: sll_bsplines_t
  use sll_m_spline_1d, only: sll_spline_1d_t

  implicit none

  public :: &
    sll_spline_interpolator_1d_t, &
    sll_spline_1d_compute_num_cells, &
    sll_spline_interpolator_1d__new, &
    sll_spline_interpolator_1d__free, &
    sll_spline_interpolator_1d__get_interp_points, &
    sll_spline_interpolator_1d__compute_interpolant

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type, bind(c) :: sll_spline_interpolator_1d_t
     type(c_ptr) :: spl_interp_obj
     integer(c_int) :: bc_xmin
     integer(c_int) :: bc_xmax
  end type

  !-----------------------------------------------------------------------------
  !> @brief      Calculate number of cells from number of interpolation points
  !> @details    Important for parallelization: for given spline degree and BCs,
  !>             calculate the number of grid cells that yields the desired
  !>             number of interpolation points
  !>
  !> @param[in]  degree   spline degree
  !> @param[in]  bc_xmin  boundary condition type at left  boundary (x=xmin)
  !> @param[in]  bc_xmax  boundary condition type at right boundary (x=xmax)
  !> @param[in]  nipts    desired number of interpolation points
  !> @param[out] ncells   calculated number of cells in domain
  !-----------------------------------------------------------------------------
  interface
      function sll_spline_1d_compute_num_cells( &
          degree , &
          bc_xmin, &
          bc_xmax, &
          nipts  ) bind(C, name="compute_num_cells")

        use iso_C_binding, only: c_ptr, c_double, c_int, c_bool
        integer(c_int), value , intent(in   ) :: degree
        integer(c_int), value , intent(in   ) :: bc_xmin
        integer(c_int), value , intent(in   ) :: bc_xmax
        integer(c_int), value , intent(in   ) :: nipts
        integer(c_int) :: sll_spline_1d_compute_num_cells
      end function
  end interface

  interface
      function new_spline_interpolator_1d( &
      bspl    , &
      xmin_bc , &
      xmax_bc )  bind(C)

         use iso_C_binding, only: c_ptr, c_double, c_int, c_bool
         type(c_ptr) :: new_spline_interpolator_1d
         type(c_ptr)   , value , intent(in)  :: bspl
         integer(c_int), value , intent(in)  :: xmin_bc
         integer(c_int), value , intent(in)  :: xmax_bc
      end function new_spline_interpolator_1d
  end interface

  interface
      subroutine free_spline_interpolator_1d( &
              spl_interp) bind (C)

          use iso_C_binding, only: c_ptr, c_double, c_int, c_bool
          type(c_ptr), value, intent(in) :: spl_interp
      end subroutine
  end interface

  interface
      subroutine get_interp_points( &
              spl_interp, &
              tau       , &
              npts      ) bind (C)

          use iso_C_binding, only: c_ptr, c_double, c_int, c_bool
          type(c_ptr), value          , intent(in   ) :: spl_interp
          integer(c_int), value       , intent(in   ) :: npts
          real(c_double)              , intent(inout) :: tau(0:npts-1)
      end subroutine
  end interface

  interface
      function get_n_interp_points( &
              spl_interp ) bind (C)

          use iso_C_binding, only: c_ptr, c_double, c_int, c_bool
          type(c_ptr), value          , intent(in   ) :: spl_interp
          integer(c_int)              :: get_n_interp_points
      end function
  end interface

  interface
      subroutine compute_interpolant( &
              spl_interp  , &
              spline      , &
              vals        , &
              nvals       , &
              derivs_xmin , &
              nderivs_xmin, &
              derivs_xmax , &
              nderivs_xmax ) bind (C)

          use iso_C_binding, only: c_ptr, c_double, c_int, c_bool
          type(c_ptr), value         , intent(in   ) :: spl_interp
          type(c_ptr), value         , intent(in   ) :: spline
          integer(c_int), value      , intent(in   ) :: nvals
          real(c_double)             , intent(in   ) :: vals(0:nvals-1)
          type(c_ptr), value         , intent(in   ) :: derivs_xmin
          integer(c_int), value      , intent(in   ) :: nderivs_xmin
          type(c_ptr), value         , intent(in   ) :: derivs_xmax
          integer(c_int), value      , intent(in   ) :: nderivs_xmax
      end subroutine
  end interface


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


  !-----------------------------------------------------------------------------
  !> @brief      Initialize a 1D spline interpolator object
  !> @param[out] self     1D spline interpolator
  !> @param[in]  bspl     B-splines (basis)
  !> @param[in]  bc_xmin  boundary condition at xmin
  !> @param[in]  bc_xmax  boundary condition at xmax
  !-----------------------------------------------------------------------------
  function sll_spline_interpolator_1d__new( bspl, bc_xmin, bc_xmax ) result( self )

    type(sll_spline_interpolator_1d_t)                 :: self
    type(sll_bsplines_t),       target , intent(in   ) :: bspl
    integer(c_int)                     , intent(in   ) :: bc_xmin
    integer(c_int)                     , intent(in   ) :: bc_xmax

    self % bc_xmin = bc_xmin
    self % bc_xmax = bc_xmax

    self % spl_interp_obj = new_spline_interpolator_1d(bspl % bspline_obj, bc_xmin, bc_xmax)

  end function sll_spline_interpolator_1d__new

  !-----------------------------------------------------------------------------
  !> @brief        Destroy local objects and free allocated memory
  !> @param[inout] self  1D spline interpolator
  !-----------------------------------------------------------------------------
  subroutine sll_spline_interpolator_1d__free( self )

    type(sll_spline_interpolator_1d_t), intent(inout) :: self

    call free_spline_interpolator_1d(self % spl_interp_obj)
  end subroutine sll_spline_interpolator_1d__free

  !-----------------------------------------------------------------------------
  !> @brief      Get coordinates of interpolation points (1D grid)
  !> @param[in]  self  1D spline interpolator
  !> @param[out] tau   x coordinates of interpolation points
  !-----------------------------------------------------------------------------
  subroutine sll_spline_interpolator_1d__get_interp_points( self, tau )

    type(sll_spline_interpolator_1d_t) , intent(in   ) :: self
    real(c_double), allocatable        , intent(inout) :: tau(:)

    integer(c_int) :: n_tau

    n_tau = get_n_interp_points( self % spl_interp_obj )

    allocate( tau(n_tau) )

    call get_interp_points( self % spl_interp_obj, tau, n_tau )

  end subroutine sll_spline_interpolator_1d__get_interp_points

  !-----------------------------------------------------------------------------
  !> @brief        Compute interpolating 1D spline
  !> @details      Compute coefficients of 1D spline that interpolates function
  !>               values on grid. If Hermite BCs are used, function derivatives
  !>               at appropriate boundary are also needed.
  !>
  !> @param[in]    self        1D spline interpolator
  !> @param[inout] spline      1D spline
  !> @param[in]    gtau        function values of interpolation points
  !> @param[in]    derivs_xmin (optional) array with boundary conditions at xmin
  !> @param[in]    derivs_xmax (optional) array with boundary conditions at xmax
  !-----------------------------------------------------------------------------
  subroutine sll_spline_interpolator_1d__compute_interpolant( self, &
      spline, gtau, derivs_xmin, derivs_xmax )

    type (sll_spline_interpolator_1d_t), intent(in   ) :: self
    type (sll_spline_1d_t)             , intent(inout) :: spline
    real(c_double)                           , intent(in   ) :: gtau(:)
    real(c_double),                  optional, intent(in   ) :: derivs_xmin(:)
    real(c_double),                  optional, intent(in   ) :: derivs_xmax(:)

    real(c_double), allocatable, target :: derivs_xmin_a(:)
    real(c_double), allocatable, target :: derivs_xmax_a(:)

    if ( present(derivs_xmin) ) then
        allocate( derivs_xmin_a(0:size(derivs_xmin,1)-1), source=derivs_xmin )
        if ( present(derivs_xmax) ) then
            allocate( derivs_xmax_a(0:size(derivs_xmax,1)-1), source=derivs_xmax )
            call compute_interpolant(self % spl_interp_obj, spline % spline_obj, gtau, size(gtau,1), &
                c_loc(derivs_xmin_a), size(derivs_xmin,1), c_loc(derivs_xmax_a), size(derivs_xmax,1))
        else
            call compute_interpolant(self % spl_interp_obj, spline % spline_obj, gtau, size(gtau,1), &
                c_loc(derivs_xmin_a), size(derivs_xmin,1), c_null_ptr, 0)
        end if
    else
        if ( present(derivs_xmax) ) then
            allocate( derivs_xmax_a(0:size(derivs_xmax,1)-1), source=derivs_xmax )
            call compute_interpolant(self % spl_interp_obj, spline % spline_obj, gtau, size(gtau,1), &
                c_null_ptr, 0, c_loc(derivs_xmax_a), size(derivs_xmax,1))
        else
            call compute_interpolant(self % spl_interp_obj, spline % spline_obj, gtau, size(gtau,1), &
                c_null_ptr, 0, c_null_ptr, 0)
        end if
    end if


  end subroutine sll_spline_interpolator_1d__compute_interpolant


end module sll_m_spline_interpolator_1d
