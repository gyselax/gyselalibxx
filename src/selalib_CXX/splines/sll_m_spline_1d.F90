!> @ingroup splines
!> @brief   Module for 1D splines, linear combination of B-spline functions.
!>
!> @details Here we define a 1D spline type as an element of the linear space
!>          given by the span of the given B-splines (basis functions).
!>          Therefore, initialization of a 1D spline object requires an existing
!>          B-splines object, to which a private (polymorphic) pointer is
!>          associated.
!>          The B-spline coefficients are stored in a public allocatable array;
!>          at initialization the array is allocated to the proper size and all
!>          values are set to zero.
!>          In most situations the B-spline coefficients are not set directly by
!>          the end user, but are computed by some other object (e.g., a Poisson
!>          solver or a spline interpolator).
!>          Various public methods allow the user to evaluate the 1D spline S(x)
!>          and its derivative ∂S(x)/∂x any position x.
!>
!> @author  Yaman Güçlü  - IPP Garching
!> @author  Edoardo Zoni - IPP Garching

module sll_m_spline_1d
#include "logging.inc"
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  use iso_C_binding, only: c_ptr, c_double, c_null_ptr, c_f_pointer

  use sll_m_bsplines    , only: sll_bsplines_t

  implicit none

  public :: &
    sll_spline_1d_t, &
    sll_spline_1d__new, &
    sll_spline_1d__free, &
    sll_spline_1d__eval, &
    sll_spline_1d__eval_deriv, &
    sll_spline_1d__eval_array, &
    sll_spline_1d__eval_array_deriv, &
    sll_spline_1d__integrate

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !> 1D spline
  type :: sll_spline_1d_t
      type(sll_bsplines_t), pointer :: bspl
      real(c_double), pointer, contiguous :: bcoef(:)
      type(c_ptr) :: spline_obj
  end type sll_spline_1d_t

  interface
      function spline_1d_new( &
              bspl ) bind(C)
          use iso_C_binding, only: c_ptr
          type(c_ptr), value, intent(in) :: bspl
          type(c_ptr)                    :: spline_1d_new
      end function
  end interface

  interface
      subroutine spline_1d_free( &
              spl ) bind(C)
          use iso_C_binding, only: c_ptr
          type(c_ptr), value, intent(in   ) :: spl
      end subroutine
  end interface

  interface
      LOG_PURE function spline_1d_eval( &
              spl, &
              x ) bind(C)
          use iso_C_binding, only: c_ptr, c_double
          type(c_ptr), value,    intent(in) :: spl
          real(c_double), value, intent(in) :: x
          real(c_double)                    :: spline_1d_eval
      end function
  end interface

  interface
      LOG_PURE function spline_1d_eval_deriv( &
              spl, &
              x ) bind(C)
          use iso_C_binding, only: c_ptr, c_double
          type(c_ptr), value,    intent(in) :: spl
          real(c_double), value, intent(in) :: x
          real(c_double)                    :: spline_1d_eval_deriv
      end function
  end interface

  interface
      LOG_PURE subroutine spline_1d_eval_array( &
              spl   , &
              x     , &
              nx    , &
              y     , &
              ny       ) bind(C)
          use iso_C_binding, only: c_ptr, c_double, c_int
          type(c_ptr), value      , intent(in   ) :: spl
          integer(c_int), value   , intent(in   ) :: nx
          real(c_double)          , intent(in   ) :: x(0:nx-1)
          integer(c_int), value   , intent(in   ) :: ny
          real(c_double)          , intent(inout) :: y(0:ny-1)
      end subroutine
  end interface

  interface
      LOG_PURE subroutine spline_1d_eval_array_deriv( &
              spl   , &
              x     , &
              nx    , &
              y     , &
              ny       ) bind(C)
          use iso_C_binding, only: c_ptr, c_double, c_int
          type(c_ptr), value      , intent(in   ) :: spl
          integer(c_int), value   , intent(in   ) :: nx
          real(c_double)          , intent(in   ) :: x(0:nx-1)
          integer(c_int), value   , intent(in   ) :: ny
          real(c_double)          , intent(inout) :: y(0:ny-1)
      end subroutine
  end interface

  interface
      LOG_PURE function spline_1d_integrate( &
              spl    ) bind(C)
          use iso_C_binding, only: c_ptr, c_double
          type(c_ptr), value      , intent(in   ) :: spl
          real(c_double)                          :: spline_1d_integrate
      end function
  end interface

  interface
      LOG_PURE function spline_1d_get_ncoeffs( &
              spl    ) bind(C)
          use iso_C_binding, only: c_ptr, c_int
          type(c_ptr), value      , intent(in   ) :: spl
          integer(c_int)                          :: spline_1d_get_ncoeffs
      end function
  end interface

  interface
      LOG_PURE function spline_1d_get_coeffs( &
              spl    ) bind(C)
          use iso_C_binding, only: c_ptr
          type(c_ptr), value      , intent(in   ) :: spl
          type(c_ptr)                             :: spline_1d_get_coeffs
      end function
  end interface

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !-----------------------------------------------------------------------------
  !> @brief      Initialize 1D spline object as element of span(B-splines)
  !> @param[out] self      1D spline: new element of 1D spline space
  !> @param[in]  bspl  B-splines: given basis of 1D spline space
  !-----------------------------------------------------------------------------
  function sll_spline_1d__new( bspl ) result( self )


    type(sll_spline_1d_t)                        :: self
    type(sll_bsplines_t ), intent(in   ), target :: bspl

    type(c_ptr) :: coeff_ptr
    integer :: n

    self % bspl => bspl

    self % spline_obj = spline_1d_new(bspl % bspline_obj)

    n = spline_1d_get_ncoeffs( self % spline_obj )

    coeff_ptr = spline_1d_get_coeffs( self % spline_obj )

    call c_f_pointer( coeff_ptr, self % bcoef, shape=[n] )

  end function sll_spline_1d__new

  !-----------------------------------------------------------------------------
  !> @brief        Destroy 1D spline (re-initialization is possible afterwards)
  !> @param[inout] self  1D spline
  !-----------------------------------------------------------------------------
  subroutine sll_spline_1d__free( self )

    type(sll_spline_1d_t), intent(inout) :: self

    call spline_1d_free(self % spline_obj)
    self % spline_obj = c_null_ptr

  end subroutine sll_spline_1d__free

  !-----------------------------------------------------------------------------
  !> @brief     Evaluate value of 1D spline at location x: y=S(x)
  !> @param[in] self  1D spline
  !> @param[in] x     evaluation point
  !> @returns         spline value y=S(x)
  !-----------------------------------------------------------------------------
  LOG_PURE function sll_spline_1d__eval( self, x ) result( y )

    type(sll_spline_1d_t), intent(in) :: self
    real(c_double)              , intent(in) :: x
    real(c_double) :: y

    y = spline_1d_eval(self % spline_obj, x)

  end function sll_spline_1d__eval

  !-----------------------------------------------------------------------------
  !> @brief     Evaluate derivative of 1D spline at location x: y=S'(x)
  !> @param[in] self  1D spline
  !> @param[in] x     evaluation point
  !> @returns         spline derivative y=S'(x)
  !-----------------------------------------------------------------------------
  LOG_PURE function sll_spline_1d__eval_deriv( self, x ) result( y )

    type(sll_spline_1d_t), intent(in) :: self
    real(c_double)              , intent(in) :: x
    real(c_double) :: y

    y = spline_1d_eval_deriv(self % spline_obj, x)

  end function sll_spline_1d__eval_deriv

  !-----------------------------------------------------------------------------
  !> @brief      Evaluate value of 1D spline at all locations in array x
  !> @param[in]  self  1D spline
  !> @param[in]  x     array of evaluation points x[i]
  !> @param[out] y     array of spline values y[i]=S(x[i])
  !-----------------------------------------------------------------------------
  LOG_PURE subroutine sll_spline_1d__eval_array( self, x, y )

    type(sll_spline_1d_t), intent(in   ) :: self
    real(c_double)              , intent(in   ) :: x(:)
    real(c_double)              , intent(  out) :: y(:)

    call spline_1d_eval_array(self % spline_obj, x, size(x,1), y, size(y,1))

  end subroutine sll_spline_1d__eval_array

  !-----------------------------------------------------------------------------
  !> @brief      Evaluate derivative of 1D spline at all locations in array x
  !> @param[in]  self  1D spline
  !> @param[in]  x     array of evaluation points x[i]
  !> @param[out] y     array of spline derivatives y[i]=S'(x[i])
  !-----------------------------------------------------------------------------
  LOG_PURE subroutine sll_spline_1d__eval_array_deriv( self, x, y )

    type(sll_spline_1d_t), intent(in   ) :: self
    real(c_double)              , intent(in   ) :: x(:)
    real(c_double)              , intent(  out) :: y(:)

    call spline_1d_eval_array_deriv(self % spline_obj, x, size(x,1), y, size(y,1))

  end subroutine sll_spline_1d__eval_array_deriv

  !-----------------------------------------------------------------------------
  !> @brief      Evaluate integral of 1D spline across the whole domain
  !> @param[in]  self  1D spline
  !> @param[out] y     value of integral
  !-----------------------------------------------------------------------------
  LOG_PURE function sll_spline_1d__integrate( self ) result( y )

    type(sll_spline_1d_t), intent(in   ) :: self
    real(c_double)  :: y

    y = spline_1d_integrate( self % spline_obj )

  end function sll_spline_1d__integrate

end module sll_m_spline_1d
