!> @ingroup splines
!> @brief   Module for tensor-product 2D splines.
!>
!> @details Here we define a 2D tensor-product spline type as an element of the
!>          linear space given by the span of 2D basis functions, which in turn
!>          are obtained as tensor-product of 1D B-splines.
!>          A 2D tensor-product B-splines type is not implemented, because an
!>          object of this type is completely described by the combination of
!>          two separate 1D B-splines objects, but this decision may be changed.
!>          Therefore, initialization of a 2D spline object requires two existing
!>          B-splines objects, to which two private (polymorphic) pointers are
!>          associated.
!>          The B-spline coefficients are stored in a 2D public allocatable array;
!>          at initialization the array is allocated to the proper shape and all
!>          values are set to zero.
!>          In most situations the B-spline coefficients are not set directly by
!>          the end user, but are computed by some other object (e.g., a Poisson
!>          solver or a spline interpolator).
!>          Various public methods allow the user to evaluate the 2D spline
!>          S(x1,x2) and its partial derivatives ∂S(x1,x2)/∂x1 and ∂S(x1,x2)/∂x2
!>          at any position (x1,x2).
!>
!> @author  Yaman Güçlü  - IPP Garching
!> @author  Edoardo Zoni - IPP Garching

module sll_m_spline_2d
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "logging.inc"

  use iso_C_binding, only: c_ptr, c_double, c_null_ptr, c_bool
  use sll_m_bsplines         , only: sll_bsplines_t

  implicit none

  public ::                               &
    sll_spline_2d_t                       , &
    sll_spline_2d__new                  , &
    sll_spline_2d__free                 , &
    sll_spline_2d__eval                 , &
    sll_spline_2d__eval_deriv_x1        , &
    sll_spline_2d__eval_deriv_x2        , &
    sll_spline_2d__eval_deriv_x1x2      , &
    sll_spline_2d__eval_array           , &
    sll_spline_2d__eval_array_deriv_x1  , &
    sll_spline_2d__eval_array_deriv_x2  , &
    sll_spline_2d__eval_array_deriv_x1x2, &
    sll_spline_2d__integrate

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


  !> 2D spline
  type :: sll_spline_2d_t
      type(c_ptr) :: spline_obj
  end type sll_spline_2d_t

  interface
      function spline_2d_new( &
              bspl1  , &
              bspl2    ) bind(C)
          use iso_C_binding, only: c_ptr
          type(c_ptr), value, intent(in) :: bspl1
          type(c_ptr), value, intent(in) :: bspl2
          type(c_ptr)            :: spline_2d_new
      end function
  end interface

  interface
      subroutine spline_2d_free( &
              spl ) bind(C)
          use iso_C_binding, only: c_ptr
          type(c_ptr), value, intent(in) :: spl
      end subroutine
  end interface

  interface
      LOG_PURE function spline_2d_eval( &
              spl, &
              x1 , &
              x2 ) bind(C)
          use iso_C_binding, only: c_ptr, c_double
          type(c_ptr), value   , intent(in) :: spl
          real(c_double)       , intent(in) :: x1
          real(c_double)       , intent(in) :: x2
          real(c_double)                    :: spline_2d_eval
      end function
  end interface

  interface
      LOG_PURE function spline_2d_eval_deriv( &
              spl, &
              x1       , &
              x2       , &
              deriv_x1 , &
              deriv_x2 ) bind(C)
          use iso_C_binding, only: c_ptr, c_double, c_bool
          type(c_ptr), value     , intent(in) :: spl
          real(c_double)         , intent(in) :: x1
          real(c_double)         , intent(in) :: x2
          logical(c_bool), value , intent(in) :: deriv_x1
          logical(c_bool), value , intent(in) :: deriv_x2
          real(c_double)                    :: spline_2d_eval_deriv
      end function
  end interface

  interface
      LOG_PURE subroutine spline_2d_eval_array( &
              spl   , &
              x1    , &
              nx1_1 , &
              nx1_2 , &
              x2    , &
              nx2_1 , &
              nx2_2 , &
              y     , &
              ny1   , &
              ny2      ) bind(C)
          use iso_C_binding, only: c_ptr, c_double, c_int
          type(c_ptr), value      , intent(in   ) :: spl
          integer(c_int), value   , intent(in   ) :: nx1_1
          integer(c_int), value   , intent(in   ) :: nx1_2
          real(c_double)          , intent(in   ) :: x1(0:nx1_1-1,0:nx1_2-1)
          integer(c_int), value   , intent(in   ) :: nx2_1
          integer(c_int), value   , intent(in   ) :: nx2_2
          real(c_double)          , intent(in   ) :: x2(0:nx2_1-1,0:nx2_2-1)
          integer(c_int), value   , intent(in   ) :: ny1
          integer(c_int), value   , intent(in   ) :: ny2
          real(c_double)          , intent(inout) :: y(0:ny1-1,0:ny2-1)
      end subroutine
  end interface

  interface
      LOG_PURE subroutine spline_2d_eval_array_deriv( &
              spl      , &
              x1       , &
              nx1_1    , &
              nx1_2    , &
              x2       , &
              nx2_1    , &
              nx2_2    , &
              deriv_x1 , &
              deriv_x2 , &
              y        , &
              ny1      , &
              ny2      ) bind(C)
          use iso_C_binding, only: c_ptr, c_double, c_int, c_bool
          type(c_ptr), value      , intent(in   ) :: spl
          integer(c_int), value   , intent(in   ) :: nx1_1
          integer(c_int), value   , intent(in   ) :: nx1_2
          real(c_double)          , intent(in   ) :: x1(0:nx1_1-1,0:nx1_2-1)
          integer(c_int), value   , intent(in   ) :: nx2_1
          integer(c_int), value   , intent(in   ) :: nx2_2
          real(c_double)          , intent(in   ) :: x2(0:nx2_1-1,0:nx2_2-1)
          logical(c_bool), value  , intent(in   ) :: deriv_x1
          logical(c_bool), value  , intent(in   ) :: deriv_x2
          integer(c_int), value   , intent(in   ) :: ny1
          integer(c_int), value   , intent(in   ) :: ny2
          real(c_double)          , intent(inout) :: y(0:ny1-1,0:ny2-1)
      end subroutine
  end interface

  !interface
  !    subroutine spline_2d_integrate_dim( &
  !            spl , &
  !            y   , &
  !            ny  , &
  !            int_dim ) bind(C)
  !        use iso_C_binding, only: c_ptr, c_double
  !        type(c_ptr), value   , intent(in) :: spl
  !        real(c_double)       , intent(in) :: y(0:ny-1)
  !        integer(c_int), value, intent(in) :: ny
  !        integer(c_int), value, intent(in) :: int_dim
  !        real(c_double)                    :: spline_2d_integrate
  !    end subroutine
  !end interface

  interface
      LOG_PURE function spline_2d_integrate( &
              spl ) bind(C)
          use iso_C_binding, only: c_ptr, c_double
          type(c_ptr), value   , intent(in) :: spl
          real(c_double)                    :: spline_2d_integrate
      end function
  end interface

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !-----------------------------------------------------------------------------
  !> @brief      Initialize 2D spline object as element of span(B-splines)
  !> @param[out] self         2D spline: new element of tensor-product space
  !> @param[in]  bsplines_x1  B-splines along x1
  !> @param[in]  bsplines_x2  B-splines along x2
  !-----------------------------------------------------------------------------
  function sll_spline_2d__new( bsplines_x1, bsplines_x2 ) result(self)

    type(sll_spline_2d_t)                        :: self
    type(sll_bsplines_t) , intent(in   ), target :: bsplines_x1
    type(sll_bsplines_t) , intent(in   ), target :: bsplines_x2

    self % spline_obj = spline_2d_new(bsplines_x1 % bspline_obj, &
                                      bsplines_x2  % bspline_obj)
  end function sll_spline_2d__new

  !-----------------------------------------------------------------------------
  !> @brief        Destroy 2D spline (re-initialization is possible afterwards)
  !> @param[inout] self  2D spline
  !-----------------------------------------------------------------------------
  subroutine sll_spline_2d__free( self )

    type(sll_spline_2d_t), intent(inout) :: self

    call spline_2d_free(self % spline_obj)
    self % spline_obj = c_null_ptr

  end subroutine sll_spline_2d__free

  !-----------------------------------------------------------------------------
  !> @brief     Evaluate value of 2D spline at location (x1,x2)
  !> @param[in] self  2D spline
  !> @param[in] x1    x1 coordinate of evaluation point
  !> @param[in] x2    x2 coordinate of evaluation point
  !> @returns         spline value y=S(x1,x2)
  !-----------------------------------------------------------------------------
  LOG_PURE function sll_spline_2d__eval( self, x1, x2 ) result( y )

    type(sll_spline_2d_t), intent(in) :: self
    real(c_double)              , intent(in) :: x1
    real(c_double)              , intent(in) :: x2
    real(c_double) :: y

    y = spline_2d_eval(self % spline_obj, x1, x2)

  end function sll_spline_2d__eval

  !-----------------------------------------------------------------------------
  !> @brief     Evaluate x1-derivative of 2D spline at location (x1,x2)
  !> @param[in] self  2D spline
  !> @param[in] x1    x1 coordinate of evaluation point
  !> @param[in] x2    x2 coordinate of evaluation point
  !> @returns         spline partial derivative y=∂S(x1,x2)/∂x1
  !-----------------------------------------------------------------------------
  LOG_PURE function sll_spline_2d__eval_deriv_x1( self, x1, x2 ) result( y )

    type(sll_spline_2d_t), intent(in) :: self
    real(c_double)              , intent(in) :: x1
    real(c_double)              , intent(in) :: x2
    real(c_double) :: y

    logical(c_bool) :: dx1
    logical(c_bool) :: dx2
    dx1 = .true.
    dx2 = .false.

    y = spline_2d_eval_deriv(self % spline_obj, x1, x2, dx1, dx2)

  end function sll_spline_2d__eval_deriv_x1

  !-----------------------------------------------------------------------------
  !> @brief     Evaluate x2-derivative of 2D spline at location (x1,x2)
  !> @param[in] self  2D spline
  !> @param[in] x1    x1 coordinate of evaluation point
  !> @param[in] x2    x2 coordinate of evaluation point
  !> @returns         spline partial derivative y=∂S(x1,x2)/∂x2
  !-----------------------------------------------------------------------------
  LOG_PURE function sll_spline_2d__eval_deriv_x2( self, x1, x2 ) result( y )

    type(sll_spline_2d_t), intent(in) :: self
    real(c_double)              , intent(in) :: x1
    real(c_double)              , intent(in) :: x2
    real(c_double) :: y

    logical(c_bool) :: dx1
    logical(c_bool) :: dx2
    dx1 = .false.
    dx2 = .true.

    y = spline_2d_eval_deriv(self % spline_obj, x1, x2, dx1, dx2)

  end function sll_spline_2d__eval_deriv_x2

  !-----------------------------------------------------------------------------
  !> @brief     Evaluate x1-x2 mixed derivative of 2D spline at location (x1,x2)
  !> @param[in] self  2D spline
  !> @param[in] x1    x1 coordinate of evaluation point
  !> @param[in] x2    x2 coordinate of evaluation point
  !> @returns         spline partial derivative y = ∂^2/(∂x1 ∂x2) S(x1,x2)
  !-----------------------------------------------------------------------------
  LOG_PURE function sll_spline_2d__eval_deriv_x1x2( self, x1, x2 ) result( y )

    type(sll_spline_2d_t), intent(in) :: self
    real(c_double)              , intent(in) :: x1
    real(c_double)              , intent(in) :: x2
    real(c_double) :: y

    logical(c_bool) :: dx1
    logical(c_bool) :: dx2
    dx1 = .true.
    dx2 = .true.

    y = spline_2d_eval_deriv(self % spline_obj, x1, x2, dx1, dx2)

  end function sll_spline_2d__eval_deriv_x1x2

  !-----------------------------------------------------------------------------
  !> @brief      Evaluate value of 2D spline at multiple locations
  !> @param[in]  self  2D spline
  !> @param[in]  x1    2D array of x1 coordinates of evaluation points
  !> @param[in]  x2    2D array of x2 coordinates of evaluation points
  !> @param[out] y     2D array of spline values y[i,j]=S(x1[i,j],x2[i,j])
  !-----------------------------------------------------------------------------
  LOG_PURE subroutine sll_spline_2d__eval_array( self, x1, x2, y )

    type(sll_spline_2d_t), intent(in   ) :: self
    real(c_double)              , intent(in   ) :: x1(:,:)
    real(c_double)              , intent(in   ) :: x2(:,:)
    real(c_double)              , intent(  out) :: y (:,:)

    call spline_2d_eval_array(self % spline_obj, x1, size(x1,1), size(x1,2), &
        x2, size(x2,1), size(x2,2), y, size(y,1), size(y,2))

  end subroutine sll_spline_2d__eval_array

  !-----------------------------------------------------------------------------
  !> @brief      Evaluate x1-derivative of 2D spline at multiple locations
  !> @param[in]  self  2D spline
  !> @param[in]  x1    2D array of x1 coordinates of evaluation points
  !> @param[in]  x2    2D array of x2 coordinates of evaluation points
  !> @param[out] y     2D array of spline partial derivatives
  !>                   y[i,j]=∂S(x1[i,j],x2[i,j])/∂x1
  !-----------------------------------------------------------------------------
  LOG_PURE subroutine sll_spline_2d__eval_array_deriv_x1( self, x1, x2, y )

    type(sll_spline_2d_t), intent(in   ) :: self
    real(c_double)              , intent(in   ) :: x1(:,:)
    real(c_double)              , intent(in   ) :: x2(:,:)
    real(c_double)              , intent(  out) :: y (:,:)

    logical(c_bool) :: dx1
    logical(c_bool) :: dx2
    dx1 = .true.
    dx2 = .false.

    call spline_2d_eval_array_deriv(self % spline_obj, x1, size(x1,1), size(x1,2), &
        x2, size(x2,1), size(x2,2), dx1, dx2, y, size(y,1), size(y,2))

  end subroutine sll_spline_2d__eval_array_deriv_x1

  !-----------------------------------------------------------------------------
  !> @brief      Evaluate x2-derivative of 2D spline at multiple locations
  !> @param[in]  self  2D spline
  !> @param[in]  x1    2D array of x1 coordinates of evaluation points
  !> @param[in]  x2    2D array of x2 coordinates of evaluation points
  !> @param[out] y     2D array of spline partial derivatives
  !>                   y[i,j]=∂S(x1[i,j],x2[i,j])/∂x2
  !-----------------------------------------------------------------------------
  LOG_PURE subroutine sll_spline_2d__eval_array_deriv_x2( self, x1, x2, y )

    type(sll_spline_2d_t), intent(in   ) :: self
    real(c_double)              , intent(in   ) :: x1(:,:)
    real(c_double)              , intent(in   ) :: x2(:,:)
    real(c_double)              , intent(  out) :: y (:,:)

    logical(c_bool) :: dx1
    logical(c_bool) :: dx2
    dx1 = .false.
    dx2 = .true.

    call spline_2d_eval_array_deriv(self % spline_obj, x1, size(x1,1), size(x1,2), &
        x2, size(x2,1), size(x2,2), dx1, dx2, y, size(y,1), size(y,2))

  end subroutine sll_spline_2d__eval_array_deriv_x2

  !-----------------------------------------------------------------------------
  !> @brief      Evaluate x1-x2 mixed derivative of 2D spline at multiple locs.
  !> @param[in]  self  2D spline
  !> @param[in]  x1    2D array of x1 coordinates of evaluation points
  !> @param[in]  x2    2D array of x2 coordinates of evaluation points
  !> @param[out] y     2D array of spline partial derivatives
  !>                   y[i,j] = ∂^2/(∂x1 ∂x2) S(x1[i,j],x2[i,j])
  !-----------------------------------------------------------------------------
  LOG_PURE subroutine sll_spline_2d__eval_array_deriv_x1x2( self, x1, x2, y )

    type(sll_spline_2d_t), intent(in   ) :: self
    real(c_double)              , intent(in   ) :: x1(:,:)
    real(c_double)              , intent(in   ) :: x2(:,:)
    real(c_double)              , intent(  out) :: y (:,:)

    logical(c_bool) :: dx1
    logical(c_bool) :: dx2
    dx1 = .true.
    dx2 = .true.

    call spline_2d_eval_array_deriv(self % spline_obj, x1, size(x1,1), size(x1,2), &
        x2, size(x2,1), size(x2,2), dx1, dx2, y, size(y,1), size(y,2))

  end subroutine sll_spline_2d__eval_array_deriv_x1x2

  !-----------------------------------------------------------------------------
  !> @brief      Integrate the spline over the chosen dimension
  !> @param[in]  self  2D spline
  !> @param[in]  y      Integral of spline
  !>                    y = ∫∫ S(x1,x2) dx1 dx2 
  !-----------------------------------------------------------------------------
  LOG_PURE function sll_spline_2d__integrate( self ) result( y )

    type(sll_spline_2d_t), intent(in   ) :: self
    real(c_double)                     :: y

    y = spline_2d_integrate(self % spline_obj)

  end function sll_spline_2d__integrate

end module sll_m_spline_2d
