!> @ingroup boundary_condition_descriptors
!> @brief Describe different boundary conditions
!> @details
!> The intent of this module is to provide a single, library-wide definition
!> of the names used to describe different boundary conditions. One should 
!> ALWAYS refer to specific boundary conditions by their
!> names and not through their integer representation, which could be changed.
!>
!> <b> How to use-it </b>
!>
!> Just add the line
!> @code
!> #include "sll_m_boundary_condition_descriptors.h"
!> @endcode
!
! To be considered here is to include also BC combinations, which may help
! save some coding instead of managing this internally within routines, for
! example a flag like SLL_DIRICHLET_NEUMANN could indicate two BC's along
! a particular dimension...
module sll_m_boundary_condition_descriptors
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  implicit none

  public :: &
    sll_p_conductor, &
    sll_p_dirichlet, &
    sll_p_halo, &
    sll_p_one_sided, &
    sll_p_hermite, &
    sll_p_interior, &
    sll_p_neumann, &
    sll_p_neumann_mode_0, &
    sll_p_polar_origin, &
    sll_p_periodic, &
    sll_p_open, &
    sll_p_mirror, &
    sll_p_greville, &
    sll_p_set_to_limit, &
    sll_p_silver_muller, &
    sll_p_user_defined

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !> User defined boundary condition
  integer, parameter :: sll_p_user_defined   = -1 
  !> Periodic boundary condition u(1)=u(n)
  integer, parameter :: sll_p_periodic       = 0 
  !> Dirichlet boundary condition 
  integer, parameter :: sll_p_dirichlet      = 1 
  !> Neumann boundary condition 
  integer, parameter :: sll_p_neumann        = 2
  !> Hermite boundary condition
  integer, parameter :: sll_p_hermite        = 3
  !> Neumann boundary condition
  integer, parameter :: sll_p_neumann_mode_0 = 4
  !> PLEASE ADD DOCUMENTATION
  integer, parameter :: sll_p_set_to_limit   = 5
  !> Interior of domain
  integer, parameter :: sll_p_interior       = 6
  !> Incoming wave boundar condition for Maxwell
  integer, parameter :: SLL_INCOMING_WAVE  = 7
  !> Metallic boundary condition for Maxwell
  integer, parameter :: sll_p_conductor      = 8
  !> Absorbing boundary condition fro Maxwell
  integer, parameter :: sll_p_silver_muller  = 9
  !> Use a one-sided stencil at the boundary
  integer, parameter :: sll_p_one_sided    = 10
  !> Values outside the domain are provided as halo cells (for domain decomposition)
  integer, parameter :: sll_p_halo         = 11
  !> Use Greville points instead of conditions on derivative for B-Spline interpolation
  integer, parameter :: sll_p_greville      = 12
  !> Duplicate boundary points to define knots for B-splines from grid
  integer, parameter :: sll_p_open         = 13
  !> Mirror points around boundary to define knots for B-splines from grid
  integer, parameter :: sll_p_mirror       = 14
  !> Treat boundary point as origin of polar coordinates (= center of a circle)
  integer, parameter :: sll_p_polar_origin = 15


end module sll_m_boundary_condition_descriptors
