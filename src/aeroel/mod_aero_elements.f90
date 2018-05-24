!!=====================================================================
!!
!! Copyright (C) 2018 Politecnico di Milano
!!
!! This file is part of DUST, an aerodynamic solver for complex
!! configurations.
!! 
!! Permission is hereby granted, free of charge, to any person
!! obtaining a copy of this software and associated documentation
!! files (the "Software"), to deal in the Software without
!! restriction, including without limitation the rights to use,
!! copy, modify, merge, publish, distribute, sublicense, and/or sell
!! copies of the Software, and to permit persons to whom the
!! Software is furnished to do so, subject to the following
!! conditions:
!! 
!! The above copyright notice and this permission notice shall be
!! included in all copies or substantial portions of the Software.
!! 
!! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
!! EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
!! OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
!! NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
!! HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
!! WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
!! FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
!! OTHER DEALINGS IN THE SOFTWARE.
!! 
!! Authors: 
!!          Federico Fonte             <federico.fonte@polimi.it>
!!          Davide Montagnani       <davide.montagnani@polimi.it>
!!          Matteo Tugnoli             <matteo.tugnoli@polimi.it>
!!=====================================================================

!> Module containing the upper level definitions of the aerodynamic
!! elements


module mod_aero_elements

use mod_linsys_vars, only: t_linsys

use mod_param, only: &
  wp

implicit none

public :: c_elem, t_elem_p

private

!----------------------------------------------------------------------

type :: t_elem_p
  class(c_elem), pointer :: p
end type

!> Abstract type defining a generic aerodynamic element
!!
!! It needs to be extended into more specific elements
type, abstract :: c_elem

  !> Intensity of the doublets/vortexes
  real(wp), pointer :: idou => null()
  real(wp)          :: didou_dt
 
  !> Element id
  integer :: id

  !> id of the component to which it belongs
  integer :: comp_id

  !> Number of vertexes
  integer :: n_ver  
  !> Vertexes coordinates
  !real(wp), pointer :: ver(:,:)
  real(wp), allocatable :: ver(:,:)

  !> Global index of vertexes
  integer,  allocatable :: i_ver(:)
  !> Element area
  real(wp)              :: area
  !> Element centre coordinates
  real(wp), allocatable :: cen(:)
  !> Element normal vector
  real(wp), allocatable :: nor(:)
  !> Element tangential vectors
  real(wp), allocatable :: tang(:,:) ! tangent unit vectors as in PANAIR
  real(wp), allocatable :: cosTi(:) , sinTi(:)
  real(wp), allocatable :: verp(:,:)
  real(wp), allocatable :: edge_vec(:,:)
  real(wp), allocatable :: edge_len(:)
  real(wp), allocatable :: edge_uni(:,:)
  !> Body velocity at the centre
  real(wp), allocatable :: ub(:)
  !> Is the element moving during simulation?
  logical :: moving
  !> Element neighbours global index (in the elements vector)
  !integer, allocatable :: i_neigh(:)
  type(t_elem_p), allocatable :: neigh(:)

  !> Coefficients to compute local velocity from the velocity potential
  !! on a stencil of neighboring elements (for surfpan only)
  real(wp), allocatable :: pot_vel_stencil(:,:)

  !> Panel width (= strip width) (for vortring only)
  real(wp)              :: dy
  !> Previous element in a stripe (for vortring only)
  !integer               :: stripe_1
  type(t_elem_p)        :: stripe_1
  !> Element indices in the component%strip_elem array (for vortring only)
  !integer, allocatable  :: stripe_elem(:)
  type(t_elem_p), allocatable :: stripe_elem(:)

  !> Fluid velocity at center for boundary condition (U_inf-rel vel)
  real(wp), allocatable :: vel(:)
  !> Fluid velocity at center for boundary condition (U_inf-rel vel)
  real(wp)              :: cp

  contains
! procedure(i_velocity_calc), deferred, pass(this) :: velocity_calc
! procedure(i_potential_calc), deferred, pass(this) :: potential_calc

! linear system ------
! new
  procedure(i_build_row)  , deferred, pass(this)      :: build_row
  procedure(i_build_row_static), deferred, pass(this) :: build_row_static
  procedure(i_add_wake),    deferred, pass(this)      :: add_wake
  procedure(i_add_liftlin), deferred, pass(this)      :: add_liftlin
  procedure(i_add_actdisk), deferred, pass(this)      :: add_actdisk
  procedure(i_compute_pot), deferred, pass(this)      :: compute_pot
  procedure(i_compute_vel), deferred, pass(this)      :: compute_vel
  procedure(i_compute_psi), deferred, pass(this)      :: compute_psi
! loads --------------
  procedure(i_compute_cp ), deferred, pass(this)      :: compute_cp
  

end type


!----------------------------------------------------------------------

abstract interface
  subroutine i_velocity_calc(this, vel, pos)
    import :: c_elem , wp
    implicit none
    class(c_elem), intent(inout) :: this
    real(wp), intent(out) :: vel
    real(wp), intent(in) :: pos(:)
  end subroutine
end interface

!----------------------------------------------------------------------

abstract interface
  subroutine i_potential_calc(this, pot, pos)
    import :: c_elem , wp
    implicit none
    class(c_elem), intent(inout) :: this
    real(wp), intent(out) :: pot
    real(wp), intent(in) :: pos(:)
  end subroutine
end interface

!----------------------------------------------------------------------

!> Build a complete row of the linear system, including the right hand
!! side
!!
!! Since the static part of the sistem should already be initialized, 
!! the function updates just the coefficients from ista to iend, 
!! which should be the moving elements ordered at the end of the 
!! elements array
!!
!! The static part od the rhs should already be initialized, and the 
!! moving contribution is just added to the static one.
abstract interface
  subroutine i_build_row(this, elems, linsys, uinf, ie, ista, iend)
    import :: wp
    import :: c_elem
    import :: t_elem_p
    import :: t_linsys
    implicit none
    class(c_elem), intent(inout)  :: this
    type(t_elem_p), intent(in)    :: elems(:)
    type(t_linsys), intent(inout) :: linsys
    real(wp), intent(in)          :: uinf(:)
    integer, intent(in)           :: ie
    integer, intent(in)           :: ista, iend
  end subroutine
end interface

!----------------------------------------------------------------------

!> Build the static part of a row of the linear system
!!
!! It is necessary to calculate the aerodynamic influence coefficients as
!! well as the rhs contribution of the static elements just once at the 
!! beginning of the simulation. 
!!
!! For this reason this subroutine is called to calculate the static aic 
!! and the static contribution to the rhs. The static contribution to the 
!! rhs is stored into a 3 components * number of elements array so that 
!! at each timestep the actual rhs contribution can be calculated
!! multiplying the 3 components for the velocity vector, which in principle
!! could vary.
!!
!! The function operates only on the components from ista to iend, which 
!! should be the static ones ordered at the beginning of the array
abstract interface
  subroutine i_build_row_static(this, elems, ll_elems, ad_elems, linsys, uinf, ie, ista, iend)
    import :: wp
    import :: c_elem
    import :: t_elem_p
    import :: t_linsys
    implicit none
    class(c_elem), intent(inout)  :: this
    type(t_elem_p), intent(in)    :: elems(:)
    type(t_elem_p), intent(in)    :: ll_elems(:)
    type(t_elem_p), intent(in)    :: ad_elems(:)
    type(t_linsys), intent(inout) :: linsys
    real(wp), intent(in)          :: uinf(:)
    integer, intent(in)           :: ie
    integer, intent(in)           :: ista, iend
  end subroutine
end interface

!----------------------------------------------------------------------

abstract interface
  subroutine i_add_wake(this, wake_elems, impl_wake_ind, linsys, uinf, &
                        ie, ista, iend)
    import :: wp
    import :: c_elem
    import :: t_elem_p
    import :: t_linsys
    implicit none
    class(c_elem), intent(inout)  :: this
    type(t_elem_p), intent(in)    :: wake_elems(:)
    integer, intent(in)           :: impl_wake_ind(:,:)
    type(t_linsys), intent(inout) :: linsys
    real(wp), intent(in)          :: uinf(:)
    integer, intent(in)           :: ie
    integer, intent(in)           :: ista
    integer, intent(in)           :: iend
  end subroutine
end interface

!----------------------------------------------------------------------

abstract interface
  subroutine i_add_liftlin(this, ll_elems, linsys, uinf, &
                        ie, ista, iend)
    import :: wp
    import :: c_elem
    import :: t_elem_p
    import :: t_linsys
    implicit none
    class(c_elem), intent(inout)  :: this
    type(t_elem_p), intent(in)    :: ll_elems(:)
    type(t_linsys), intent(inout) :: linsys
    real(wp), intent(in)          :: uinf(:)
    integer, intent(in)           :: ie
    integer, intent(in)           :: ista
    integer, intent(in)           :: iend
  end subroutine
end interface

!----------------------------------------------------------------------

abstract interface
  subroutine i_add_actdisk(this, ad_elems, linsys, uinf, &
                        ie, ista, iend)
    import :: wp
    import :: c_elem
    import :: t_elem_p
    import :: t_linsys
    implicit none
    class(c_elem), intent(inout)  :: this
    type(t_elem_p), intent(in)    :: ad_elems(:)
    type(t_linsys), intent(inout) :: linsys
    real(wp), intent(in)          :: uinf(:)
    integer, intent(in)           :: ie
    integer, intent(in)           :: ista
    integer, intent(in)           :: iend
  end subroutine
end interface

!----------------------------------------------------------------------

abstract interface
  subroutine i_compute_pot(this, A, b, pos,i,j)
    import :: c_elem , wp , t_linsys
    implicit none
    class(c_elem), intent(inout) :: this
    real(wp), intent(out) :: A
    real(wp), intent(out) :: b(3)
    real(wp), intent(in) :: pos(:)
    integer , intent(in) :: i,j
  end subroutine
end interface

!----------------------------------------------------------------------

abstract interface
  subroutine i_compute_vel(this, pos, uinf, vel)
    import :: c_elem , wp 
    implicit none
    class(c_elem), intent(inout) :: this
    real(wp), intent(in) :: pos(:)
    real(wp), intent(in) :: uinf(3)
    real(wp), intent(out) :: vel(3)
  end subroutine
end interface

!----------------------------------------------------------------------

abstract interface
  subroutine i_compute_psi(this, A, b, pos, nor,i,j)
    import :: c_elem , wp , t_linsys
    implicit none
    class(c_elem), intent(inout) :: this
    real(wp), intent(out) :: A
    real(wp), intent(out) :: b(3)
    real(wp), intent(in) :: pos(:)
    real(wp), intent(in) :: nor(:)
    integer , intent(in) :: i,j
  end subroutine
end interface

!----------------------------------------------------------------------

abstract interface
  subroutine i_compute_cp (this, elems, uinf)
    import :: c_elem , t_elem_p , wp 
    implicit none
    class(c_elem), intent(inout) :: this
    type(t_elem_p), intent(in)   :: elems(:)
    real(wp), intent(in)         :: uinf(:)
  end subroutine
end interface

!----------------------------------------------------------------------

abstract interface
  subroutine i_build_aic(this, elems, linsys, ie)
    import :: c_elem
    import :: t_elem_p
    import :: t_linsys
    implicit none
    class(c_elem), intent(inout)  :: this
    type(t_elem_p), intent(in)    :: elems(:)
    type(t_linsys), intent(inout) :: linsys
    integer                       :: ie
  end subroutine
end interface

!----------------------------------------------------------------------

abstract interface
  subroutine i_build_aic_dyn(this, elems, linsys, ie)
    import :: c_elem
    import :: t_elem_p
    import :: t_linsys
    implicit none
    class(c_elem), intent(inout)  :: this
    type(t_elem_p), intent(in)    :: elems(:)
    type(t_linsys), intent(inout) :: linsys
    integer                       :: ie
  end subroutine
end interface

!----------------------------------------------------------------------

abstract interface
  subroutine i_build_rhs(this, elems, linsys, ie)
    import :: c_elem
    import :: t_elem_p
    import :: t_linsys
    implicit none
    class(c_elem), intent(inout)  :: this
    type(t_elem_p), intent(in)    :: elems(:)
    type(t_linsys), intent(inout) :: linsys
    integer                       :: ie
  end subroutine
end interface

!----------------------------------------------------------------------


end module mod_aero_elements
