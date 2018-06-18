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


!> Module containing the specific subroutines for the vortex lattice
!! type of aerodynamic elements
module mod_vortlatt

!use mod_aero_elements, only: &
!  c_elem, t_elem_p

use mod_doublet, only: &
  potential_calc_doublet , &
  velocity_calc_doublet

use mod_linsys_vars, only: &
  t_linsys

use mod_sim_param, only: &
  t_sim_param

use mod_param, only: &
  wp, pi

use mod_aeroel, only: &
  c_elem, c_pot_elem, c_vort_elem, c_impl_elem, c_expl_elem, &
  t_elem_p, t_pot_elem_p, t_vort_elem_p, t_impl_elem_p, t_expl_elem_p
!----------------------------------------------------------------------

implicit none

public :: t_vortlatt


!----------------------------------------------------------------------

!> Planar aerodynamic elements with a surface distribution of doublets
!!
!! Aerodynamic elements characterized by a uniform surface distribution
!! of doublets, employed to model a planar surface.
!! Thanks to the equivalence of a surface distribution of doublets and a
!! closed vortex ring it can be also viewed as ring vortex laying on the
!! edges of the element.
type, extends(c_impl_elem) :: t_vortlatt

  !TODO: consider applying the correct element pointer here
contains

procedure, pass(this) :: build_row        => build_row_vortlatt
procedure, pass(this) :: build_row_static => build_row_static_vortlatt
procedure, pass(this) :: add_wake         => add_wake_vortlatt
procedure, pass(this) :: add_expl         => add_expl_vortlatt
procedure, pass(this) :: compute_pot      => compute_pot_vortlatt
procedure, pass(this) :: compute_vel      => compute_vel_vortlatt
procedure, pass(this) :: compute_psi      => compute_psi_vortlatt
procedure, pass(this) :: compute_pres     => compute_pres_vortlatt
procedure, pass(this) :: compute_dforce   => compute_dforce_vortlatt
end type


character(len=*), parameter :: this_mod_name='mod_vortlatt'

!----------------------------------------------------------------------
contains
!----------------------------------------------------------------------

!> Build a row of the linear system for a vortex ring
!!
!! Only the dynamic part of the linear system is actually built here:
!! the rest of the system was already built in the \ref build_row_static
!! subroutine.
subroutine build_row_vortlatt(this, elems, linsys, uinf, ie, ista, iend)
class(t_vortlatt), intent(inout) :: this
type(t_impl_elem_p), intent(in)       :: elems(:)
type(t_linsys), intent(inout)    :: linsys
real(wp), intent(in)             :: uinf(:)
integer, intent(in)              :: ie
integer, intent(in)              :: ista, iend

integer :: j1
real(wp) :: b1(3)
!Not moving components, in the rhs contribution there is no body velocity
linsys%b(ie) = sum(linsys%b_static(:,ie) * (-uinf))

! ista and iend will be the end of the unknowns vector, containing
! the moving elements
do j1 = ista , iend

  call elems(j1)%p%compute_psi( linsys%A(ie,j1), b1,  &
                                this%cen, this%nor, ie, j1 )

  if (ie .eq. j1) then
    !diagonal, we are certainly employing vortrin, enforce the b.c. on ie
    linsys%b(ie) = linsys%b(ie) + sum(b1*(elems(ie)%p%ub-uinf))
  else
    ! off-diagonal: if it is a vortrin b1 is zero, if it is a surfpan
    ! enforce the boundary condition on it (j1)
    linsys%b(ie) = linsys%b(ie) + sum(b1*(elems(j1)%p%ub-uinf))
  endif

end do

end subroutine build_row_vortlatt

!----------------------------------------------------------------------

!> Build a static row of the linear system for a vortex ring
!!
!! In this subroutine only the static part of the equations is built. It is
!! called just once at the beginning of the simulation, and saves the AIC
!! coefficients for te static part and the static contribution to the rhs
subroutine build_row_static_vortlatt(this, elems, expl_elems, linsys, &
                                   uinf, ie, ista, iend)
class(t_vortlatt), intent(inout) :: this
type(t_impl_elem_p), intent(in)       :: elems(:)
type(t_expl_elem_p), intent(in)       :: expl_elems(:)
type(t_linsys), intent(inout)    :: linsys
real(wp), intent(in)             :: uinf(:)
integer, intent(in)              :: ie
integer, intent(in)              :: ista, iend

integer :: j1
real(wp) :: b1(3)

linsys%b(ie) = 0.0_wp
linsys%b_static(:,ie) = 0.0_wp

!Cycle just all the static elements, ista and iend will be the beginning of
!the result vector. Then save the rhs in b_static
do j1 = ista , iend

  call elems(j1)%p%compute_psi( linsys%A(ie,j1), b1,  &
                                this%cen, this%nor, ie, j1 )

  linsys%b_static(:,ie) = linsys%b_static(:,ie) + b1

end do

!!Now build the static contribution from the lifting line elements
!do j1 = 1,linsys%nstatic_ll
!  call ll_elems(j1)%p%compute_psi( linsys%L_static(ie,j1), b1,  &
!                                this%cen, this%nor,  1, 2 )
!enddo
!
!!Now build the static contribution from the lifting line elements
!do j1 = 1,linsys%nstatic_ad
!  call ad_elems(j1)%p%compute_psi( linsys%D_static(ie,j1), b1,  &
!                                this%cen, this%nor,  1, 2 )
!enddo

!Now build the static contribution from the lifting line elements
do j1 = 1,linsys%nstatic_expl
  call expl_elems(j1)%p%compute_psi( linsys%L_static(ie,j1), b1,  &
                                this%cen, this%nor,  1, 2 )
enddo
!The rest of the dynamic part will be completed during the first
! iteration of the assempling

end subroutine build_row_static_vortlatt

!----------------------------------------------------------------------

!> Add the contribution of the lifting lines to one equation for a vortex ring
!!
!! The rhs of the equation for a vortex ring is updated  adding the
!! the contribution of velocity due to the lifting lines
subroutine add_expl_vortlatt(this, expl_elems, linsys, uinf, &
                           ie, ista, iend)
class(t_vortlatt), intent(inout) :: this
type(t_expl_elem_p), intent(in)       :: expl_elems(:)
type(t_linsys), intent(inout)    :: linsys
real(wp), intent(in)             :: uinf(:)
integer, intent(in)              :: ie
integer, intent(in)             :: ista
integer, intent(in)             :: iend

integer :: j1
real(wp) :: a, b(3)


!Static part: take what was already computed
do  j1 = 1, ista-1
  linsys%b(ie) = linsys%b(ie) - linsys%L_static(ie,j1)*expl_elems(j1)%p%mag
enddo

! Add the explicit vortex panel wake contribution to the rhs
do j1 = ista, iend

  call expl_elems(j1)%p%compute_psi( a, b, this%cen, this%nor, 1, 2 )


  linsys%b(ie) = linsys%b(ie) - a*expl_elems(j1)%p%mag

end do

end subroutine add_expl_vortlatt

!----------------------------------------------------------------------

!> Add the contribution of the wake to one equation for a vortex ring
!!
!! The rhs of the equation for a surface panel is updated  adding the
!! the contribution of velocity due to the wake
subroutine add_wake_vortlatt(this, wake_elems, impl_wake_ind, linsys, uinf, &
                           ie, ista, iend)
class(t_vortlatt), intent(inout) :: this
type(t_pot_elem_p), intent(in)       :: wake_elems(:)
integer, intent(in)             :: impl_wake_ind(:,:)
type(t_linsys), intent(inout)    :: linsys
real(wp), intent(in)             :: uinf(:)
integer, intent(in)              :: ie
integer, intent(in)             :: ista
integer, intent(in)             :: iend

integer :: j1, ind1, ind2
real(wp) :: a, b(3)
integer :: n_impl

!Count the number of implicit wake contributions
n_impl = size(impl_wake_ind,2)

!Add the contribution of the implicit wake panels to the linear system
!Implicitly we assume that the first set of wake panels are the implicit
!ones since are at the beginning of the list
do j1 = 1 , n_impl
  ind1 = impl_wake_ind(1,j1); ind2 = impl_wake_ind(2,j1)
!   if ((ind1.ge.ista .and. ind1.le.iend) .and. &
!       (ind2.ge.ista .and. ind2.le.iend)) then
  if ((ind1.ge.ista .and. ind1.le.iend)) then

    call wake_elems(j1)%p%compute_psi( a, b, this%cen, this%nor, 1, 2 )

    linsys%A(ie,ind1) = linsys%A(ie,ind1) + a
    if ( ind2 .ne. 0 ) linsys%A(ie,ind2) = linsys%A(ie,ind2) - a

  endif
end do

! Add the explicit vortex panel wake contribution to the rhs
do j1 = n_impl+1 , size(wake_elems)

  call wake_elems(j1)%p%compute_psi( a, b, this%cen, this%nor, 1, 2 )


  linsys%b(ie) = linsys%b(ie) - a*wake_elems(j1)%p%mag

end do

end subroutine add_wake_vortlatt

!----------------------------------------------------------------------

!> Compute the potential due to a vortex ring
!!
!! this subroutine employs doublets basic subroutines to calculate
!! the AIC of a vortex ring on a surface panel. The contribution to its rhs
!! is zero since there are no sources (and no b.c. enforcing)
subroutine compute_pot_vortlatt(this, A, b, pos,i,j)
class(t_vortlatt), intent(inout) :: this
real(wp), intent(out) :: A
real(wp), intent(out) :: b(3)
real(wp), intent(in) :: pos(:)
integer , intent(in) :: i,j

real(wp) :: dou

if ( i .ne. j ) then
  call potential_calc_doublet(this, dou, pos)
else
!   AIC (doublets) = 0.0   -> dou = 0
  dou = -2.0_wp*pi
end if

A = -dou

b=0.0_wp

end subroutine compute_pot_vortlatt

!----------------------------------------------------------------------

!> Compute the velocity due to a vortex ring
!!
!! This subroutine employs basic doublet subroutines to calculate the AIC of
!! a vortex ring on a vortex ring. The contribution to the rhs is the boundary
!! condition and it is enforced only if the influencing and influenced
!! vortex rings are the same
subroutine compute_psi_vortlatt(this, A, b, pos, nor, i, j )
class(t_vortlatt), intent(inout) :: this
real(wp), intent(out) :: A
real(wp), intent(out) :: b(3)
real(wp), intent(in) :: pos(:)
real(wp), intent(in) :: nor(:)
integer , intent(in) :: i , j

real(wp) :: vdou(3)

call velocity_calc_doublet(this, vdou, pos)

A = sum(vdou * nor)


! b = ... (from boundary conditions)
!TODO: consider moving this outside
if ( i .eq. j ) then
  b =  4.0_wp*pi*this%nor
else
  b = 0.0_wp
end if

end subroutine compute_psi_vortlatt

!----------------------------------------------------------------------

!> Compute the velocity induced by a vortex ring in a prescribed position
!!
!! The velocity in the position is calculated considering the influece of
!! doublets
!!
!! WARNING: the velocity calculated, to be consistent with the formulation of
!! the equations is multiplied by 4*pi, to obtain the actual velocity the
!! result of the present subroutine MUST be DIVIDED by 4*pi
subroutine compute_vel_vortlatt(this, pos, uinf, vel)
class(t_vortlatt), intent(inout) :: this
real(wp), intent(in) :: pos(:)
real(wp), intent(in) :: uinf(3)
real(wp), intent(out) :: vel(3)

real(wp) :: vdou(3)


! doublet ---
call velocity_calc_doublet(this, vdou, pos)

vel = vdou*this%mag

end subroutine compute_vel_vortlatt

!----------------------------------------------------------------------

!> Compute an approximate value of the mean DELTA pressure on the actual
!! element
!!
!! pres = " DELTA pressure = ( pressure lower - pressure upper ) "
!!  s.t. vec{dforce} = pres * vec{n}  ( since vec{n} = vec{n_upper} )
!!
!! see compute_dforce_vortlatt
subroutine compute_pres_vortlatt(this, elems, sim_param)
class(t_vortlatt), intent(inout) :: this
type(t_elem_p), intent(in) :: elems(:)
type(t_sim_param), intent(in) :: sim_param

integer  :: i_stripe
real(wp) :: dG_dt

this%pres = 0.0_wp

i_stripe = size(this%stripe_elem)

dG_dt = this%didou_dt

if ( i_stripe .gt. 1 ) then
  !this%pres = - sim_param%rho_inf * &
  !      ( norm2(sim_param%u_inf - this%ub) * this%dy / this%area * &
  !      ( elems(this%id)%p%mag - this%stripe_elem(i_stripe-1)%p%mag ) + &
  !           dG_dt )
  ! dunno why elems(this%id)%p%mag was used
  this%pres = - sim_param%rho_inf * &
        ( norm2(sim_param%u_inf - this%ub) * this%dy / this%area * &
        ( this%mag - this%stripe_elem(i_stripe-1)%p%mag ) + &
             dG_dt )
else
  !this%pres = - sim_param%rho_inf * &
  !      ( norm2(sim_param%u_inf - this%ub) * this%dy / this%area * &
  !             elems(this%id)%p%mag + &
  !           dG_dt )
  ! the same...
  this%pres = - sim_param%rho_inf * &
        ( norm2(sim_param%u_inf - this%ub) * this%dy / this%area * &
               this%mag + &
             dG_dt )
end if


end subroutine compute_pres_vortlatt

!----------------------------------------------------------------------

!>  Compute the elementary force on the on the actual element
subroutine compute_dforce_vortlatt(this, elems, sim_param)
class(t_vortlatt), intent(inout) :: this
type(t_elem_p), intent(in) :: elems(:)
type(t_sim_param), intent(in) :: sim_param

this%dforce = this%pres * this%area * this%nor


end subroutine compute_dforce_vortlatt

!----------------------------------------------------------------------

end module mod_vortlatt
