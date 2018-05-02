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


!> Module containing the specific subroutines for the vortex ring
!! type of aerodynamic elements
module mod_vortring

use mod_aero_elements, only: &
  c_elem, t_elem_p

use mod_doublet, only: &
  potential_calc_doublet , &
  velocity_calc_doublet

use mod_linsys_vars, only: &
  t_linsys

use mod_param, only: &
  wp, pi

!----------------------------------------------------------------------

implicit none

public :: t_vortring


!----------------------------------------------------------------------

type, extends(c_elem) :: t_vortring

  real, allocatable :: vert(:,:)
  real, allocatable :: bar(:)
contains

! linear system ------
! new
  procedure, pass(this) :: build_row        => build_row_vortring
  procedure, pass(this) :: build_row_static => build_row_static_vortring
  procedure, pass(this) :: add_wake         => add_wake_vortring
  procedure, pass(this) :: add_liftlin      => add_liftlin_vortring
  procedure, pass(this) :: compute_pot      => compute_pot_vortring
  procedure, pass(this) :: compute_vel      => compute_vel_vortring
  procedure, pass(this) :: compute_psi      => compute_psi_vortring
  procedure, pass(this) :: compute_cp       => compute_cp_vortring
end type


character(len=*), parameter :: this_mod_name='mod_vortring'

!----------------------------------------------------------------------
contains
!----------------------------------------------------------------------

!> Build a row of the linear system for a vortex ring
!!
!! Only the dynamic part of the linear system is actually built here:
!! the rest of the system was already built in the \ref build_row_static 
!! subroutine. 
subroutine build_row_vortring(this, elems, linsys, uinf, ie, ista, iend)
 class(t_vortring), intent(inout) :: this
 type(t_elem_p), intent(in)       :: elems(:)
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

end subroutine build_row_vortring

!----------------------------------------------------------------------

!> Build a static row of the linear system for a vortex ring
!!
!! In this subroutine only the static part of the equations is built. It is
!! called just once at the beginning of the simulation, and saves the AIC 
!! coefficients for te static part and the static contribution to the rhs
subroutine build_row_static_vortring(this, elems, ll_elems, linsys, uinf, ie, ista, iend)
 class(t_vortring), intent(inout) :: this
 type(t_elem_p), intent(in)       :: elems(:)
 type(t_elem_p), intent(in)       :: ll_elems(:)
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

  !Now build the static contribution from the lifting line elements
  do j1 = 1,linsys%nstatic_ll
    call ll_elems(j1)%p%compute_psi( linsys%L_static(ie,j1), b1,  &
                                  this%cen, this%nor,  1, 2 )
  enddo
  
  !The rest of the dynamic part will be completed during the first 
  ! iteration of the assempling

end subroutine build_row_static_vortring

!----------------------------------------------------------------------

!> Add the contribution of the lifting lines to one equation for a vortex ring
!!
!! The rhs of the equation for a vortex ring is updated  adding the 
!! the contribution of velocity due to the lifting lines
subroutine add_liftlin_vortring(this, ll_elems, linsys, uinf, &
                             ie, ista, iend)
 class(t_vortring), intent(inout) :: this
 type(t_elem_p), intent(in)       :: ll_elems(:)
 type(t_linsys), intent(inout)    :: linsys
 real(wp), intent(in)             :: uinf(:)
 integer, intent(in)              :: ie
 integer, intent(in)             :: ista
 integer, intent(in)             :: iend

 integer :: j1, ind1, ind2
 real(wp) :: a, b(3)
 integer :: n_impl
  

  !Static part: take what was already computed
  do  j1 = 1, ista-1
    linsys%b(ie) = linsys%b(ie) - linsys%L_static(ie,j1)*ll_elems(j1)%p%idou
  enddo

  ! Add the explicit vortex panel wake contribution to the rhs
  do j1 = ista, iend

    call ll_elems(j1)%p%compute_psi( a, b, this%cen, this%nor, 1, 2 )


    linsys%b(ie) = linsys%b(ie) - a*ll_elems(j1)%p%idou

  end do

end subroutine add_liftlin_vortring

!----------------------------------------------------------------------

!> Add the contribution of the wake to one equation for a vortex ring
!!
!! The rhs of the equation for a surface panel is updated  adding the 
!! the contribution of velocity due to the wake
subroutine add_wake_vortring(this, wake_elems, impl_wake_ind, linsys, uinf, &
                             ie, ista, iend)
 class(t_vortring), intent(inout) :: this
 type(t_elem_p), intent(in)       :: wake_elems(:)
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


    linsys%b(ie) = linsys%b(ie) - a*wake_elems(j1)%p%idou

  end do

end subroutine add_wake_vortring

!----------------------------------------------------------------------

!> Compute the potential due to a vortex ring
!!
!! this subroutine employs doublets basic subroutines to calculate
!! the AIC of a vortex ring on a surface panel. The contribution to its rhs
!! is zero since there are no sources (and no b.c. enforcing) 
subroutine compute_pot_vortring(this, A, b, pos,i,j)
  class(t_vortring), intent(inout) :: this
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

  ! TODO: check coefficients 1/4*pi, ...
  A = -dou

  b=0.0_wp


end subroutine compute_pot_vortring

!----------------------------------------------------------------------

!> Compute the velocity due to a vortex ring
!!
!! This subroutine employs basic doublet subroutines to calculate the AIC of 
!! a vortex ring on a vortex ring. The contribution to the rhs is the boundary
!! condition and it is enforced only if the influencing and influenced
!! vortex rings are the same
subroutine compute_psi_vortring(this, A, b, pos, nor, i, j )
  class(t_vortring), intent(inout) :: this
  real(wp), intent(out) :: A
  real(wp), intent(out) :: b(3)
  real(wp), intent(in) :: pos(:)
  real(wp), intent(in) :: nor(:)
  integer , intent(in) :: i , j

  real(wp) :: vdou(3) 

  call velocity_calc_doublet(this, vdou, pos)

  ! TODO: check coefficients 1/4*pi, ...
  A = sum(vdou * nor)


! b = ... (from boundary conditions)
!TODO: consider moving this outside
  if ( i .eq. j ) then
    b =  4.0_wp*pi*this%nor
  else
    b = 0.0_wp
  end if

end subroutine compute_psi_vortring

!----------------------------------------------------------------------

!> Compute the velocity induced by a vortex ring in a prescribed position
!!
!! The velocity in the position is calculated considering the influece of
!! doublets
!!
!! WARNING: the velocity calculated, to be consistent with the formulation of 
!! the equations is multiplied by 4*pi, to obtain the actual velocity the 
!! result of the present subroutine MUST be DIVIDED by 4*pi
subroutine compute_vel_vortring(this, pos, uinf, vel)
  class(t_vortring), intent(inout) :: this
  real(wp), intent(in) :: pos(:)
  real(wp), intent(in) :: uinf(3)
  real(wp), intent(out) :: vel(3)

  real(wp) :: vdou(3)


  ! doublet ---
  call velocity_calc_doublet(this, vdou, pos)

  vel = vdou*this%idou

end subroutine compute_vel_vortring

!----------------------------------------------------------------------

subroutine compute_cp_vortring(this, elems, uinf)
  class(t_vortring), intent(inout) :: this
  type(t_elem_p), intent(in) :: elems(:)
  real(wp), intent(in) :: uinf(:)

  integer  :: i_stripe , i_c
  real(wp) :: dG_dt

  this%cp = 0.0_wp

  i_stripe = size(this%stripe_elem)
! WRONG(?): effective circulation (Gamma_i^{(e)} = Gamma_i - Gamma_{i-1})
!  must be used, s.t.  dG_dt = sum_{j=1}^{i} dGamma_j^{(e)} = dGamma_i
! dG_dt = 0.0_wp
! do i_c = 1 , i_stripe
!   dG_dt = dG_dt + elems(this%stripe_elem(i_c))%p%didou_dt 
! end do 

  dG_dt = this%didou_dt

  if ( i_stripe .gt. 1 ) then
    this%cp =   2.0_wp / norm2(uinf)**2.0_wp * &
          ( norm2(uinf - this%ub) * this%dy / this%area * &
               ( elems(this%id)%p%idou - this%stripe_elem(i_stripe-1)%p%idou ) + &
               dG_dt )
  else
    this%cp =   2.0_wp / norm2(uinf)**2.0_wp * &
          ( norm2(uinf - this%ub) * this%dy / this%area * &
                 elems(this%id)%p%idou + &
               dG_dt )
  end if

! old structures (it should be working as well). For stationary problems
! if ( this%stripe_1 .ne. 0 ) then
!   this%cp = 2.0_wp * norm2(uinf - this%ub) * this%dy / &
!           ( norm2(uinf)**2.0_wp * this%area ) * ( elems(this%id)%p%idou - elems(this%stripe_1)%p%idou )
! else
!   this%cp = 2.0_wp * norm2(uinf - this%ub) * this%dy / &
!           ( norm2(uinf)**2.0_wp * this%area ) *   elems(this%id)%p%idou
! end if

end subroutine compute_cp_vortring

!----------------------------------------------------------------------

end module mod_vortring
