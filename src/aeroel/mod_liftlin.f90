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

!> Module containing the specific subroutines for the lifting line 
!! type of aerodynamic elements
module mod_liftlin

use mod_param, only: &
  wp, pi

use mod_handling, only: &
  error

use mod_doublet, only: &
  potential_calc_doublet , &
  velocity_calc_doublet

use mod_linsys_vars, only: &
  t_linsys

use mod_c81, only: &
  t_aero_tab, interp_aero_coeff

use mod_aero_elements, only: &
  c_elem, t_elem_p
!----------------------------------------------------------------------

implicit none

public :: t_liftlin, update_liftlin, solve_liftlin


!----------------------------------------------------------------------

type, extends(c_elem) :: t_liftlin
  real(wp), allocatable :: tang_cen(:)
  real(wp), allocatable :: bnorm_cen(:)
  real(wp)              :: csi_cen
  integer               :: i_airfoil(2)
  real(wp)              :: chord
contains

  procedure, pass(this) :: build_row        => build_row_liftlin
  procedure, pass(this) :: build_row_static => build_row_static_liftlin
  procedure, pass(this) :: add_wake         => add_wake_liftlin
  procedure, pass(this) :: add_liftlin      => add_liftlin_liftlin
  procedure, pass(this) :: add_actdisk      => add_actdisk_liftlin
  procedure, pass(this) :: compute_pot      => compute_pot_liftlin
  procedure, pass(this) :: compute_vel      => compute_vel_liftlin
  procedure, pass(this) :: compute_psi      => compute_psi_liftlin
  procedure, pass(this) :: compute_cp       => compute_cp_liftlin
end type

character(len=*), parameter :: this_mod_name='mod_vortring'

integer :: it=0

!----------------------------------------------------------------------
contains
!----------------------------------------------------------------------

subroutine build_row_liftlin (this, elems, linsys, uinf, ie, ista, iend)
 class(t_liftlin), intent(inout) :: this
 type(t_elem_p), intent(in)       :: elems(:)
 type(t_linsys), intent(inout)    :: linsys
 real(wp), intent(in)             :: uinf(:)
 integer, intent(in)              :: ie
 integer, intent(in)              :: ista, iend
 character(len=*), parameter      :: this_sub_name='build_row_liftlin'
 
  call error(this_sub_name, this_mod_name, 'This was not supposed to &
  &happen, a team of professionals is underway to remove the evidence')
 
end subroutine build_row_liftlin

!----------------------------------------------------------------------

subroutine build_row_static_liftlin (this, elems, ll_elems, ad_elems, linsys, uinf, ie, ista, iend)
 class(t_liftlin), intent(inout) :: this
 type(t_elem_p), intent(in)       :: elems(:)
 type(t_elem_p), intent(in)       :: ll_elems(:)
 type(t_elem_p), intent(in)       :: ad_elems(:)
 type(t_linsys), intent(inout)    :: linsys
 real(wp), intent(in)             :: uinf(:)
 integer, intent(in)              :: ie
 integer, intent(in)              :: ista, iend
 character(len=*), parameter      :: this_sub_name='build_row_static_liftlin'
 
  call error(this_sub_name, this_mod_name, 'This was not supposed to &
  &happen, a team of professionals is underway to remove the evidence')

end subroutine build_row_static_liftlin

!----------------------------------------------------------------------

subroutine add_wake_liftlin (this, wake_elems, impl_wake_ind, linsys, uinf, &
                             ie, ista, iend)
 class(t_liftlin), intent(inout) :: this
 type(t_elem_p), intent(in)      :: wake_elems(:)
 integer, intent(in)             :: impl_wake_ind(:,:)
 type(t_linsys), intent(inout)   :: linsys
 real(wp), intent(in)            :: uinf(:)
 integer, intent(in)             :: ie
 integer, intent(in)             :: ista
 integer, intent(in)             :: iend
 character(len=*), parameter      :: this_sub_name='add_wake_liftlin'
 
  call error(this_sub_name, this_mod_name, 'This was not supposed to &
  &happen, a team of professionals is underway to remove the evidence')

end subroutine add_wake_liftlin

!----------------------------------------------------------------------

subroutine add_liftlin_liftlin (this, ll_elems, linsys, uinf, &
                             ie, ista, iend)
 class(t_liftlin), intent(inout) :: this
 type(t_elem_p), intent(in)      :: ll_elems(:)
 type(t_linsys), intent(inout)   :: linsys
 real(wp), intent(in)            :: uinf(:)
 integer, intent(in)             :: ie
 integer, intent(in)             :: ista
 integer, intent(in)             :: iend
 character(len=*), parameter      :: this_sub_name='add_liftlin_liftlin'
 
  call error(this_sub_name, this_mod_name, 'This was not supposed to &
  &happen, a team of professionals is underway to remove the evidence')

end subroutine add_liftlin_liftlin

!----------------------------------------------------------------------

subroutine add_actdisk_liftlin (this, ad_elems, linsys, uinf, &
                             ie, ista, iend)
 class(t_liftlin), intent(inout) :: this
 type(t_elem_p), intent(in)      :: ad_elems(:)
 type(t_linsys), intent(inout)   :: linsys
 real(wp), intent(in)            :: uinf(:)
 integer, intent(in)             :: ie
 integer, intent(in)             :: ista
 integer, intent(in)             :: iend
 character(len=*), parameter      :: this_sub_name='add_actdisk_liftlin'
 
  call error(this_sub_name, this_mod_name, 'This was not supposed to &
  &happen, a team of professionals is underway to remove the evidence')

end subroutine add_actdisk_liftlin

!----------------------------------------------------------------------

subroutine compute_pot_liftlin (this, A, b, pos,i,j)
 class(t_liftlin), intent(inout) :: this
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

end subroutine compute_pot_liftlin

!----------------------------------------------------------------------

subroutine compute_psi_liftlin (this, A, b, pos, nor, i, j )
 class(t_liftlin), intent(inout) :: this
 real(wp), intent(out) :: A
 real(wp), intent(out) :: b(3)
 real(wp), intent(in) :: pos(:)
 real(wp), intent(in) :: nor(:)
 integer , intent(in) :: i , j

 real(wp) :: vdou(3) 

  call velocity_calc_doublet(this, vdou, pos)

  A = sum(vdou * nor)


  !  b = ... (from boundary conditions)
  !TODO: consider moving this outside
  if ( i .eq. j ) then
    b =  4.0_wp*pi*this%nor
  else
    b = 0.0_wp
  end if

end subroutine compute_psi_liftlin

!----------------------------------------------------------------------
subroutine compute_vel_liftlin (this, pos, uinf, vel)
 class(t_liftlin), intent(inout) :: this
 real(wp), intent(in) :: pos(:)
 real(wp), intent(in) :: uinf(3)
 real(wp), intent(out) :: vel(3)

 real(wp) :: vdou(3)


  ! doublet ---
  call velocity_calc_doublet(this, vdou, pos)
 

  vel = vdou*this%idou


end subroutine compute_vel_liftlin

!----------------------------------------------------------------------
subroutine compute_cp_liftlin (this, elems, uinf)
 class(t_liftlin), intent(inout) :: this
 type(t_elem_p), intent(in) :: elems(:)
 real(wp), intent(in) :: uinf(:)

 character(len=*), parameter      :: this_sub_name='add_liftlin_liftlin'
 
  call error(this_sub_name, this_mod_name, 'This was not supposed to &
  &happen, a team of professionals is underway to remove the evidence')

end subroutine compute_cp_liftlin 

!----------------------------------------------------------------------

subroutine update_liftlin(elems_ll, linsys)
 type(t_elem_p), intent(inout) :: elems_ll(:)
 type(t_linsys), intent(inout) :: linsys

 real(wp), allocatable :: res_temp(:)

  it = it + 1
  !DEBUG
  write(*,*) 'iteration ',it
  
  !HERE extrapolate the solution before the linear system
  if (it .gt. 2) then
    allocate(res_temp(size(linsys%res_expl,1)))
    res_temp = linsys%res_expl(:,1)
    linsys%res_expl(:,1) = 2.0_wp*res_temp - linsys%res_expl(:,2)
    linsys%res_expl(:,2) = res_temp
    deallocate(res_temp)
  else
    linsys%res_expl(:,2) = linsys%res_expl(:,1)
  endif

end subroutine update_liftlin

!----------------------------------------------------------------------

subroutine solve_liftlin(elems_ll, elems_tot, elems_wake,  uinf, airfoil_data)
 type(t_elem_p), intent(inout) :: elems_ll(:)
 type(t_elem_p), intent(in)    :: elems_tot(:)
 type(t_elem_p), intent(in)    :: elems_wake(:)
 real(wp), intent(in)          :: uinf(3)
 type(t_aero_tab),  intent(in) :: airfoil_data(:)

 integer :: i_l, j, ic
 real(wp) :: vel(3), v(3), up(3)
 real(wp), allocatable :: vel_w(:,:)
 real(wp) :: unorm, alpha, mach, re
 real(wp) :: cl
 real(wp), allocatable :: aero_coeff(:)
 real(wp), allocatable :: dou_temp(:)
 real(wp) :: damp=2.0_wp
 real(wp), parameter :: toll=1e-6_wp
 real(wp) :: diff
 
 !TODO: is missing the velocity of the wake
 allocate(dou_temp(size(elems_ll)))
 allocate(vel_w(3,size(elems_ll)))
 
 do i_l = 1,size(elems_ll)
   vel_w(:,i_l) = 0.0_wp
   do j = 1,size(elems_wake)
     call elems_wake(j)%p%compute_vel(elems_ll(i_l)%p%cen,uinf,v)
     vel_w(:,i_l) = vel_w(:,i_l) + v
   enddo
 enddo
 
 vel_w = vel_w/(4.0_wp*pi)

 !Crack the starting guess
 !do i_l = 1,size(elems_ll)
 !  elems_ll(i_l)%p%idou = -1
 !enddo


 !Calculate the induced velocity on the airfoil
 do ic = 1,100
   diff = 0.0_wp
   do i_l = 1,size(elems_ll)
     vel = 0.0_wp
     do j = 1,size(elems_tot)
       call elems_tot(j)%p%compute_vel(elems_ll(i_l)%p%cen,uinf,v)
       vel = vel + v
     enddo
     select type(el => elems_ll(i_l)%p)
     type is(t_liftlin)
       vel = vel/(4.0_wp*pi) + uinf - el%ub +vel_w(:,i_l)
       !vel = uinf - el%ub
       up = vel-el%bnorm_cen*sum(el%bnorm_cen*vel)
       !Employing the free stream velocity to get into tables
       unorm = norm2((uinf-el%ub) - el%bnorm_cen*sum(el%bnorm_cen*(uinf-el%ub)))
       !unorm = norm2(up)
       alpha = atan2(sum(up*el%nor), sum(up*el%tang_cen))
       alpha = alpha * 180.0_wp/pi
       !TODO: fix these parameters which are still hard-coded
       mach = 0.0_wp
       re = 1000000.0_wp
       !if (ic .le. 1) then
       !  cl = 2*pi*alpha*pi/180.0_wp
       !else
       call interp_aero_coeff ( airfoil_data ,  &
                      el%csi_cen, el%i_airfoil , (/alpha, mach, re/) , aero_coeff )
       cl = aero_coeff(1)
       !endif
       dou_temp(i_l) = - 0.5_wp * unorm * cl * el%chord
       diff = max(diff,abs(elems_ll(i_l)%p%idou-dou_temp(i_l))) 
     end select
   enddo
   damp = 5.0_wp
  
   do i_l = 1,size(elems_ll)
     elems_ll(i_l)%p%idou = ( dou_temp(i_l)+ damp*elems_ll(i_l)%p%idou )/(1.0_wp+damp)
     !elems_ll(i_l)%p%idou = ( (damp)*dou_temp(i_l) + (1.0_wp-damp)*elems_ll(i_l)%p%idou )
     !elems_ll(i_l)%p%idou = dou_temp(i_l)
   enddo 
   
   !DEBUG:
   write(*,*) 'diff', diff
   if(diff .le. toll) exit

 enddo

 !Rough cp compuation ! Only steady ...
 do i_l = 1,size(elems_ll)
   elems_ll(i_l)%p%cp =  2.0_wp / norm2(uinf)**2.0_wp * &
            ( norm2(uinf - elems_ll(i_l)%p%ub) ) * elems_ll(i_l)%p%idou 
 end do

 !DEBUG:
 write(*,*) 'iterations: ',ic
 write(*,*) 'diff',diff
 deallocate(dou_temp, vel_w)


 !Get the angle of attack, as well as the other parameters

 !Get into the tables to obtain the loads

 !Get the vorticity of the element

end subroutine solve_liftlin

!----------------------------------------------------------------------

end module mod_liftlin
