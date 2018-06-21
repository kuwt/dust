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

!> Module containing the specific subroutines for a single vortex line
module mod_vortline

use mod_param, only: &
  wp, pi, max_char_len, prev_tri, next_tri, prev_qua, next_qua

use mod_handling, only: &
  error, printout

use mod_doublet, only: &
  potential_calc_doublet , &
  velocity_calc_doublet

use mod_linsys_vars, only: &
  t_linsys

use mod_sim_param, only: &
  t_sim_param

use mod_math, only: &
  cross

use mod_c81, only: &
  t_aero_tab, interp_aero_coeff

!use mod_aero_elements, only: &
!  c_elem, t_elem_p

use mod_aeroel, only: &
  c_elem, c_pot_elem, c_vort_elem, c_impl_elem, c_expl_elem, &
  t_elem_p, t_pot_elem_p, t_vort_elem_p, t_impl_elem_p, t_expl_elem_p
!----------------------------------------------------------------------

implicit none

public :: t_vortline, initialize_vortline


!----------------------------------------------------------------------

type, extends(c_vort_elem) :: t_vortline
  real(wp) :: ver(3,2)
  real(wp) :: edge_vec(3)
  real(wp) :: edge_uni(3)
  real(wp) :: edge_len
contains

  procedure, pass(this) :: compute_vel      => compute_vel_vortline
  procedure, pass(this) :: calc_geo_data    => calc_geo_data_vortline

end type

character(len=*), parameter :: this_mod_name='mod_vortline'

real(wp) :: r_Rankine
real(wp) :: r_cutoff

!----------------------------------------------------------------------
contains
!----------------------------------------------------------------------

!> Initialize vortex line 
subroutine initialize_vortline(r_Rankine_in, r_cutoff_in)
 real(wp), intent(in) :: r_Rankine_in, r_cutoff_in

  r_Rankine = r_Rankine_in
  r_cutoff  = r_cutoff_in

end subroutine initialize_vortline

!----------------------------------------------------------------------

!> Compute the velocity induced by a vortex line in a prescribed position
!!
!! WARNING: the velocity calculated, to be consistent with the formulation of
!! the equations is multiplied by 4*pi, to obtain the actual velocity the
!! result of the present subroutine MUST be DIVIDED by 4*pi
subroutine compute_vel_vortline (this, pos, uinf, vel)
 class(t_vortline), intent(in) :: this
 real(wp), intent(in) :: pos(:)
 real(wp), intent(in) :: uinf(3)
 real(wp), intent(out) :: vel(3)

 real(wp) :: vdou(3)
 real(wp) :: av(3) , hv(3)
 real(wp) :: ai    , hi
 real(wp) :: R1 , R2
 real(wp) :: r_Ran

  !TODO: add far field approximations
  !radius_v = pos-this%cen
  !radius   = norm2(radius_v)
  if(associated(this%mag)) then

    ! use this%ver instead of its projection this%verp
    av = pos-this%ver(:,1)
    ai = sum(av*this%edge_uni(:))
    R1 = norm2(av)
    R2 = norm2(pos-this%ver(:,2))
    hv = av - ai*this%edge_uni(:)
    hi = norm2(hv)
    if ( hi .gt. this%edge_len*r_Rankine ) then
      ! (a/r+(s-a)/r)/h
      vdou = ( (this%edge_len-ai)/r2 + ai/r1 )/(hi**2.0_wp) * &
                      cross(this%edge_uni(:),hv)
    else
!     ! (a/r+(s-a)/r)* h/r_Rankine^2.
      if ( ( R1 .gt. this%edge_len*r_cutoff ) .and. &! avoid singularity
           ( R2 .gt. this%edge_len*r_cutoff )   ) then
        r_Ran = r_Rankine * this%edge_len
        vdou = ((this%edge_len-ai)/R2 + ai/R1)/(r_Ran**2.0_wp)* &
                         cross(this%edge_uni(:),hv)
!     else

      end if
    end if

    vel = vdou*this%mag

  else
    vel = 0.0_wp
  endif

end subroutine compute_vel_vortline

!----------------------------------------------------------------------

subroutine calc_geo_data_vortline(this, vert)
 class(t_vortline), intent(inout) :: this
 real(wp), intent(in) :: vert(:,:)

 integer :: sides, is, nsides
 real(wp):: nor(3), tanl(3)

  this%ver = vert


  ! center, for lines  is the mid-point
  this%cen =  sum ( this%ver(:,1:2),2 ) / 2.0_wp


  !! local tangent unit vector as in PANAIR
  !tanl = 0.5_wp * ( this%ver(:,nsides) + this%ver(:,1) ) - this%cen

  !this%tang(:,1) = tanl / norm2(tanl)
  !this%tang(:,2) = cross( this%nor, this%tang(:,1)  )


  !! vector connecting the two consecutive vertices:
  this%edge_vec(:) = this%ver(:,2) - this%ver(:,1)

  this%edge_len = norm2(this%edge_vec(:))

  this%edge_uni(:) = this%edge_vec(:) / this%edge_len

end subroutine calc_geo_data_vortline
!----------------------------------------------------------------------

end module mod_vortline
