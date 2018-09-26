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


!> Module to handle the octree grid
module mod_multipole

use mod_param, only: &
  wp, nl, pi, max_char_len

use mod_sim_param, only: &
  t_sim_param

use mod_handling, only: &
  error, warning, info, printout, dust_time, t_realtime

use mod_vortpart, only: &
  t_vortpart, t_vortpart_p

!----------------------------------------------------------------------

implicit none

public :: t_multipole, t_polyexp

private

!----------------------------------------------------------------------

!> Type containing all the data relative to a single octree cell
type :: t_multipole

  real(wp), allocatable :: a(:,:)

  real(wp), allocatable :: b(:,:)

  contains

  procedure, pass(this) :: init => init_multipole

  procedure, pass(this) :: leaf_M => leaf_M_multipole

  procedure, pass(this) :: M2M => M2M_multipole

  procedure, pass(this) :: M2L => M2L_multipole

  procedure, pass(this) :: L2L => L2L_multipole

end type

!----------------------------------------------------------------------

type :: t_polyexp

  integer :: degree

  integer :: n_mon

  integer, allocatable :: idx(:,:,:) 
  
  integer, allocatable :: pwr(:,:)

contains

  procedure, pass(this) :: set_degree => set_degree_polyexp  

end type

!----------------------------------------------------------------------
character(len=*), parameter :: this_mod_name='mod_multipole'
character(len=max_char_len) :: msg

!----------------------------------------------------------------------
contains
!----------------------------------------------------------------------

subroutine init_multipole(this, polyexp)
 class(t_multipole) :: this
 type(t_polyexp), intent(in) :: polyexp

  allocate(this%a(1,polyexp%n_mon))
  allocate(this%b(1,polyexp%n_mon))


end subroutine

!----------------------------------------------------------------------

subroutine leaf_M_multipole(this, cen, parts, pexp)
 class(t_multipole) :: this
 real(wp), intent(in) :: cen(3)
 type(t_vortpart_p), intent(in) :: parts(:)
 type(t_polyexp), intent(in) :: pexp

 integer :: i, m

  this%a = 0.0_wp
  do m=1,size(this%a,2) 
    do i=1,size(parts)
      this%a(1,m) = this%a(1,m) + &
      parts(i)%p%mag*product((parts(i)%p%cen-cen)**pexp%pwr(:,m))
    enddo
  enddo

end subroutine

!----------------------------------------------------------------------

subroutine M2M_multipole(this)
 class(t_multipole) :: this

end subroutine

!----------------------------------------------------------------------

subroutine M2L_multipole(this)
 class(t_multipole) :: this

end subroutine

!----------------------------------------------------------------------

subroutine L2L_multipole(this)
 class(t_multipole) :: this

end subroutine

!----------------------------------------------------------------------

subroutine set_degree_polyexp(this, deg)
 class(t_polyexp) :: this
 integer, intent(in) :: deg

 integer :: i, j, k, ipol, o

  this%degree = deg

  allocate(this%idx(0:deg,0:deg,0:deg))

  ipol = 0
!  do o = 0,deg
!  !do i=0,deg
!  !  do j=0,deg-i
!  !    do k=0,deg-i-j
!  do k=0,o
!    do j=0,o-k
!      !do i=0,o-k-j
!       i=o-k-j
!        ipol = ipol+1 
!        this%idx(i,j,k) = ipol
!  enddo; enddo; enddo;
!  enddo
  do o = 0,deg
  !do i=0,deg
  !  do j=0,deg-i
  !    do k=0,deg-i-j
  write(*,*) '  order ',o
  do i=o,0,-1
    do j=o-i,0,-1
      !do i=0,o-k-j
        k=o-j-i
        ipol = ipol+1 
        this%idx(i,j,k) = ipol
  write(*,*) i,j,k
  enddo; enddo;
  enddo
  
  this%n_mon = ipol

  allocate(this%pwr(3,this%n_mon))

  ipol = 0
  do i=0,deg
    do j=0,deg-i
      do k=0,deg-i-j
        ipol = ipol+1 
        this%pwr(:,ipol) = (/i,j,k/)
  enddo; enddo; enddo;


end subroutine set_degree_polyexp

!----------------------------------------------------------------------

end module mod_multipole
