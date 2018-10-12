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


!> Module to handle the multipole and local expansions, and their 
!! manipulation
module mod_multipole

use mod_param, only: &
  wp, nl, pi, max_char_len

use mod_sim_param, only: &
  t_sim_param

use mod_math, only: &
  cross

use mod_handling, only: &
  error, warning, info, printout, dust_time, t_realtime

use mod_vortpart, only: &
  t_vortpart, t_vortpart_p

!----------------------------------------------------------------------

implicit none

public :: t_multipole, t_polyexp, t_ker_der

private

!----------------------------------------------------------------------

!> Type containing the multipole expansion relative to a single cell
type :: t_multipole

  !> Multipole expansion coefficients (velocity)
  real(wp), allocatable :: a(:,:)

  !> Local expansion coefficients (velocity)
  real(wp), allocatable :: b(:,:)

  !> Local expansion coefficients (gradient)
  real(wp), allocatable :: c(:,:,:)

contains
  
  !> Initialize the data
  procedure, pass(this) :: init => init_multipole
  
  !> Calculate the coefficients of the multipole in the leaf
  procedure, pass(this) :: leaf_M => leaf_M_multipole

  !> Pass the multipole coefficients to the parents
  procedure, pass(this) :: M2M => M2M_multipole

  !> Convert the multipole coefficients of the interacting cells into
  !! local expansion coefficients
  procedure, pass(this) :: M2L => M2L_multipole

  !> Inherit the local expansion coefficients from parents cells
  procedure, pass(this) :: L2L => L2L_multipole

end type

!----------------------------------------------------------------------

!> Type containing the derivatives of kernel evalued at a certain distance
type :: t_ker_der
  
  !> Store the derivatives of the kernel evalued between two points
  real(wp), allocatable :: D(:,:)

  real(wp), allocatable :: Dc(:,:,:)

contains

  procedure, pass(this) :: compute_der => compute_der_ker

end type 

!----------------------------------------------------------------------

!> Polynomial expansion tools
type :: t_polyexp
  
  !> Maximum degree of the polynomial expansion
  integer :: degree

  !> Number of monomials of the polynomial expansion
  integer :: n_mon
  
  !> Number of monomials, at all the different degrees
  integer, allocatable :: n_mon_d(:)
  
  !> For a given set of degrees in the three dimensions returns the index
  !! in the unrolled 
  integer, allocatable :: idx(:,:,:) 
  
  integer, allocatable :: pwr(:,:)

  integer, allocatable :: fact(:)

  integer, allocatable :: nfact(:,:,:)

contains

  procedure, pass(this) :: set_degree => set_degree_polyexp  

  procedure, pass(this) :: nbinom => nbinom_polyexp

  !procedure, pass(this) :: nfact => nfact_polyexp

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

  allocate(this%a(3,polyexp%n_mon))
  allocate(this%b(3,polyexp%n_mon))
  allocate(this%c(3,3,polyexp%n_mon))


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
      this%a(:,m) = this%a(:,m) + &
      parts(i)%p%mag*product((parts(i)%p%cen-cen)**pexp%pwr(:,m))* &
      parts(i)%p%dir
    enddo
  enddo

end subroutine

!----------------------------------------------------------------------

!this should be called for all the children
subroutine M2M_multipole(this, cen, child, child_cen, pexp)
 class(t_multipole) :: this
 real(wp), intent(in) :: cen(3)
 type(t_multipole), intent(in) :: child
 real(wp), intent(in) :: child_cen(3)
 type(t_polyexp), intent(in) :: pexp

 integer :: m, idx(3), s
 integer :: is, js, ks

  do m=1,size(this%a,2) 
    idx = pexp%pwr(:,m)
    do ks = 0,idx(3); do js = 0,idx(2); do is = 0,idx(1)
      s = pexp%idx(is,js,ks)

      this%a(:,m) = this%a(:,m) + &
      pexp%nbinom(idx,(/is,js,ks/))* &
      product((child_cen-cen)**(idx-(/is, js, ks/)))*child%a(:,s)

    enddo; enddo; enddo
  enddo

end subroutine

!----------------------------------------------------------------------

subroutine M2L_multipole(this, ker_der, pexp, pexp_der, multipol_int)
 class(t_multipole) :: this
 type(t_ker_der), intent(in) :: ker_der
 type(t_polyexp), intent(in) :: pexp
 type(t_polyexp), intent(in) :: pexp_der
 type(t_multipole), intent(in) :: multipol_int

 real(wp) :: sum_v(3), sum_g(3,3), mult

 integer :: n, m, idx(3), idx_der
 
  !Subdivide all the degrees, for the first ones expand both the velocity
  !and the gradient, for the last one only the velocity

  do m = 1, pexp%n_mon
    sum_g = 0.0_wp
    sum_v = 0.0_wp
    do n = 1, pexp%n_mon
      idx = pexp%pwr(:,m)+pexp%pwr(:,n)
      idx_der = pexp_der%idx(idx(1),idx(2),idx(3))

      mult = pexp_der%nfact(idx(1),idx(2),idx(3))/( &
                pexp%nfact(pexp%pwr(1,m),pexp%pwr(2,m),pexp%pwr(3,m))*&
                pexp%nfact(pexp%pwr(1,n),pexp%pwr(2,n),pexp%pwr(3,n)))

      sum_v = sum_v + mult * cross(ker_der%D(:,idx_der), multipol_int%a(:,n))
      
      sum_g(:,1) = sum_g(:,1) + mult * &
                   cross(ker_der%Dc(:,1,idx_der), multipol_int%a(:,n))
      sum_g(:,2) = sum_g(:,2) + mult * &
                   cross(ker_der%Dc(:,2,idx_der), multipol_int%a(:,n))
      sum_g(:,3) = sum_g(:,3) + mult * &
                   cross(ker_der%Dc(:,3,idx_der), multipol_int%a(:,n))
      !sum1 = sum1 +  &
      !pexp_der%nfact(idx(1),idx(2),idx(3))/( &
      !pexp%nfact(pexp%pwr(1,m),pexp%pwr(2,m),pexp%pwr(3,m))*&
      !pexp%nfact(pexp%pwr(1,n),pexp%pwr(2,n),pexp%pwr(3,n))) * &
      ! cross(ker_der%D(:,idx_der), multipol_int%a(:,n))
    enddo
    this%b(:,m) = this%b(:,m) + real((-1)**(sum(pexp%pwr(:,m))),wp)*sum_v
    this%c(:,:,m) = this%c(:,:,m) + real((-1)**(sum(pexp%pwr(:,m))),wp)*sum_g
  enddo 

end subroutine

!----------------------------------------------------------------------

subroutine L2L_multipole(this, cen, parent, parent_cen, pexp) 
 class(t_multipole) :: this
 real(wp), intent(in) :: cen(3)
 type(t_multipole), intent(in) :: parent
 real(wp), intent(in) :: parent_cen(3)
 type(t_polyexp), intent(in) :: pexp

 integer :: m, idx(3), s
 integer :: is, js, ks
 real(wp) :: mult
  
  do m=1,size(this%b,2) 
    idx = pexp%pwr(:,m)
    do ks = idx(3),pexp%degree; do js = idx(2),pexp%degree-ks; do is = idx(1),pexp%degree-ks-js
      s = pexp%idx(is,js,ks)

      mult = pexp%nbinom((/is,js,ks/),idx)* &
          product((cen-parent_cen)**((/is, js, ks/) - idx))

      this%b(:,m) = this%b(:,m) + mult*parent%b(:,s)
      this%c(:,:,m) = this%c(:,:,m) + mult*parent%c(:,:,s)
    enddo; enddo; enddo
  enddo

end subroutine

!----------------------------------------------------------------------

subroutine compute_der_ker(this,diff,delta,pexp)
 class(t_ker_der) :: this
 real(wp), intent(in) :: diff(3)
 real(wp), intent(in) :: delta
 type(t_polyexp), intent(in) :: pexp
 
 integer :: i, j, k
 real(wp) :: kmod, Rnorm2
 real(wp) :: sum1, sum2
 real(wp), allocatable :: bk(:)


  allocate(this%D(3,pexp%n_mon_d(pexp%degree-1)))
  allocate(this%Dc(3,3,pexp%n_mon_d(pexp%degree-2)))
  allocate(bk(pexp%n_mon))
  this%D = 0.0_wp
  this%Dc = 0.0_wp
  bk = 0.0_wp
  Rnorm2 = sum(diff**2) + delta**2 

  do k=0,pexp%degree;  do j=0,pexp%degree-k; do i=0,pexp%degree-j-k;
    kmod = real(k+j+i,wp)

    if (all((/i,j,k/).eq.0)) then
      bk(pexp%idx(i,j,k)) = 1.0_wp/(4.0_wp*pi*sqrt(Rnorm2))
    else

        sum1 = 0.0_wp; sum2 = 0.0_wp

        if(i-1 .ge. 0) sum1 = sum1 + diff(1)*bk(pexp%idx(i-1,j,k))
        if(j-1 .ge. 0) sum1 = sum1 + diff(2)*bk(pexp%idx(i,j-1,k))
        if(k-1 .ge. 0) sum1 = sum1 + diff(3)*bk(pexp%idx(i,j,k-1))

        if(i-2 .ge. 0) sum2 = sum2 + bk(pexp%idx(i-2,j,k))
        if(j-2 .ge. 0) sum2 = sum2 + bk(pexp%idx(i,j-2,k))
        if(k-2 .ge. 0) sum2 = sum2 + bk(pexp%idx(i,j,k-2))


        bk(pexp%idx(i,j,k)) = 1.0_wp/(Rnorm2*kmod)*&
                     ((2.0_wp*kmod-1.0_wp)*sum1 - (kmod-1.0_wp)*sum2) 
    endif
  enddo; enddo; enddo

  do k=0,pexp%degree-1;  do j=0,pexp%degree-1-k; do i=0,pexp%degree-1-j-k;
    this%D(1,pexp%idx(i,j,k)) = - (1+i)*bk(pexp%idx(i+1,j,k))
    this%D(2,pexp%idx(i,j,k)) = - (1+j)*bk(pexp%idx(i,j+1,k))
    this%D(3,pexp%idx(i,j,k)) = - (1+k)*bk(pexp%idx(i,j,k+1))
  enddo; enddo; enddo

  do k=0,pexp%degree-2;  do j=0,pexp%degree-2-k; do i=0,pexp%degree-2-j-k;
    this%Dc(:,1,pexp%idx(i,j,k)) = - (1+i)*this%D(:,pexp%idx(i+1,j,k))
    this%Dc(:,2,pexp%idx(i,j,k)) = - (1+j)*this%D(:,pexp%idx(i,j+1,k))
    this%Dc(:,3,pexp%idx(i,j,k)) = - (1+k)*this%D(:,pexp%idx(i,j,k+1))
  enddo; enddo; enddo
end subroutine

!----------------------------------------------------------------------

subroutine set_degree_polyexp(this, deg)
 class(t_polyexp) :: this
 integer, intent(in) :: deg

 integer :: i, j, k, ipol, o

  this%degree = deg

  allocate(this%idx(0:deg,0:deg,0:deg))
  allocate(this%n_mon_d(0:deg))

  ipol = 0
  do o = 0,deg
   do i=o,0,-1
    do j=o-i,0,-1
        k=o-j-i
        ipol = ipol+1 
        this%idx(i,j,k) = ipol
  enddo; enddo;
    this%n_mon_d(o) = ipol
  enddo
  
  this%n_mon = ipol

  allocate(this%pwr(3,this%n_mon))

  ipol = 0
  do o = 0,deg
   do i=o,0,-1
    do j=o-i,0,-1
        k=o-j-i
        ipol = ipol+1 
        this%pwr(:,ipol) = (/i,j,k/)
  enddo; enddo;
  enddo
  
  allocate(this%fact(0:deg))
  this%fact(0) = 1
  do o = 1,deg
    this%fact(o) = o*this%fact(o-1)
  enddo

  allocate(this%nfact(0:deg,0:deg,0:deg))
  do i=0,deg; do j=0,deg; do k=0,deg
    this%nfact(i,j,k) = nfact_polyexp(this,(/i,j,k/))
  enddo; enddo; enddo

end subroutine set_degree_polyexp

!----------------------------------------------------------------------

function nbinom_polyexp(this, m, s) result(nbinom)
 class(t_polyexp) :: this
 integer, intent(in) :: m(3), s(3)

 integer :: nbinom
 integer :: i, d

  nbinom = 1
  do d=1,3
    nbinom = nbinom* &
             (this%fact(m(d)) / &
             ( this%fact(s(d))*this%fact(m(d)-s(d)) ) ) 
  enddo

end function nbinom_polyexp

!----------------------------------------------------------------------


function nfact_polyexp(this, m) result(nfact)
 class(t_polyexp) :: this
 integer, intent(in) :: m(3)

 integer :: nfact
 integer :: i, d

  nfact = 1
  do d=1,3
    nfact = nfact * this%fact(m(d))  
  enddo

end function nfact_polyexp

!----------------------------------------------------------------------

end module mod_multipole
