!./\\\\\\\\\\\...../\\\......./\\\..../\\\\\\\\\..../\\\\\\\\\\\\\. 
!.\/\\\///////\\\..\/\\\......\/\\\../\\\///////\\\.\//////\\\////..
!..\/\\\.....\//\\\.\/\\\......\/\\\.\//\\\....\///.......\/\\\......
!...\/\\\......\/\\\.\/\\\......\/\\\..\////\\.............\/\\\......
!....\/\\\......\/\\\.\/\\\......\/\\\.....\///\\...........\/\\\......
!.....\/\\\......\/\\\.\/\\\......\/\\\.......\///\\\........\/\\\......
!......\/\\\....../\\\..\//\\\...../\\\../\\\....\//\\\.......\/\\\......
!.......\/\\\\\\\\\\\/....\///\\\\\\\\/..\///\\\\\\\\\/........\/\\\......
!........\///////////........\////////......\/////////..........\///.......
!!=========================================================================
!!
!! Copyright (C) 2018-2019 Davide   Montagnani, 
!!                         Matteo   Tugnoli, 
!!                         Federico Fonte
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
!!          Federico Fonte             <federico.fonte@outlook.com>
!!          Davide Montagnani       <davide.montagnani@gmail.com>
!!          Matteo Tugnoli                <tugnoli.teo@gmail.com>
!!=========================================================================

module mod_spline

use mod_param, only: &
  wp, max_char_len

implicit none

public :: t_spline , hermite_spline , deallocate_spline

private

! ----------------------------------------------------------------------
!> spline type
type :: t_spline
 real(wp) , allocatable  :: rr(:,:)  ! points to be interpolate ( n_d , n_points)
 real(wp) , allocatable  :: dp(:,:)  ! derivative at int. points( n_d , n_points)
 real(wp) , allocatable  :: t(:)     ! array of the param describing the spline
 real(wp) , allocatable  :: dt(:)    ! dt between two points of the spline
 character(max_char_len) :: e_bc(2)  ! bc: 'der', 'free'
 real(wp) , allocatable  :: d0(:)    ! derivative at the first point
 real(wp) , allocatable  :: d1(:)    ! derivative at the end point
end type t_spline

! ----------------------------------------------------------------------
contains
! ----------------------------------------------------------------------
!> build hermite_spline
subroutine hermite_spline ( spl , nelems , type_span , rr_spl )
  type(t_spline)               , intent(inout) :: spl
  integer                      , intent(in)    :: nelems
  character(max_char_len)      , intent(in)    :: type_span
  real(wp)       , allocatable , intent(out)   :: rr_spl(:,:)


  integer :: n_d , n_points , n_splines
  real(wp) , allocatable :: ll(:)  ! useless with this def of %t
  real(wp) , allocatable :: A(:,:) , b(:)

  ! lapack
  integer , allocatable :: ipiv(:)
  integer :: info

  ! spline reconstruction 
  integer , parameter :: n_points1 = 100  ! now hardcoded (pass as a default input)
  real(wp) , allocatable :: tv1(:)
  real(wp) , allocatable :: rr(:,:)

  integer :: i_r , i , j , n

  if ( allocated(spl% t) ) deallocate(spl% t)
  if ( allocated(spl%dt) ) deallocate(spl%dt)
  if ( allocated(spl%dp) ) deallocate(spl%dp)


  !> dimensions
  n_d      = size(spl%rr,2)  
  n_points = size(spl%rr,1) ; n_splines = n_points - 1

  ! === parameter vectors %t, %dt ===
  ! todo: check the sensitivity of the spline w.r.t. the param array t
  allocate( spl%t ( n_points  ) ) ; spl%t  = 0.0_wp
  allocate( spl%dt( n_splines ) ) ; spl%dt = 0.0_wp
  allocate( ll    ( n_splines ) ) ; ll     = 0.0_wp
  spl%t( 1) = 0.0_wp
  ll(    1) = norm2( spl%rr(2,:) - spl%rr(1,:) )
  spl%t( 2) = ll(1)
  spl%dt(1) = spl%t(2) - spl%t(1)
  do i = 2 , n_splines
    ll(    i  ) = norm2( spl%rr(i+1,:) - spl%rr(i,:) ) + ll(i-1)
    spl%t( i+1) = ll(i)
    spl%dt(i  ) = spl%t(i+1) - spl%t(i) 
  end do
  write(*,*) ' ll     : ' , ll 
  deallocate(ll)
 
  write(*,*) ' spl%t  : ' , spl%t
  write(*,*) ' spl%dt : ' , spl%dt

  ! === Linear system to solve for the derivatives of the spline ===
  !> one linear system per coorindate since the systems are not coupled
  allocate(spl%dp(n_points,n_d)) ; spl%dp = 0.0_wp
  allocate(A(n_points,n_points)) ; A = 0.0_wp
  allocate(b(n_points))          ; b = 0.0_wp

  do n = 1 , n_d

    !> Initialise (reset) A,b 
    A = 0.0_wp ; b = 0.0_wp

    !> inner points
    do i = 2 , n_points - 1

      A(i,i-1:i+1) = (/ spl%dt(i) ,  & 
                2.0_wp*(spl%dt(i-1)+spl%dt(i)) , &
                        spl%dt(i-1) /)
      b(i) = 3.0_wp * spl%dt(i-1)/spl%dt(i  ) * &
                    ( spl%rr(i+1,n) - spl%rr(i  ,n) ) + &
             3.0_wp * spl%dt(i  )/spl%dt(i-1) * & 
                    ( spl%rr(i  ,n) - spl%rr(i-1,n) ) 

    end do
 
    !> end points: boundary conditions 
    if ( trim( spl%e_bc(1) ) .eq. 'derivative' ) then
      A(1,1) = 1.0_wp
      b(1  ) = spl%d0(n)
    else
      write(*,*) ' hermite_spline. Wrong b.c. '
      write(*,*) ' Only "derivative" implemented so far. Stop.' ; stop
    end if 
    if ( trim( spl%e_bc(2) ) .eq. 'derivative' ) then
      A(n_points,n_points) = 1.0_wp
      b(n_points         ) = spl%d1(n)
    else
      write(*,*) ' hermite_spline. Wrong b.c. '
      write(*,*) ' Only "derivative" implemented so far. Stop.' ; stop
    end if 

    ! check ---
    do i = 1 , size(A,1)
      write(*,*) A(i,:) , '           ' , b(i)
    end do
    write(*,*)


    !> solve the linear system ( call to lapack dgesv )
    spl%dp(:,n) = b   ! intent(inout) in dgesv
    allocate(ipiv(n_points))
    call dgesv( n_points, 1 , A , n_points , &
             ipiv , spl%dp(:,n) , n_points , info ) 
    deallocate(ipiv)

  end do 

  ! check ---
  write(*,*)
  do i = 1 , size(spl%dp,1)
    write(*,*) spl%dp(i,:)
  end do

  ! === spline (non-uniform param) ===
  allocate(rr(n_points1*n_splines + 1,n_d))
  allocate(tv1(n_points1))
  i_r = 0
  do i = 1 , n_splines

    if ( i .ne. n_splines ) then
      tv1 = (/ ( dble(j) , j = 0,n_points1-1 ) /) / dble(n_points1)
    else
      if ( allocated(tv1) ) deallocate(tv1)
      allocate(tv1(n_points1+1))
      tv1 = (/ ( dble(j) , j = 0,n_points1   ) /) / dble(n_points1)
    end if
   
    do j = 1 , size(tv1)
      i_r = i_r + 1  ! update index
      rr(i_r,:) = hermite_p1(tv1(j)) * spl%rr(i  ,:) + &
                  hermite_p2(tv1(j)) * spl%rr(i+1,:) + &
                  hermite_d1(tv1(j)) * spl%dp(i  ,:) + &
                  hermite_d2(tv1(j)) * spl%dp(i+1,:)
    end do 

  end do

  deallocate(tv1)

  ! === spline ===
  allocate(rr_spl(nelems+1,3))








end subroutine hermite_spline

! ----------------------------------------------------------------------
!> hermite functions
!> h1: h1(0)=1 , h1(1)=0 , h1'(0)=0 , h1'(1)=0
function hermite_p1(t) result(h)
  real(wp) , intent(in)  :: t
  real(wp) :: h

  h =   2.0_wp * t**3.0_wp - 3.0_wp * t**2.0_wp + 1.0_wp

end function hermite_p1

!> h2: h2(0)=0 , h2(1)=1 , h2'(0)=0 , h2'(1)=0
function hermite_p2(t) result(h)
  real(wp) , intent(in)  :: t
  real(wp) :: h

  h = - 2.0_wp * t**3.0_wp + 3.0_wp * t**2.0_wp 

end function hermite_p2

!> h3: h3(0)=0 , h3(1)=0 , h3'(0)=1 , h3'(1)=0
function hermite_d1(t) result(h)
  real(wp) , intent(in)  :: t
  real(wp) :: h

  h = t**3.0_wp - 2.0_wp * t**2.0_wp + t

end function hermite_d1

!> h4: h4(0)=0 , h4(1)=0 , h4'(0)=0 , h4'(1)=1
function hermite_d2(t) result(h)
  real(wp) , intent(in)  :: t
  real(wp) :: h

  h = t**3.0_wp - t**2.0_wp

end function hermite_d2

! ----------------------------------------------------------------------
!> deallocate t_spline structure
subroutine deallocate_spline( spl )
  type(t_spline) , intent(inout) :: spl

  if ( allocated(spl%rr  ) ) deallocate(spl%rr  )
  if ( allocated(spl%dp  ) ) deallocate(spl%dp  )
  if ( allocated(spl%t   ) ) deallocate(spl%t   )
  if ( allocated(spl%dt  ) ) deallocate(spl%dt  )
  if ( allocated(spl%d0  ) ) deallocate(spl%d0  )
  if ( allocated(spl%d1  ) ) deallocate(spl%d1  )

end subroutine deallocate_spline

! ----------------------------------------------------------------------


end module mod_spline
