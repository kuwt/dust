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
!! Copyright (C) 2018-2022 Politecnico di Milano,
!!                           with support from A^3 from Airbus
!!                    and  Davide   Montagnani,
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
!!          Federico Fonte
!!          Davide Montagnani
!!          Matteo Tugnoli
!!=========================================================================

!> Module containing some mathematical utilities
module mod_math

use mod_param, only: &
  wp

use mod_handling, only: &
  error, warning

implicit none

public :: dot, cross , linear_interp , compute_qr, rotation_vector_combination, sort_vector_real

private

interface linear_interp
  module procedure  linear_interp_vector, &
                    linear_interp_array
end interface linear_interp

character(len=*), parameter :: this_mod_name='mod_math'

contains

pure function dot(a, b)
  real(wp)                           :: dot
  real(wp), dimension(:), intent(in) :: a, b

  dot = sum(a * b)

end function dot


! ----------------------------------------------------------------------

pure function cross(a, b)
  real(wp), dimension(3) :: cross
  real(wp), dimension(3), intent(in) :: a, b

  cross(1) = a(2) * b(3) - a(3) * b(2)
  cross(2) = a(3) * b(1) - a(1) * b(3)
  cross(3) = a(1) * b(2) - a(2) * b(1)

end function cross

! ----------------------------------------------------------------------

subroutine linear_interp_vector(val_vec , t_vec , t , val)
  real(wp), intent(in)        :: val_vec(:)
  real(wp), intent(in)        :: t_vec(:)
  real(wp), intent(in)        :: t
  real(wp), intent(out)       :: val
  integer                     :: it , nt
  character(len=*), parameter :: this_sub_name='linear_interp_vector'

  nt = size(t_vec)
  ! Check dimensions -----
  if ( size(val_vec) .ne. nt ) then
    call error(this_sub_name, this_mod_name, 'Different sizes for x and y &
                                  &data vector provided for interpolation')
  end if

  ! Check if t \in [ minval(t_vec) , maxval(t_vec) ]
  if ( t .lt. minval(t_vec) ) then
    call error(this_sub_name, this_mod_name, 'x value requested to be &
          &interpolated is lower than the minimum of the interpolation data')
  end if
  if ( t .gt. maxval(t_vec) ) then
    call error(this_sub_name, this_mod_name, 'x value requested to be &
          &interpolated is higher than the minimum of the interpolation data')
  end if

  do it = 1 , nt-1

    if ( ( t .ge. t_vec( it ) ) .and. ( t .le. t_vec( it+1 ) ) ) then

      val = val_vec(it) + (t-t_vec(it))/(t_vec(it+1)-t_vec(it)) * &
                                  ( val_vec(it+1)-val_vec(it) )

    end if

  end do

end subroutine linear_interp_vector

! ----------------------------------------------------------------------

subroutine linear_interp_array( val_arr , t_vec , t , val )
  real(wp) , intent(in) :: val_arr(:,:)
  real(wp) , intent(in) :: t_vec(:)
  real(wp) , intent(in) :: t
  real(wp) , allocatable , intent(out) :: val(:)

  integer :: it , nt
  character(len=*), parameter :: this_sub_name='linear_interp_array'

  nt = size(t_vec)
  allocate(val(size(val_arr,1)))

  ! Check dimensions -----
  if ( size(val_arr,2) .ne. nt ) then
    call error(this_sub_name, this_mod_name, 'Different sizes for x and y &
                                  &data vector provided for interpolation')
  end if

  ! Check if t \in [ minval(t_vec) , maxval(t_vec) ]
  if ( t .lt. minval(t_vec) ) then
    call error(this_sub_name, this_mod_name, 'x value requested to be &
          &interpolated is lower than the minimum of the interpolation data')
  end if
  if ( t .gt. maxval(t_vec) ) then
    call error(this_sub_name, this_mod_name, 'x value requested to be &
          &interpolated is higher than the minimum of the interpolation data')
  end if

  do it = 1 , nt-1

    if ( ( t .ge. t_vec( it ) ) .and. ( t .le. t_vec( it+1 ) ) ) then

      val = val_arr(:,it) + (t-t_vec(it))/(t_vec(it+1)-t_vec(it)) * &
                                  ( val_arr(:,it+1)-val_arr(:,it) )

    end if

  end do

end subroutine linear_interp_array

! ----------------------------------------------------------------------

subroutine compute_qr ( A , Q , R )
 real(wp) , intent(inout) ::  A(:,:)
 real(wp) , allocatable , intent(inout) :: Q(:,:) , R(:,:)

 integer :: m , n , i_m

 ! lapack dgeqrf routine ----
 real(wp) , allocatable :: tau(:) , work(:)
 integer :: lwork , info

 ! tmp matrices to get Q matrix -----
 real(wp) , allocatable :: H(:,:) , v(:,:) , eye(:,:)

 character(len=*), parameter :: this_sub_name='compute_qr'


  ! input check and warnings
  if ( allocated(Q) ) then
    call warning(this_sub_name, this_mod_name, ' Q was already allocated.&
                                          & Deallocated and re-allocated')
    deallocate(Q)
  end if
  if ( allocated(R) ) then
    call warning(this_sub_name, this_mod_name, ' R was already allocated.&
                                          & Deallocated and re-allocated')
    deallocate(R)
  end if

  !
  m = size(A,1)
  n = size(A,2)

  ! qr factorisation of matrix B
  allocate(   tau(min(m,n)) ) ;    tau = 0.0_wp
  allocate(  work(    n   ) ) ;   work = 0.0_wp
  lwork = n       ! <-- its size should be .ge. n*nb
                  ! with nb = optimal blocksize (???)

#if (DUST_PRECISION==1)
  call sgeqrf( m , n , A , m , tau , work , lwork , info )
#elif(DUST_PRECISION==2)
  call dgeqrf( m , n , A , m , tau , work , lwork , info )
#endif /*DUST_PRECISION*/

  ! build Q , R matrices
  ! allocataion and initialisation to 0.0 (useless for Q)
  allocate( Q(m,m) , R(m,n) ) ; R = 0.0_wp

  ! R upper triangular matrix (initialised to zero!)
  do i_m = 1 , m
    R( i_m , i_m : size(A,2) ) = A( i_m , i_m : size(A,2) )
  end do

  ! Q square unitary matrix
  allocate( H(m,m) , v(m,1) ) ! ; H = 0.0_wp  ; v(:,1) = 0.0_wp
  allocate( eye(m,m) ) ; eye = 0.0_wp ; do i_m = 1,m ; eye(i_m,i_m) = 1.0_wp ; end do
  Q = eye ! Initialisation
  do i_m = 1 , min(m,n)
    v(:,1) = 0.0_wp ; v(i_m,1) = 1.0_wp ; v(i_m+1:m,1) = A(i_m+1:m,i_m) ;
    H = eye - tau(i_m) * matmul( v , transpose(v) )
    Q = matmul( Q , H )
  end do

  deallocate(H,v,eye)
  deallocate(tau,work)

end subroutine compute_qr

!----------------------------------------------------------------
!> Combination of two rotations, provided as rotation vectors r1, r2
! x: original vector
! y: rotated vector
! R1: rotation matrix from rotation vector r1, r1 -> R1
! R2: rotation matrix from rotation vector r2, r2 -> R2
!   y = R1 * R2 * x = R * x
! Inputs: r1, r2
! Outputs: r, th, n
! R <- r = th * n, with |n| = 1
!
subroutine rotation_vector_combination( r1, r2, r, th, n )
  real(wp),              intent(in)  :: r1(3), r2(3)
  real(wp),              intent(out) :: r(3)
  real(wp),              intent(out) :: n(3)
  real(wp),              intent(out) :: th

  real(wp) :: n1(3), n2(3), q1(4), q2(4), qs(4)
  real(wp) :: th1, th2, sin_t_2, cos_t_2

  if ( ( size(r1) .ne. 3 ) .or. ( size(r2) .ne. 3 ) ) then
    write(*,*) ' Wrong input dimensions in rotation_vector_combinations. Stop '
  end if

  th1 = norm2(r1)
  if ( th1 .ne. 0.0_wp ) then
    n1 = r1 / th1
    q1 = (/ cos(th1/2), n1(1)*sin(th1/2), &
                        n1(2)*sin(th1/2), &
                        n1(3)*sin(th1/2) /)
  else
    q1 = (/ 1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp /)
  end if
  th2 = norm2(r2)
  if ( th2 .ne. 0.0_wp ) then
    n2 = r2 / th2
    q2 = (/ cos(th2/2), n2(1)*sin(th2/2), &
                        n2(2)*sin(th2/2), &
                        n2(3)*sin(th2/2) /)
  else
    q2 = (/ 1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp /)
  end if

  !> Quaternion qs = q1 * q2 = cos(th/2) + sin(th/2) n · (i,j,k)
  qs = (/ &
       q1(1)*q2(1) - q1(2)*q2(2) - q1(3)*q2(3) - q1(4)*q2(4) , &
       q1(1)*q2(2) + q1(2)*q2(1) + q1(3)*q2(4) - q1(4)*q2(3) , &
       q1(1)*q2(3) + q1(3)*q2(1) + q1(4)*q2(2) - q1(2)*q2(4) , &
       q1(1)*q2(4) + q1(4)*q2(1) + q1(2)*q2(3) - q1(3)*q2(2)   &
       /)

  !> Rotation vector
  sin_t_2 = norm2(qs(2:4))
  if ( sin_t_2 .eq. 0.0_wp ) then
    n = (/ 1.0_wp, 0.0_wp, 0.0_wp /)
    cos_t_2 = 1.0_wp
    th = 0.0_wp
  else
    n = qs(2:4) / sin_t_2
    cos_t_2 = qs(1)
    th = 2.0_wp * atan2( sin_t_2, cos_t_2 )
  end if

  r = th * n


end subroutine rotation_vector_combination


! ----------------------------------------------------------------------

subroutine sort_vector_real( vec, nel, sor, ind )
  real(wp), intent(inout)               :: vec(:)
  integer , intent(in)                  :: nel
  real(wp), allocatable, intent(out)    :: sor(:)
  integer , allocatable, intent(out)    :: ind(:)

  real(wp)                              :: maxv
  integer                               :: i

  allocate(sor(nel)); sor = 0.0_wp
  allocate(ind(nel)); ind = 0

  maxv = maxval(vec)
  do i = 1, nel
    sor(i) = minval(vec, 1)
    ind(i) = minloc(vec, 1)
    vec(ind(i)) = maxv + 0.1_wp 
  end do

end subroutine sort_vector_real

! ----------------------------------------------------------------------

end module mod_math
