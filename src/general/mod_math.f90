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

!> Module containing some mathematical utilities
module mod_math

use mod_param, only: &
  wp

use mod_handling, only: &
  error, warning

implicit none

public :: cross , linear_interp , compute_qr

private

interface linear_interp
  module procedure linear_interp_vector, &
                   linear_interp_array
end interface linear_interp

character(len=*), parameter :: this_mod_name='mod_math'

contains

! ----------------------------------------------------------------------

pure function cross(a, b)
  real(wp), dimension(3) :: cross
  real(wp), dimension(3), intent(in) :: a, b

  cross(1) = a(2) * b(3) - a(3) * b(2)
  cross(2) = a(3) * b(1) - a(1) * b(3)
  cross(3) = a(1) * b(2) - a(2) * b(1)

end function cross

! ----------------------------------------------------------------------

subroutine linear_interp_vector( val_vec , t_vec , t , val ) 
  real(wp) , intent(in) :: val_vec(:)
  real(wp) , intent(in) :: t_vec(:)
  real(wp) , intent(in) :: t
  real(wp) , intent(out) :: val

  integer :: it , nt
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

! ----------------------------------------------------------------------


end module mod_math
