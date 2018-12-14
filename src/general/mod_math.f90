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


module mod_math

use mod_param, only: &
  wp

implicit none

public :: cross , linear_interp , compute_qr

private

interface linear_interp
  module procedure linear_interp_vector, &
                   linear_interp_array
end interface linear_interp

contains

! ----------------------------------------------------------------------

function cross(a, b)
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

  nt = size(t_vec)
  ! Check dimensions -----
  if ( size(val_vec) .ne. nt ) then
    write(*,*) ' ERROR: wrong sizes: size(val_vec) .ne. size(t_vec).'
    stop
  end if

  ! Check if t \in [ minval(t_vec) , maxval(t_vec) ]
  if ( t .lt. minval(t_vec) ) then
    write(*,*) ' ERROR: t .lt. t_vec ' 
  end if
  if ( t .gt. maxval(t_vec) ) then
    write(*,*) ' ERROR: t .gt. t_vec ' 
    write(*,*) ' t             = ' , t 
    write(*,*) ' maxval(t_vec) = ' , maxval(t_vec)
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

  nt = size(t_vec)
  allocate(val(size(val_arr,1)))

! check ----
! write(*,*) ' size(val_arr) : ' , size(val_arr,1) , size(val_arr,2)
! write(*,*) ' nt            : ' , nt
! write(*,*) ' tv : ' 
! check ----

  ! Check dimensions -----
  if ( size(val_arr,2) .ne. nt ) then
    write(*,*) ' ERROR: wrong sizes: size(val_vec,2) .ne. size(t_vec).'
    stop
  end if

  ! Check if t \in [ minval(t_vec) , maxval(t_vec) ]
  if ( t .lt. minval(t_vec) ) then
    write(*,*) ' ERROR: t .lt. t_vec ' 
  end if
  if ( t .gt. maxval(t_vec) ) then
    write(*,*) ' ERROR: t .gt. t_vec '
    write(*,*) ' t             = ' , t 
    write(*,*) ' maxval(t_vec) = ' , maxval(t_vec)
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
 
  integer :: m , n , i_m , i_n
 
  ! lapack dgeqrf routine ----
  real(wp) , allocatable :: tau(:) , work(:)
  integer :: lwork , info 
 
  ! tmp matrices to get Q matrix -----
  real(wp) , allocatable :: H(:,:) , v(:,:) , eye(:,:)
 
 
  ! input check and warnings
  if ( allocated(Q) ) then
    write(*,*) ' Warning. In compute_qr(), Q was already allocated.'
    write(*,*) ' Deallocated and re-allocated. '
    deallocate(Q)
  end if
  if ( allocated(R) ) then
    write(*,*) ' Warning. In compute_qr(), Q was already allocated.'
    write(*,*) ' Deallocated and re-allocated. '
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
! ! debug ----
! write(*,*) ' In compute_qr(), A : ' 
! do i_m = 1 , size(A,1) ; write(*,*) '  ' , A(i_m,:) ; end do
! ! debug ----

  call dgeqrf( m , n , A , m , tau , work , lwork , info )

! ! debug ----
! write(*,*) ' In compute_qr(), after calling dgeqrf(). A : ' 
! do i_m = 1 , size(A,1) ; write(*,*) '  ' , A(i_m,:) ; end do
! ! debug ----
 
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
