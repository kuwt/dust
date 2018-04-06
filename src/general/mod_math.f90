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

public :: cross , linear_interp

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

  write(*,*) ' size(val_arr) : ' , size(val_arr,1) , size(val_arr,2)
  write(*,*) ' nt            : ' , nt
  write(*,*) ' tv : ' 
  do it = 1 , nt
    write(*,*) t_vec(it) , val_arr(:,it)
  end do

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


end module mod_math
