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
!!          Davide Montagnani
!!=========================================================================

module mod_spline

use mod_param, only: &
  wp, max_char_len, pi

use mod_handling, only: &
  error, warning, info, printout, new_file_unit, check_file_exists
  
implicit none

public :: t_spline , hermite_spline , deallocate_spline

private

character(len=*), parameter :: this_mod_name = 'mod_spline' 
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
subroutine hermite_spline ( mesh_file, spl, nelems , par_tension , par_bias      , &
                            type_span, rr_spl , nn_spl , &
                            ip_spl, ss_spl , &
                            spl_spl, &
                            leng, s_in , nor_in )
  type(t_spline)               , intent(inout) :: spl
  integer                      , intent(in)    :: nelems
  real(wp)                     , intent(in)    :: par_tension
  real(wp)                     , intent(in)    :: par_bias
  character(len=*)             , intent(in)    :: type_span
  character(len=*)             , intent(in)    :: mesh_file
  real(wp)                     , intent(out)   :: rr_spl(:,:)
  real(wp)                     , intent(out)   :: nn_spl(:,:)
  integer                      , intent(out)   :: ip_spl(:,:)
  real(wp)                     , intent(out)   :: ss_spl(:)
  real(wp)                     , intent(out)   :: spl_spl(:)
  real(wp)                     , intent(out)   :: leng
  real(wp)                     , intent(out)   ::   s_in(:)
  real(wp)                     , intent(out)   :: nor_in(:,:)
  character(len=*), parameter                  :: this_sub_name = 'hermite_spline'
  integer :: n_d , n_points , n_splines
  real(wp) , allocatable :: ll(:)  ! useless with this def of %t

  real(wp) , allocatable :: spl_s(:)

  ! spline reconstruction
  integer , parameter :: n_points1 = 100  ! now hardcoded (pass as a default input)
  real(wp) , allocatable :: tv1(:) , s(:) , s_spl(:)
  real(wp) , allocatable :: rr(:,:)

  integer :: i_r , i , j

  !> scale derivatives
  spl%d0 = spl%d0
  spl%d1 = spl%d1


  !> check input consistency
  if ( nelems .ne. size(rr_spl,1)-1 ) then
    write(*,*) ' error in hermite_spline : '
    write(*,*) '  nelems .ne. size(rr_spl,1)-1 . stop ' ; stop
  end if
  if ( nelems .ne. size(nn_spl,1)-1 ) then
    write(*,*) ' error in hermite_spline : '
    write(*,*) '  nelems .ne. size(nn_spl,1)-1 . stop ' ; stop
  end if

  !> ...
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
  deallocate(ll)

  allocate(spl%dp(n_points,n_d)) ; spl%dp = 0.0_wp

  !> following Paul Bourke description of splines
  !> http://paulbourke.net/miscellaneous/interpolation/

  !> === inner points ===
  do i = 2 , n_points - 1

    spl%dp(i,:) = 0.5_wp * ( 1.0_wp - par_tension ) * &
               ( ( spl%rr(i  ,:) -spl%rr(i-1,:) ) * ( 1.0_wp + par_bias ) + &
                 ( spl%rr(i+1,:) -spl%rr(i  ,:) ) * ( 1.0_wp - par_bias ) )

  end do

  !> === end points ===
  !-> first point
  if (     trim( spl%e_bc(1) ) .eq. 'derivative' ) then

    spl%dp(1,:) = spl%d0 / norm2(spl%d0) * norm2( spl%dp(2,:) )

  elseif ( trim( spl%e_bc(1) ) .eq. 'free'       ) then
    spl%dp(1,:) = 0.5_wp * ( - spl%dp(2,:) + 3.0_wp / spl%dt(1) * &
                            ( spl%rr(2,:) - spl%rr(1,:) ) ) ;
  else
    write(*,*) ' hermite_spline. Wrong b.c. '
    write(*,*) ' Only "derivative" implemented so far. Stop.' ; stop
  end if

  !-> last point
  if (     trim( spl%e_bc(2) ) .eq. 'derivative' ) then
    spl%dp(n_points,:) = spl%d1 / norm2(spl%d1) * norm2(spl%dp(n_points-1,:))
  elseif ( trim( spl%e_bc(2) ) .eq. 'free'       ) then
    spl%dp(n_points,:) = 0.5_wp * ( - spl%dp(n_points-1,:) + 3.0_wp / spl%dt(n_splines) * &
                                    ( spl%rr(n_points  ,:) - spl%rr(n_points-1,:) ) )
  else
    write(*,*) ' hermite_spline. Wrong b.c. '
    write(*,*) ' Only "derivative" implemented so far. Stop.' ; stop
  end if


  ! === tangent vector to the reference line ===
  nor_in = 0.0_wp
  do i = 1 , n_points
    nor_in(i,:) = spl%dp(i,:)
    nor_in(i,:) = nor_in(i,:) / norm2(nor_in(i,:))
  end do


  ! === spline (non-uniform param) ===
  allocate(rr( n_points1*n_splines + 1,n_d))
  allocate(tv1(n_points1))
  allocate( s( n_points1*n_splines + 1 ) )   ! overall curvilinear coord.
  allocate(spl_s(n_points)) ; spl_s = 0.0_wp ! overall curv. coord of the input points
  i_r = 0
  do i = 1 , n_splines

    if ( i .ne. n_splines ) then
      tv1 = (/ ( real(j,wp) , j = 0,n_points1-1 ) /) / real(n_points1,wp)
    else
      if ( allocated(tv1) ) deallocate(tv1)
      allocate(tv1(n_points1+1))
      tv1 = (/ ( real(j,wp) , j = 0,n_points1   ) /) / real(n_points1,wp)
    end if

    do j = 1 , size(tv1)
      i_r = i_r + 1  ! update index
      !> point coords
      rr(i_r,:) = hermite_p1(tv1(j)) * spl%rr(i  ,:) + &
                  hermite_p2(tv1(j)) * spl%rr(i+1,:) + &
                  hermite_d1(tv1(j)) * spl%dp(i  ,:) + &
                  hermite_d2(tv1(j)) * spl%dp(i+1,:)
      !> curvilinear coord.
      if ( i_r .ne. 1 ) then
        s(i_r) = s(i_r-1) + norm2(rr(i_r,:) - rr(i_r-1,:) )
      else
        s(i_r) = 0.0_wp  ! first point: s = 0
      end if
    end do

    spl_s(i+1) = s(i_r)

  end do

  !> s \in [ 0 , 1 ]
  leng = s(size(s))
  s     = s / s(size(s))

  s_in = spl_s           ! curvilinear coord of the input points
  spl_s = spl_s / spl_s(size(spl_s)) ! and its non-dimensionalisation


  deallocate(tv1)

  ! === spline ===
  !  allocate(rr_spl(nelems+1,3)) <- passed as an input
  allocate(s_spl(nelems+1)) ! curvilinear coord of the output points
  if ( trim(type_span) .eq. 'uniform' ) then
    s_spl = (/ ( real(j,wp) , j = 0, nelems ) /) / real(nelems,wp)
  elseif ( trim(type_span) .eq. 'cosine' ) then
    s_spl = 0.5_wp * ( 1.0_wp - cos( (/ (real(j,wp) , j = 0, nelems) /)*pi/ real(nelems,wp) ) )
  elseif ( trim(type_span) .eq. 'cosineOB' ) then
    s_spl = sin( 0.5_wp*(/(real(j,wp) , j = 0, nelems)/)*pi/ real(nelems,wp) )
  elseif ( trim(type_span) .eq. 'cosineIB' ) then
    s_spl = 1.0_wp - cos(0.5_wp*(/(real(j,wp) , j = 0, nelems)/)*pi/ real(nelems,wp) )
  elseif ( trim(type_span) .eq. 'equalarea' ) then
    s_spl = sqrt((/(real(j,wp) , j = 0, nelems)/)/real(nelems,wp)) 
  else 
    write(*,*) ' Mesh file   : ' , trim(mesh_file)
    write(*,*) ' type_span   : ' , trim(type_span)
    call error(this_sub_name, this_mod_name, 'Incorrect input: &
          & type_span must be equal to uniform, cosine, cosineIB, cosineOB.')
  end if

  !> first and last points
  rr_spl(       1,:) = spl%rr(       1,:)
  rr_spl(nelems+1,:) = spl%rr(n_points,:)

  nn_spl(       1,:) = rr(2,:) - rr(1,:)
  nn_spl(       1,:) = nn_spl(1,:) / norm2(nn_spl(1,:))
  nn_spl(nelems+1,:) = rr(size(rr,1),:) - rr(size(rr,1)-1,:)
  nn_spl(nelems+1,:) = nn_spl(nelems+1,:) / norm2(nn_spl(nelems+1,:))

  ip_spl(       1,:) = (/ 1 , 2 /)
  ip_spl(nelems+1,:) = (/ n_points-1, n_points /)
  ss_spl(       1)   = 0.0_wp
  ss_spl(nelems+1)   = 1.0_wp

  do i = 2 , nelems

    ! interp from the finest discretisation
    do j = 1 , size(s)
      if ( s(j) .gt. s_spl(i)  ) then

        if ( j .eq. 1 ) then
          write(*,*) ' error in hermite_spline, during the definition &
                      &of the points on the spline. Stop ' ; stop
        end if

        rr_spl(i,:) = ( s(j) - s_spl(i)   )/( s(j)-s(j-1) ) * rr(j-1,:) + &
                      ( s_spl(i) - s(j-1) )/( s(j)-s(j-1) ) * rr(j  ,:)

        nn_spl(i,:) = rr(j,:) - rr(j-1,:)
        nn_spl(i,:) = nn_spl(i,:) / norm2(nn_spl(i,:))

        exit

      end if
    end do

    ! find neighbouring input points (for input interpolation)
    do j = 1 , n_points
      if ( spl_s(j) .gt. s_spl(i) ) then

        ip_spl(i,:) = (/ j-1 , j /)
        ss_spl(i)   = ( s_spl(i) - spl_s(j-1) )/( spl_s(j)-spl_s(j-1) )

        exit

      end if
    end do

    ! spl_spl
    spl_spl(i) = ( s_spl(i) - spl_s(1) ) / ( spl_s(size(spl_s)) - spl_s(1) )

  end do


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



! !> build hermite_spline
! subroutine hermite_spline ( spl , nelems , der_factor , &
!                                            type_span , rr_spl , nn_spl , &
!                                                        ip_spl , ss_spl , &
!                                                        spl_spl  , &
!                                                        leng , s_in , nor_in )
!   type(t_spline)               , intent(inout) :: spl
!   integer                      , intent(in)    :: nelems
!   real(wp)                     , intent(in)    :: der_factor
!   character(max_char_len)      , intent(in)    :: type_span
!   real(wp)                     , intent(out)   :: rr_spl(:,:)
!   real(wp)                     , intent(out)   :: nn_spl(:,:)
!   integer                      , intent(out)   :: ip_spl(:,:)
!   real(wp)                     , intent(out)   :: ss_spl(:)
!   real(wp)                     , intent(out)   :: spl_spl(:)
!   real(wp)                     , intent(out)   :: leng
!   real(wp)                     , intent(out)   ::   s_in(:)
!   real(wp)                     , intent(out)   :: nor_in(:,:)
!
!   integer :: n_d , n_points , n_splines
!   real(wp) , allocatable :: ll(:)  ! useless with this def of %t
!   real(wp) , allocatable :: A(:,:) , b(:)
!
!   real(wp) , allocatable :: spl_s(:)
!
!   ! lapack
!   integer , allocatable :: ipiv(:)
!   integer :: info
!
!   ! spline reconstruction
!   integer , parameter :: n_points1 = 100  ! now hardcoded (pass as a default input)
!   real(wp) , allocatable :: tv1(:) , s(:) , s_spl(:)
!   real(wp) , allocatable :: rr(:,:)
!
!   integer :: i_r , i , j , n
!
!   !> scale derivatives
!   spl%d0 = spl%d0
!   spl%d1 = spl%d1
!
!
!   !> check input consistency
!   if ( nelems .ne. size(rr_spl,1)-1 ) then
!     write(*,*) ' error in hermite_spline : '
!     write(*,*) '  nelems .ne. size(rr_spl,1)-1 . stop ' ; stop
!   end if
!   if ( nelems .ne. size(nn_spl,1)-1 ) then
!     write(*,*) ' error in hermite_spline : '
!     write(*,*) '  nelems .ne. size(nn_spl,1)-1 . stop ' ; stop
!   end if
!
!   !> ...
!   if ( allocated(spl% t) ) deallocate(spl% t)
!   if ( allocated(spl%dt) ) deallocate(spl%dt)
!   if ( allocated(spl%dp) ) deallocate(spl%dp)
!
!
!   !> dimensions
!   n_d      = size(spl%rr,2)
!   n_points = size(spl%rr,1) ; n_splines = n_points - 1
!
!   ! === parameter vectors %t, %dt ===
!   ! todo: check the sensitivity of the spline w.r.t. the param array t
!   allocate( spl%t ( n_points  ) ) ; spl%t  = 0.0_wp
!   allocate( spl%dt( n_splines ) ) ; spl%dt = 0.0_wp
!   allocate( ll    ( n_splines ) ) ; ll     = 0.0_wp
!   spl%t( 1) = 0.0_wp
!   ll(    1) = norm2( spl%rr(2,:) - spl%rr(1,:) )
!   spl%t( 2) = ll(1)
!   spl%dt(1) = spl%t(2) - spl%t(1)
!   do i = 2 , n_splines
!     ll(    i  ) = norm2( spl%rr(i+1,:) - spl%rr(i,:) ) + ll(i-1)
!     spl%t( i+1) = ll(i)
!     spl%dt(i  ) = spl%t(i+1) - spl%t(i)
!   end do
!   deallocate(ll)
!
!
!   ! === Linear system to solve for the derivatives of the spline ===
!   !> one linear system per coorindate since the systems are not coupled
!   allocate(spl%dp(n_points,n_d)) ; spl%dp = 0.0_wp
!   allocate(A(n_points,n_points)) ; A = 0.0_wp
!   allocate(b(n_points))          ; b = 0.0_wp
!
!   do n = 1 , n_d
!
!     !> Initialise (reset) A,b
!     A = 0.0_wp ; b = 0.0_wp
!
!     !> inner points
!     do i = 2 , n_points - 1
!
!       A(i,i-1:i+1) = (/ spl%dt(i) ,  &
!                 2.0_wp*(spl%dt(i-1)+spl%dt(i)) , &
!                         spl%dt(i-1) /)
!       b(i) = 3.0_wp * spl%dt(i-1)/spl%dt(i  ) * &
!                     ( spl%rr(i+1,n) - spl%rr(i  ,n) ) + &
!              3.0_wp * spl%dt(i  )/spl%dt(i-1) * &
!                     ( spl%rr(i  ,n) - spl%rr(i-1,n) )
!
!     end do
!
!     !> end points: boundary conditions
!     if ( trim( spl%e_bc(1) ) .eq. 'derivative' ) then
!       A(1,1) = 1.0_wp
!       b(1  ) = spl%d0(n)
!     else
!       write(*,*) ' hermite_spline. Wrong b.c. '
!       write(*,*) ' Only "derivative" implemented so far. Stop.' ; stop
!     end if
!     if ( trim( spl%e_bc(2) ) .eq. 'derivative' ) then
!       A(n_points,n_points) = 1.0_wp
!       b(n_points         ) = spl%d1(n)
!     else
!       write(*,*) ' hermite_spline. Wrong b.c. '
!       write(*,*) ' Only "derivative" implemented so far. Stop.' ; stop
!     end if
!
!     ! ! check ---
!     ! do i = 1 , size(A,1)
!     !   write(*,*) A(i,:) , '           ' , b(i)
!     ! end do
!     ! write(*,*)
!
!     !> add a multiplication factor
!     A = A * der_factor
!
!     !> solve the linear system ( call to lapack dgesv )
!     spl%dp(:,n) = b   ! intent(inout) in dgesv
!     allocate(ipiv(n_points))
!     call dgesv( n_points, 1 , A , n_points , &
!              ipiv , spl%dp(:,n) , n_points , info )
!     deallocate(ipiv)
!
!   end do
!
!   ! === tangent vector to the reference line ===
!   nor_in = 0.0_wp
!   do i = 1 , n_points
!     nor_in(i,:) = spl%dp(i,:)
!     nor_in(i,:) = nor_in(i,:) / norm2(nor_in(i,:))
!   end do
!
!   ! === spline (non-uniform param) ===
!   allocate(rr( n_points1*n_splines + 1,n_d))
!   allocate(tv1(n_points1))
!   allocate( s( n_points1*n_splines + 1 ) )   ! overall curvilinear coord.
!   allocate(spl_s(n_points)) ; spl_s = 0.0_wp ! overall curv. coord of the input points
!   i_r = 0
!   do i = 1 , n_splines
!
!     if ( i .ne. n_splines ) then
!       tv1 = (/ ( dble(j) , j = 0,n_points1-1 ) /) / dble(n_points1)
!     else
!       if ( allocated(tv1) ) deallocate(tv1)
!       allocate(tv1(n_points1+1))
!       tv1 = (/ ( dble(j) , j = 0,n_points1   ) /) / dble(n_points1)
!     end if
!
!     do j = 1 , size(tv1)
!       i_r = i_r + 1  ! update index
!       !> point coords
!       rr(i_r,:) = hermite_p1(tv1(j)) * spl%rr(i  ,:) + &
!                   hermite_p2(tv1(j)) * spl%rr(i+1,:) + &
!                   hermite_d1(tv1(j)) * spl%dp(i  ,:) + &
!                   hermite_d2(tv1(j)) * spl%dp(i+1,:)
!       !> curvilinear coord.
!       if ( i_r .ne. 1 ) then
!         s(i_r) = s(i_r-1) + norm2(rr(i_r,:) - rr(i_r-1,:) )
!       else
!         s(i_r) = 0.0_wp  ! first point: s = 0
!       end if
!     end do
!
!     spl_s(i+1) = s(i_r)
!
!   end do
!
!   !> s \in [ 0 , 1 ]
!   leng = s(size(s))
!   s     = s / s(size(s))
!
!   s_in = spl_s           ! curvilinear coord of the input points
!   spl_s = spl_s / spl_s(size(spl_s)) ! and its non-dimensionalisation
!
! ! ! check ---
! ! do i = 1 , size(spl_s)
! !   write(*,*) spl_s(i)
! ! end do
! !
! ! ! check ---
! ! do i = 1 , size(rr,1)
! !   write(*,*) ' rr,  s : ' , rr(i,:) , '   ' , s(i)
! ! end do
!
!   deallocate(tv1)
!
!   ! === spline ===
!   !  allocate(rr_spl(nelems+1,3)) <- passed as an input
!   allocate(s_spl(nelems+1)) ! curvilinear coord of the output points
!   if ( trim(type_span) .eq. 'uniform' ) then
!     s_spl = (/ ( dble(j) , j = 0,nelems ) /) / dble(nelems)
!   else
!     write(*,*) ' error in hermite_spline. Only type_span = uniform '
!     write(*,*) ' implemented so far. Stop. ' ; stop
!   end if
!
!   !> first and last points
!   rr_spl(       1,:) = spl%rr(       1,:)
!   rr_spl(nelems+1,:) = spl%rr(n_points,:)
!
!   nn_spl(       1,:) = rr(2,:) - rr(1,:)
!   nn_spl(       1,:) = nn_spl(1,:) / norm2(nn_spl(1,:))
!   nn_spl(nelems+1,:) = rr(size(rr,1),:) - rr(size(rr,1)-1,:)
!   nn_spl(nelems+1,:) = nn_spl(nelems+1,:) / norm2(nn_spl(nelems+1,:))
!
!   ip_spl(       1,:) = (/ 1 , 2 /)
!   ip_spl(nelems+1,:) = (/ n_points-1, n_points /)
!   ss_spl(       1)   = 0.0_wp
!   ss_spl(nelems+1)   = 1.0_wp
!
!   do i = 2 , nelems
!
!     ! interp from the finest discretisation
!     do j = 1 , size(s)
!       if ( s(j) .gt. s_spl(i)  ) then
!
!         if ( j .eq. 1 ) then
!           write(*,*) ' error in hermite_spline, during the definition &
!                       &of the points on the spline. Stop ' ; stop
!         end if
!
!         rr_spl(i,:) = ( s(j) - s_spl(i)   )/( s(j)-s(j-1) ) * rr(j-1,:) + &
!                       ( s_spl(i) - s(j-1) )/( s(j)-s(j-1) ) * rr(j  ,:)
!
!         nn_spl(i,:) = rr(j,:) - rr(j-1,:)
!         nn_spl(i,:) = nn_spl(i,:) / norm2(nn_spl(i,:))
!
!         exit
!
!       end if
!     end do
!
!     ! find neighbouring input points (for input interpolation)
!     do j = 1 , n_points
!       if ( spl_s(j) .gt. s_spl(i) ) then
!
!         ip_spl(i,:) = (/ j-1 , j /)
!         ss_spl(i)   = ( s_spl(i) - spl_s(j-1) )/( spl_s(j)-spl_s(j-1) )
!
!         exit
!
!       end if
!     end do
!
!     ! spl_spl
!     spl_spl(i) = ( s_spl(i) - spl_s(1) ) / ( spl_s(size(spl_s)) - spl_s(1) )
!
!   end do
!
! ! ! check
! ! write(*,*) ' in hermite_spline : ip , ss '
! ! do i = 1 , size(ss_spl)
! !   write(*,*) ip_spl(i,:) , '     ' , ss_spl(i)
! ! end do
!
! end subroutine hermite_spline
