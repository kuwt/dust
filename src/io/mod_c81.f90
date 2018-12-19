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


! SUBROUTINEs:
!   read_c81_table
!   interp2d

!> Module to treat the c81 tables
module mod_c81

use mod_param, only: &
  wp, max_char_len

use mod_handling, only: &
  error, warning

implicit none

!-----------------------------------

! public and private

!-----------------------------------

type t_aero_2d_par
  real(wp) , allocatable :: par1(:)
  real(wp) , allocatable :: par2(:)
  real(wp) , allocatable :: cf(:,:)
end type t_aero_2d_par

!-----------------------------------

type t_aero_2d_tab

 !> value of the third parameter, Re
 real(wp) :: Re

 type(t_aero_2d_par) , allocatable :: coeff(:)

!!> cl(a,M,Re) 3d-data (extrapolation on Re or no dependence from Re)
!real(wp) , allocatable :: cl(:,:)
!type(t_aero_par) :: cl_par
!
!!> cd(a,M,Re) 3d-data
!real(wp) , allocatable :: cd(:,:)
!type(t_aero_par) :: cd_par
!
!!> cm(a,M,Re) 3d-data
!real(wp) , allocatable :: cm(:,:)
!type(t_aero_par) :: cm_par

end type t_aero_2d_tab

!-----------------------------------

!> Tables containing aerodynamic coefficients
!!
!! Tabulated data for lifting line elements 
type t_aero_tab

 !> airfoil file
 character(len=max_char_len) :: airfoil_file 

 !> airfoil data id
 integer :: id

 !> tables cl(a,M) for each Re number
 type(t_aero_2d_tab) , allocatable :: aero_coeff(:)

end type t_aero_tab

!-----------------------------------

character(len=*), parameter :: this_mod_name='mod_c81'

contains

!-----------------------------------
! Vahana input file style '../vahana_input/vahana4.inp'
! HP: the aerodynamic coefficients are given for the same values
! of al_i, M_i
subroutine read_c81_table ( filen , coeff )
  character(len=max_char_len) , intent(in) :: filen
  type(t_aero_tab) , intent(inout) :: coeff

  character(len=max_char_len) :: line , string
  integer :: cl1 , cl2 , cd1 , cd2 , cm1 , cm2 !  table_dim(2)
  integer :: cp1(3) , cp2(3)
  real(wp) :: dummy
  integer :: dummy_int , iblnk
  integer :: fid , i1 , iRe , nRe , i_c
  real(wp) :: Re

  integer , parameter :: n_coeff = 3   ! cl , cd , cm

  fid = 21
  open(unit=fid,file=trim(adjustl(filen))) 
  ! First tree lines containing unintelligilbe parameters
  read(fid,*) nRe , dummy_int , dummy_int ! 2 0 0
  read(fid,*) ! 0 1
  read(fid,*) ! 0.158 0.158

  allocate( coeff%aero_coeff(nRe) )

  do iRe = 1 , nRe

    ! For each Reynolds number (1 Reynolds number, 3 coefficients)
    read(fid,*) ! COMMENT#1
    read(fid,*) Re  , dummy
    read(fid,'(A)') line
      iblnk = index(line,' ') 
      string = trim(adjustl(line(iblnk:)))
      read(string,'(6I2)') cl2 , cl1 , cd2 , cd1 , cm2 , cm1
    cp1 = (/ cl1 , cd1 , cm1 /)
    cp2 = (/ cl2 , cd2 , cm2 /)

    ! Reynolds number
    coeff%aero_coeff(iRe)%Re = Re
    allocate(coeff%aero_coeff(iRe)%coeff(n_coeff))
    ! allocate coefficient structures
    do i_c = 1 , n_coeff
      allocate(coeff%aero_coeff(iRe)%coeff(i_c)%par1(cp1(i_c)))
      allocate(coeff%aero_coeff(iRe)%coeff(i_c)%par2(cp2(i_c)))
      allocate(coeff%aero_coeff(iRe)%coeff(i_c)%cf(cp1(i_c),cp2(i_c)))
    end do

    do i_c = 1 , n_coeff
      read(fid,*) coeff%aero_coeff(iRe)%coeff(i_c)%par2
      do i1 = 1 , cp1(i_c)
        read(fid,*) coeff%aero_coeff(iRe)%coeff(i_c)%par1(i1) , & 
                    coeff%aero_coeff(iRe)%coeff(i_c)%cf(i1,:)
      end do
    end do

  end do

  close(fid)



end subroutine read_c81_table

!-----------------------------------
!> Routine to compute aerodynamic coefficients,
!>  given the geometry ( airfoil_data, csi , airfoil_id) and
!>        the aerodynamic parameters( aero_par = (/ al , M , Re /)

! TODO: add some checks : if reyn1 .eq. reyn2 to avoid singularity
!       clear implementation
!       ...
subroutine interp_aero_coeff ( airfoil_data ,  &
                         csi , airfoil_id , aero_par , aero_coeff )
  type(t_aero_tab) , intent(in) :: airfoil_data(:)
  real(wp) , intent(in) :: csi
  integer  , intent(in) :: airfoil_id(2)
  real(wp) , intent(in) :: aero_par(3) ! (/al,M,Re/)
  real(wp) , allocatable , intent(out) :: aero_coeff(:) 

  real(wp) :: al , mach , reyn
  real(wp) :: cf1(3) , cf2(3)
  real(wp) , allocatable :: coeff1(:) , coeff2(:)
  real(wp) , allocatable :: coeff_airfoil(:,:) 
  integer :: nRe, i_a, id_a

  real(wp) :: reyn1 , reyn2
  integer  :: irey

  al   = aero_par(1)
  mach = aero_par(2)
  reyn = aero_par(3)

  cf1 = 0.0_wp
  cf2 = 0.0_wp

  ! Find the aerodynamic profiles to be interpolated in order to obtain
  ! the aerodynamic characteristics of the wing section
  !DEBUG
  do i_a = 1 , 2

   id_a = airfoil_id(i_a)
   nRe = size(airfoil_data(id_a)%aero_coeff)
  
   ! Some checks ----
   if ( nRe .eq. 1 ) then

    call interp2d_aero_coeff ( airfoil_data(id_a)%aero_coeff(1)%coeff , &
                                                 aero_par(1:2) , coeff1 )

    if ( .not. allocated(coeff_airfoil) ) then
      allocate(coeff_airfoil(2,size(coeff1)))
    end if

    coeff_airfoil(i_a,:) = coeff1

   else  
   ! Some checks ----
    if ( reyn .lt. airfoil_data(id_a)%aero_coeff(1)%Re ) then
     reyn1 = airfoil_data(id_a)%aero_coeff(1)%Re
     reyn2 = airfoil_data(id_a)%aero_coeff(2)%Re
     irey  = 1
    else if ( reyn .gt. airfoil_data(id_a)%aero_coeff(nRe)%Re ) then
     reyn1 = airfoil_data(id_a)%aero_coeff(nRe-1)%Re
     reyn2 = airfoil_data(id_a)%aero_coeff(nRe)%Re
     irey  = nRe-1
    else

     irey  = 1
     reyn1 = airfoil_data(id_a)%aero_coeff(irey)%Re
     do while ( ( reyn .ge. reyn1 ) .and. ( irey .lt. nRe ) )
      irey = irey + 1
      reyn1 = airfoil_data(id_a)%aero_coeff(irey)%Re 
     end do
     irey = irey - 1
     reyn1 = airfoil_data(id_a)%aero_coeff(irey)%Re 
     reyn2 = airfoil_data(id_a)%aero_coeff(irey+1)%Re 
    end if
    ! Some checks ----

    
    call interp2d_aero_coeff ( airfoil_data(id_a)%aero_coeff(irey)%coeff , &
                                                 aero_par(1:2) , coeff1 )
    call interp2d_aero_coeff ( airfoil_data(id_a)%aero_coeff(irey+1)%coeff , &
                                                 aero_par(1:2) , coeff2 )

    if ( .not. allocated(coeff_airfoil) ) then
      allocate(coeff_airfoil(2,size(coeff1)))
    end if

    coeff_airfoil(i_a,:) = ( coeff1 * ( reyn2 - reyn ) + coeff2 * ( reyn - reyn1 ) ) /  &
                           ( reyn2 - reyn1 ) 

   end if

  end do

  allocate( aero_coeff(size(coeff_airfoil,2)) )
  aero_coeff = coeff_airfoil(1,:) * ( 1-csi ) + coeff_airfoil(2,:) * csi


  ! deallocate
  deallocate(coeff_airfoil)


end subroutine interp_aero_coeff 

!-----------------------------------
! 
subroutine interp2d_aero_coeff ( aero_coeff , x , c )
 type(t_aero_2d_par) , intent(in) :: aero_coeff(:)
 real(wp), intent(in)  :: x(:)
 real(wp), allocatable , intent(out) :: c(:)

 integer  :: n1 , n2 , nc
 real(wp) :: csi1 , csi2
 real(wp) :: phi1 , phi2 , phi3 , phi4

 integer :: ic , i1 , i2

 character(len=*), parameter :: this_sub_name='interp2d_aero_coeff'

 nc = size(aero_coeff)

 allocate(c(nc)) ; c = 0.0_wp

 ! only 2d parameters are allowed (al,M)
 if ( size(x) .ne. 2 ) then 
   call error(this_sub_name, this_mod_name, 'Attempting to 2D interpolate&
   & data with more dimensions')
 end if

 ! Do it once ...
 ! Check dimensions ---------
 do ic = 1 , nc

   n1 = size(aero_coeff(ic)%par1)
   n2 = size(aero_coeff(ic)%par2)

!  ! DEBUG
!  write(*,*) ' n1 , n2 : ' , n1 , n2

   if ( size(aero_coeff(ic)%cf,1) .ne. n1 ) then
     write(*,*) ' Error in interp2d. '
     write(*,*) ' size(coeff(',ic,')%Mat,1) .ne. ', n1,'. STOP. ' ; stop
   end if
   if ( size(aero_coeff(ic)%cf,2) .ne. n2 ) then
     write(*,*) ' Error in interp2d. '
     write(*,*) ' size(coeff(',ic,')%Mat,2) .ne. n2. STOP. ' ; stop
   end if

   ! Check range of parameters the parameters are supposed to be defined in an
   ! increasing order
   if ( x(1) .lt. aero_coeff(ic)%par1(1) ) then
     write(*,*) ' Error in interp2d. x(1) .lt. minval(aero_coeff(ic)%par1). STOP' ; stop
   end if
   if ( x(1) .gt. aero_coeff(ic)%par1(n1) ) then
     write(*,*) ' Error in interp2d. x(1) .gt. maxval(aero_coeff(ic)%par1). STOP' ; stop
   end if
   if ( x(2) .lt. aero_coeff(ic)%par2(1) ) then
     write(*,*) ' Error in interp2d. x(2) .lt. minval(aero_coeff(ic)%par2). STOP' ; stop
   end if
   if ( x(2) .gt. aero_coeff(ic)%par2(n2) ) then
     write(*,*) ' Error in interp2d. x(2) .gt. maxval(aero_coeff(ic)%par2). STOP' ; stop
   end if
   ! Check dimensions ---------

   i1 = 1 
   do while ( (aero_coeff(ic)%par1(i1) .le. x(1)) ) 
     i1 = i1 + 1 
   end do
   i2 = 1 
   do while ( (aero_coeff(ic)%par2(i2) .le. x(2)) ) 
     i2 = i2 + 1 
   end do
  
   csi1 = 2.0_wp * ( x(1) - 0.5_wp*(aero_coeff(ic)%par1(i1-1)+aero_coeff(ic)%par1(i1)) ) / &
         (aero_coeff(ic)%par1(i1)-aero_coeff(ic)%par1(i1-1))

   csi2 = 2.0_wp * ( x(2) - 0.5_wp*(aero_coeff(ic)%par2(i2-1)+aero_coeff(ic)%par2(i2)) ) / &
         (aero_coeff(ic)%par2(i2)-aero_coeff(ic)%par2(i2-1))
  
   phi1 =   0.25_wp * ( 1 + csi1 ) * ( 1 + csi2 )
   phi2 =   0.25_wp * ( 1 - csi1 ) * ( 1 + csi2 )
   phi3 =   0.25_wp * ( 1 - csi1 ) * ( 1 - csi2 )
   phi4 =   0.25_wp * ( 1 + csi1 ) * ( 1 - csi2 )
  
   c(ic) = phi1 * aero_coeff(ic)%cf(i1  ,i2  ) + &
           phi2 * aero_coeff(ic)%cf(i1-1,i2  ) + &
           phi3 * aero_coeff(ic)%cf(i1-1,i2-1) + &
           phi4 * aero_coeff(ic)%cf(i1  ,i2-1)

 end do 


!  i1 = 1 
!  do while ( (x1(i1) .lt. x(1)) ) 
!    i1 = i1 + 1 
!  end do
!  i2 = 1 
!  do while ( (x2(i2) .lt. x(2)) ) 
!    i2 = i2 + 1 
!  end do
! 
!  csi1 = 2.0_wp * ( x(1) - 0.5_wp*(x1(i1-1)+x1(i1)) ) / (x1(i1)+x1(i1-1))
!  csi2 = 2.0_wp * ( x(2) - 0.5_wp*(x2(i2-1)+x2(i2)) ) / (x2(i2)+x2(i2-1))
! 
!  phi1 =   0.25_wp * ( 1 + csi1 ) * ( 1 + csi2 )
!  phi2 =   0.25_wp * ( 1 - csi1 ) * ( 1 + csi2 )
!  phi3 =   0.25_wp * ( 1 - csi1 ) * ( 1 - csi2 )
!  phi4 =   0.25_wp * ( 1 + csi1 ) * ( 1 - csi2 )
! 
!  allocate(c(nc)) ; c = 0.0_wp
!  do ic = 1 , nc
!    c(ic) = phi1 * coeff(ic)%Mat(i1  ,i2  ) + &
!            phi2 * coeff(ic)%Mat(i1-1,i2  ) + &
!            phi3 * coeff(ic)%Mat(i1-1,i2-1) + &
!            phi4 * coeff(ic)%Mat(i1  ,i2-1)
!  end do

end subroutine interp2d_aero_coeff

!-----------------------------------


end module mod_c81
